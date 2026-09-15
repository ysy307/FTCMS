#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscviewerhdf5.h>

!> Field output through PETSc's HDF5 viewer.
!>
!> The geometry written is the DM itself, so this writer has no node ordering of
!> its own to keep in step with the mesh, and the viewer is collective, so one
!> file is produced however many ranks are running.
!>
!> PETSc lays the file out as /geometry, /viz/topology and /fields. The XDMF
!> descriptor that a viewer needs alongside it is not written here: its data
!> items carry the dataset shapes, and PETSc exposes no way to read those back
!> from Fortran. Generate it with PETSc's own lib/petsc/bin/petsc_gen_xdmf.py,
!> which reads the shapes from the file.
submodule(io_output_overall) output_overall_hdf5
    use :: petscdmplex
    use :: stdlib_logger, only:global_logger
    implicit none

contains

    module subroutine initialize_type_output_overall_hdf5(self, dir_output, config, mesh)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        character(*), intent(in) :: dir_output
        type(type_config_overall), intent(in) :: config
        type(type_mesh_plex), intent(inout), optional, target :: mesh

        self%directory = dir_output
        self%file_format = config%file_format
        self%space_dimension = 2
        if (present(mesh)) self%space_dimension = mesh%dimension
        self%num_points = config%num_points
        self%num_cells = config%num_cells

        if (present(mesh)) self%mesh => mesh
        if (.not. associated(self%mesh)) then
            call global_logger%log_error("HDF5 output needs the mesh; none was handed to it")
            error stop "HDF5 output without a mesh"
        end if
    end subroutine initialize_type_output_overall_hdf5

    module subroutine write_fields_hdf5(self, file_counts, temperature, water_content, &
                                        ice_content, vapor_content, pressure, water_flux, time)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        integer(int32), intent(in) :: file_counts
        real(real64), intent(in), optional :: temperature(:)
        real(real64), intent(in), optional :: water_content(:)
        real(real64), intent(in), optional :: ice_content(:)
        real(real64), intent(in), optional :: vapor_content(:)
        real(real64), intent(in), optional :: pressure(:)
        type(type_coordinate_array_dp), intent(in), optional :: water_flux
        real(real64), intent(in), optional :: time

        PetscErrorCode :: ierr
        PetscViewer :: viewer
        DM :: dm
        character(:), allocatable :: file_name
        integer(int32) :: num_written

        if (.not. associated(self%mesh)) return
        if (present(time)) self%current_time = time

        ! The output clone, not the solver's DM: fields written through it are
        ! placed against the geometry PETSc emits.
        call self%mesh%get_output_dm(dm)

        if (self%file_format == FILE_FORMATS%PETSC_VTU) then
            file_name = self%directory//"Out_"//pad_counter(file_counts)//".vtu"
            call PetscViewerVTKOpen(PETSC_COMM_WORLD, trim(file_name), FILE_MODE_WRITE, viewer, ierr)
            if (ierr /= 0) then
                call global_logger%log_error("Could not open "//file_name//" for VTU output")
                return
            end if
            ! The .vtu extension already selects the format; the viewer holds
            ! the fields and writes the file when it is destroyed.
        else
            file_name = self%directory//"Out_"//pad_counter(file_counts)//".h5"
            call PetscViewerHDF5Open(PETSC_COMM_WORLD, trim(file_name), FILE_MODE_WRITE, viewer, ierr)
            if (ierr /= 0) then
                call global_logger%log_error("Could not open "//file_name//" for HDF5 output")
                return
            end if
            call PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_XDMF, ierr)
            ! The mesh goes in first: the fields are described relative to it.
            call DMView(dm, viewer, ierr)
        end if

        num_written = 0
        call write_nodal_field(self, dm, viewer, "temperature", temperature, num_written)
        call write_nodal_field(self, dm, viewer, "water_content", water_content, num_written)
        call write_nodal_field(self, dm, viewer, "ice_content", ice_content, num_written)
        call write_nodal_field(self, dm, viewer, "vapor_content", vapor_content, num_written)
        call write_nodal_field(self, dm, viewer, "pressure", pressure, num_written)

        if (self%file_format /= FILE_FORMATS%PETSC_VTU) call PetscViewerPopFormat(viewer, ierr)
        ! A VTK viewer holds the fields until it is destroyed, then writes the
        ! file in one pass.
        call PetscViewerDestroy(viewer, ierr)
        if (ierr /= 0) call global_logger%log_error("Writing "//file_name//" failed")

        call record_timestep(self, file_counts, file_name)

    end subroutine write_fields_hdf5

    !> One nodal field, in a Vec created from the output DM so PETSc can place
    !> it against the geometry it writes.
    subroutine write_nodal_field(self, dm, viewer, name, values, num_written)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        DM, intent(in) :: dm
        PetscViewer, intent(in) :: viewer
        character(*), intent(in) :: name
        real(real64), intent(in), optional :: values(:)
        integer(int32), intent(inout) :: num_written

        PetscErrorCode :: ierr
        Vec :: field, local
        PetscScalar, pointer :: raw(:)
        integer(int32), allocatable :: vertex_node(:)
        integer(int32) :: i, n

        if (.not. present(values)) return
        if (size(values) <= 0) return

        call DMCreateLocalVector(dm, local, ierr)
        if (ierr /= 0) return

        ! The output vector carries one value per vertex, in the DM's order;
        ! our arrays are indexed by our own node numbering, so the values are
        ! gathered rather than copied.
        call self%mesh%get_vertex_nodes(vertex_node, n)
        call VecGetArray(local, raw, ierr)
        if (ierr == 0) then
            n = min(n, size(raw))
            do i = 1, n
                if (vertex_node(i) >= 1 .and. vertex_node(i) <= size(values)) then
                    raw(i) = values(vertex_node(i))
                else
                    raw(i) = 0.0d0
                end if
            end do
            call VecRestoreArray(local, raw, ierr)
        end if

        call DMCreateGlobalVector(dm, field, ierr)
        if (ierr == 0) call PetscObjectSetName(PetscObjectCast(field), name, ierr)
        ! Local to global drops the ghost copies, keeping one value per node.
        if (ierr == 0) call DMLocalToGlobalBegin(dm, local, INSERT_VALUES, field, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(dm, local, INSERT_VALUES, field, ierr)
        if (ierr == 0) call VecView(field, viewer, ierr)
        if (ierr /= 0) call global_logger%log_error("Could not write field "//name)
        if (ierr == 0) then
            num_written = num_written + 1
            call remember_field(self, name)
        end if

        call VecDestroy(field, ierr)
        call VecDestroy(local, ierr)
    end subroutine write_nodal_field

    !> Append this step to the ParaView collection.
    !>
    !> ParaView reads a numbered series on its own, but only as frame indices.
    !> The collection carries the physical time of each frame, which is what
    !> makes the animation show the simulated time.
    subroutine record_timestep(self, file_counts, data_file)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        integer(int32), intent(in) :: file_counts
        character(*), intent(in) :: data_file

        integer(int32) :: unit_pvd, ios, rank, mpi_ierr, step
        character(:), allocatable :: collection, entry_file

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpi_ierr)
        if (rank /= 0) return


        self%num_recorded = max(self%num_recorded, file_counts + 1)
        if (allocated(self%recorded_time)) then
            if (size(self%recorded_time) < self%num_recorded) call grow_times(self)
        else
            allocate (self%recorded_time(64))
            self%recorded_time = 0.0d0
        end if
        self%recorded_time(file_counts + 1) = self%current_time

        if (self%file_format /= FILE_FORMATS%PETSC_VTU) then
            call write_xdmf_series(self)
            return
        end if

        collection = self%directory//"Out.pvd"
        open (newunit=unit_pvd, file=collection, status="replace", action="write", iostat=ios)
        if (ios /= 0) return

        write (unit_pvd, '(A)') '<?xml version="1.0"?>'
        write (unit_pvd, '(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
        write (unit_pvd, '(A)') '  <Collection>'
        do step = 1, self%num_recorded
            entry_file = "Out_"//pad_counter(step - 1)//".vtu"
            write (unit_pvd, '(A,ES16.8,A)') '    <DataSet timestep="', self%recorded_time(step), &
                '" group="" part="0" file="'//entry_file//'"/>'
        end do
        write (unit_pvd, '(A)') '  </Collection>'
        write (unit_pvd, '(A)') '</VTKFile>'
        close (unit_pvd)
    end subroutine record_timestep

    !> Note that a field reached the file, so the descriptor names only what is
    !> actually there.
    subroutine remember_field(self, name)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        character(*), intent(in) :: name

        integer(int32) :: i

        if (.not. allocated(self%written_fields)) then
            allocate (self%written_fields(8))
            self%written_fields = ""
            self%num_written_fields = 0
        end if
        do i = 1, self%num_written_fields
            if (trim(self%written_fields(i)) == name) return
        end do
        if (self%num_written_fields >= size(self%written_fields)) return
        self%num_written_fields = self%num_written_fields + 1
        self%written_fields(self%num_written_fields) = name
    end subroutine remember_field

    !> The XDMF descriptor for the HDF5 series.
    !>
    !> HDF5 alone is not readable by a viewer: the descriptor is what says which
    !> dataset is the geometry, which is the topology and which are the fields,
    !> and what shape each has. The shapes come from the mesh rather than from
    !> the file, because PETSc exposes no way to read a dataset's shape back
    !> from Fortran.
    subroutine write_xdmf_series(self)
        implicit none
        class(type_output_overall_hdf5), intent(in) :: self

        integer(int32) :: unit_xdmf, ios, step
        integer(int32) :: num_interior_cells, num_corners
        character(:), allocatable :: topology_name, descriptor, data_file

        if (.not. associated(self%mesh)) return
        call self%mesh%get_visualisation_shape(num_interior_cells, num_corners, topology_name)
        if (num_corners <= 0) return

        descriptor = self%directory//"Out.xmf"
        open (newunit=unit_xdmf, file=descriptor, status="replace", action="write", iostat=ios)
        if (ios /= 0) return

        write (unit_xdmf, '(A)') '<?xml version="1.0" ?>'
        write (unit_xdmf, '(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        write (unit_xdmf, '(A)') '<Xdmf Version="2.0">'
        write (unit_xdmf, '(A)') '  <Domain>'
        write (unit_xdmf, '(A)') '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'

        do step = 1, self%num_recorded
            data_file = "Out_"//pad_counter(step - 1)//".h5"
            write (unit_xdmf, '(A)') '      <Grid Name="domain" GridType="Uniform">'
            write (unit_xdmf, '(A,ES16.8,A)') '        <Time Value="', self%recorded_time(step), '"/>'

            write (unit_xdmf, '(A,A,A,I0,A)') '        <Topology TopologyType="', topology_name, &
                '" NumberOfElements="', num_interior_cells, '">'
            write (unit_xdmf, '(A,I0,1X,I0,A)') '          <DataItem Dimensions="', &
                num_interior_cells, num_corners, '" NumberType="Int" Format="HDF">'
            write (unit_xdmf, '(A)') '            '//data_file//':/viz/topology/cells'
            write (unit_xdmf, '(A)') '          </DataItem>'
            write (unit_xdmf, '(A)') '        </Topology>'

            write (unit_xdmf, '(A)') '        <Geometry GeometryType="'//geometry_kind(self%space_dimension)//'">'
            write (unit_xdmf, '(A,I0,1X,I0,A)') '          <DataItem Dimensions="', &
                self%num_points, self%space_dimension, '" Format="HDF">'
            write (unit_xdmf, '(A)') '            '//data_file//':/geometry/vertices'
            write (unit_xdmf, '(A)') '          </DataItem>'
            write (unit_xdmf, '(A)') '        </Geometry>'

            call write_attribute(self, unit_xdmf, data_file, "temperature")
            call write_attribute(self, unit_xdmf, data_file, "water_content")
            call write_attribute(self, unit_xdmf, data_file, "ice_content")
            call write_attribute(self, unit_xdmf, data_file, "vapor_content")
            call write_attribute(self, unit_xdmf, data_file, "pressure")

            write (unit_xdmf, '(A)') '      </Grid>'
        end do

        write (unit_xdmf, '(A)') '    </Grid>'
        write (unit_xdmf, '(A)') '  </Domain>'
        write (unit_xdmf, '(A)') '</Xdmf>'
        close (unit_xdmf)
    end subroutine write_xdmf_series

    !> One nodal attribute, named as it was written into the file.
    subroutine write_attribute(self, unit_xdmf, data_file, name)
        implicit none
        class(type_output_overall_hdf5), intent(in) :: self
        integer(int32), intent(in) :: unit_xdmf
        character(*), intent(in) :: data_file
        character(*), intent(in) :: name

        if (.not. was_written(self, name)) return

        write (unit_xdmf, '(A)') '        <Attribute Name="'//name//'" AttributeType="Scalar" Center="Node">'
        write (unit_xdmf, '(A,I0,A)') '          <DataItem Dimensions="', self%num_points, '" Format="HDF">'
        write (unit_xdmf, '(A)') '            '//data_file//':/fields/'//name
        write (unit_xdmf, '(A)') '          </DataItem>'
        write (unit_xdmf, '(A)') '        </Attribute>'
    end subroutine write_attribute

    !> Whether a field reached the file. A run that does not solve a physics
    !> never writes its field, and naming it in the descriptor would break the
    !> whole series for a viewer.
    pure logical function was_written(self, name)
        implicit none
        class(type_output_overall_hdf5), intent(in) :: self
        character(*), intent(in) :: name

        integer(int32) :: i

        was_written = .false.
        if (.not. allocated(self%written_fields)) return
        do i = 1, self%num_written_fields
            if (trim(self%written_fields(i)) == name) then
                was_written = .true.
                return
            end if
        end do
    end function was_written


    pure function geometry_kind(dimension) result(kind_name)
        implicit none
        integer(int32), intent(in) :: dimension
        character(:), allocatable :: kind_name

        if (dimension >= 3) then
            kind_name = "XYZ"
        else
            kind_name = "XY"
        end if
    end function geometry_kind

    subroutine grow_times(self)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        real(real64), allocatable :: grown(:)
        integer(int32) :: capacity

        capacity = size(self%recorded_time)
        do while (capacity < self%num_recorded)
            capacity = 2 * capacity
        end do
        allocate (grown(capacity))
        grown = 0.0d0
        grown(1:size(self%recorded_time)) = self%recorded_time
        call move_alloc(grown, self%recorded_time)
    end subroutine grow_times

    module subroutine write_cell_hdf5(self, file_name, variable_name, variable)
        implicit none
        class(type_output_overall_hdf5), intent(inout) :: self
        character(*), intent(in) :: file_name
        character(*), intent(in) :: variable_name
        integer(int32), intent(in) :: variable(:)

        ! Cell data is not written yet: the quantities this code produces are
        ! nodal, and averaging them onto cells would invent a value.
    end subroutine write_cell_hdf5


    !> Zero-padded file counter, matching the other writers.
    pure function pad_counter(file_counts) result(text)
        implicit none
        integer(int32), intent(in) :: file_counts
        character(len=3) :: text

        write (text, '(I3.3)') file_counts
    end function pad_counter

end submodule output_overall_hdf5
