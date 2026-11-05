module io_input_geometry
    use, intrinsic :: iso_fortran_env
    use :: core_parallel_mpi
    use :: stdlib_strings, only:to_string, strip, ends_with
    use :: stdlib_logger
    use :: json_module, only:json_file
    use :: module_core
    use :: domain_mesh_plex, only:type_mesh_plex
    use :: io_input_base, only:get_json_value, abst_input
    implicit none
    private

    public :: type_input_geometry

    character(*), parameter :: thermal = "thermal"
    character(*), parameter :: calculate_thermal = "calculate_thermal"
    character(*), parameter :: hydraulic = "hydraulic"
    character(*), parameter :: calculate_hydraulic = "calculate_hydraulic"
    character(*), parameter :: mechanical = "mechanical"
    character(*), parameter :: calculate_mechanical = "calculate_mechanical"

    type :: type_input_geometry
        class(abst_input), pointer :: parent => null()
        !> The mesh itself, held as a PETSc DM. Nothing is copied out of it.
        type(type_mesh_plex) :: mesh
        character(:), allocatable :: point_data_names(:)
    contains
        procedure, pass(self), public :: initialize => initialize_type_input_geometry
        procedure, pass(self), private :: collect_fields_from_conditions
    end type type_input_geometry

    interface
        module subroutine initialize_type_input_geometry(self)
            implicit none
            class(type_input_geometry), intent(inout) :: self
        end subroutine initialize_type_input_geometry

        module subroutine collect_fields_from_conditions(self, field_list)
            implicit none
            class(type_input_geometry), intent(inout) :: self
            character(len=256), allocatable, intent(out) :: field_list(:)
        end subroutine collect_fields_from_conditions

        ! module subroutine display_input_geometry(self)
        !     implicit none
        !     class(type_input_geometry), intent(in) :: self
        ! end subroutine display_input_geometry
    end interface

end module io_input_geometry
