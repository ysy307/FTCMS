submodule(io_input_basic) input_basic_base
    implicit none
contains

    module subroutine initialize_type_input_basic(self)
        !< Load the input parameters from the JSON file
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file) :: json

        call json%initialize()

        call json%load(filename=self%file_name)
        call json%print_error_message(output_unit)

        call read_simulation_settings(self, json)
        call read_analysis_controls(self, json)
        call read_geometry_settings(self, json)
        call read_materials(self, json)
        call read_solver_settings(self, json)
        call read_discretization_settings(self, json)

        call json%destroy()
        call json%print_error_message(output_unit)

    end subroutine initialize_type_input_basic

    module subroutine display_input_basic(self)
        implicit none
        class(type_input_basic), intent(in) :: self

        integer(int32) :: ierr, myrank
        integer(int32) :: i

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
        if (myrank == 0) then
            write (*, '(A)') "=== Simulation Settings ==="
            call self%simulation_settings%display()
            write (*, '(A)') "=== Analysis Controls ==="
            call self%analysis_controls%display()
            write (*, '(A)') "=== Geometry Settings ==="
            call self%geometry_settings%display()
            write (*, '(A)') "=== Material Settings ==="
            do i = 1, self%num_materials
                call self%materials(i)%display()
            end do
            write (*, '(A)') "=== Solver Settings ==="
            call self%solver_settings%display()
            write (*, '(A)') "=== Discretization Settings ==="
            call self%discretization_settings%display(self%analysis_controls%is_active)
        end if
    end subroutine display_input_basic

end submodule input_basic_base
