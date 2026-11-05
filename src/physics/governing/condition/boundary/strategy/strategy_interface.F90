module boundary_strategy
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: boundary_data_provider
    use :: boundary_data_provider_factory
    implicit none
    private

    public :: abst_bc
    public :: type_bc_result
    public :: type_bc_dirichlet
    public :: type_bc_flux
    public :: type_bc_robin
    public :: type_bc_atmospheric
    public :: type_bc_radiation
    public :: type_bc_convective
    public :: type_bc_seepage

    type :: type_bc_result
        logical :: is_dirichlet
        real(real64) :: prescribed_value
        real(real64) :: flux_value
        real(real64) :: flux_derivative
    contains
        procedure, pass(self) :: initialize => initialize_type_bc_result
    end type type_bc_result

    type, abstract :: abst_bc
        type(type_constant_id) :: physics_type = type_constant_id("", "", -1)
        type(type_constant_id) :: bc_kind = type_constant_id("", "", -1)
        class(abst_bc_data), allocatable, private :: provider
    contains
        procedure, public, pass(self) :: initialize => initialize_abst_bc
        procedure, public, pass(self) :: destroy => destroy_abst_bc
        procedure, public, pass(self) :: update_provider_data => update_provider_data_abst_bc
        procedure, public, pass(self) :: get_data => get_data_abst_bc
        procedure(abst_evaluate), public, pass(self), deferred :: evaluate
    end type abst_bc

    abstract interface
        subroutine abst_evaluate(self, current_time, u_curr, result)
            import :: abst_bc, real64, type_bc_result
            implicit none
            class(abst_bc), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine abst_evaluate
    end interface

    type, extends(abst_bc) :: type_bc_dirichlet
    contains
        procedure, public, pass(self) :: evaluate => evaluate_dirichlet_bc
    end type type_bc_dirichlet

    type, extends(abst_bc) :: type_bc_flux
    contains
        procedure, public, pass(self) :: evaluate => evaluate_flux_bc
    end type type_bc_flux

    type, extends(abst_bc) :: type_bc_robin
    contains
        procedure, public, pass(self) :: evaluate => evaluate_robin_bc
    end type type_bc_robin

    type, extends(abst_bc) :: type_bc_atmospheric
    contains
        procedure, public, pass(self) :: evaluate => evaluate_atmospheric_bc
    end type type_bc_atmospheric

    type, extends(abst_bc) :: type_bc_radiation
    contains
        procedure, public, pass(self) :: evaluate => evaluate_radiation_bc
    end type type_bc_radiation

    type, extends(abst_bc) :: type_bc_convective
    contains
        procedure, public, pass(self) :: evaluate => evaluate_convective_bc
    end type type_bc_convective

    type, extends(abst_bc) :: type_bc_seepage
    contains
        procedure, public, pass(self) :: evaluate => evaluate_seepage_bc
    end type type_bc_seepage

    interface
        module pure elemental subroutine initialize_type_bc_result(self)
            implicit none
            class(type_bc_result), intent(inout) :: self

        end subroutine initialize_type_bc_result

        module subroutine initialize_abst_bc(self, physics_type, config)
            implicit none
            class(abst_bc), intent(inout), target :: self
            type(type_constant_id), intent(in) :: physics_type
            type(type_config_bc), intent(in) :: config
        end subroutine initialize_abst_bc

        module subroutine destroy_abst_bc(self)
            implicit none
            class(abst_bc), intent(inout) :: self
        end subroutine destroy_abst_bc

        module subroutine update_provider_data_abst_bc(self, new_values)
            implicit none
            class(abst_bc), intent(inout) :: self
            real(real64), intent(in) :: new_values(:)
        end subroutine update_provider_data_abst_bc

        module subroutine get_data_abst_bc(self, current_time, values)
            implicit none
            class(abst_bc), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(inout) :: values(BC_DATA_SIZE)
        end subroutine get_data_abst_bc

        module subroutine evaluate_dirichlet_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_dirichlet), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_dirichlet_bc

        module subroutine evaluate_flux_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_flux), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_flux_bc

        module subroutine evaluate_robin_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_robin), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_robin_bc

        module subroutine evaluate_atmospheric_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_atmospheric), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_atmospheric_bc

        module subroutine evaluate_radiation_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_radiation), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_radiation_bc

        module subroutine evaluate_convective_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_convective), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_convective_bc

        module subroutine evaluate_seepage_bc(self, current_time, u_curr, result)
            implicit none
            class(type_bc_seepage), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: u_curr
            type(type_bc_result), intent(inout) :: result
        end subroutine evaluate_seepage_bc
    end interface

end module boundary_strategy
