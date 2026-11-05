module materials_thermal_conductivity
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core, only:type_state, type_config_constitutive, allocate_array, type_coordinate_dp, type_state_thc
    use :: constitutive_constants, only:TtoK => celsius_to_kelvin
    use :: materials_base, only:abst_material
    implicit none
    private

    public :: holder_thcs
    public :: abst_thc
    public :: type_thc_1phase
    public :: type_thc_2phase
    public :: type_thc_3phase
    public :: type_thc_4phase
    public :: type_state_thc

    type :: holder_thcs
        class(abst_thc), allocatable :: p
    contains
        procedure, pass(self) :: initialize => initialize_holder_thcs
    end type holder_thcs

    interface
        module subroutine initialize_holder_thcs(self, material_id, constitutive_info, water, ice)
            implicit none
            class(holder_thcs), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_config_constitutive), intent(in) :: constitutive_info
            type(type_iapws97), intent(in), target :: water
            type(type_iapws06), intent(in), target :: ice

        end subroutine initialize_holder_thcs
    end interface


    type, extends(abst_material), abstract :: abst_thc
        real(real64), allocatable :: dispersivity(:)
        logical :: use_dispersivity = .false.
        real(real64), allocatable :: params(:)
    contains
        procedure, pass(self), public :: initialize => initialize_abst_thc
        procedure, pass(self) :: calc_lambda_dispersivity => calc_lambda_dispersivity_abst_thc
        procedure(abst_calc_thc_gp), pass(self), public, deferred :: calc_lambda_0_gp
        procedure(abst_calc_thc_dispersivity_gp), pass(self), public, deferred :: calc_lambda_dispersivity_gp
        generic, public :: calc => calc_lambda_0_gp, calc_lambda_dispersivity_gp
    end type abst_thc

    abstract interface
        !> Bulk thermal conductivity \(\lambda\) of a state, W/(m K), and
        !> optionally its tangent in \((p_w,T)\) through the phase contents.
        subroutine abst_calc_thc_gp(self, state, lambda, tangent)
            import :: abst_thc, type_state, real64
            implicit none
            class(abst_thc), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: lambda
            !> \((\partial\lambda/\partial p_w,\partial\lambda/\partial T)\); overwritten on exit when present
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine abst_calc_thc_gp

        subroutine abst_calc_thc_dispersivity_gp(self, state, lambda)
            import :: abst_thc, type_state, type_state_thc
            implicit none
            class(abst_thc), intent(in) :: self
            type(type_state), intent(in) :: state
            type(type_state_thc), intent(inout) :: lambda

        end subroutine abst_calc_thc_dispersivity_gp
    end interface

    interface
        module subroutine initialize_abst_thc(self, material_id, constitutive_info, water, ice)
            implicit none
            class(abst_thc), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_config_constitutive), intent(in) :: constitutive_info
            type(type_iapws97), intent(in), target :: water
            type(type_iapws06), intent(in), target :: ice

        end subroutine initialize_abst_thc
    end interface

    !--------------------------------------------------------------------------------
    type, extends(abst_thc) :: type_thc_1phase
    contains
        procedure, pass(self) :: calc_lambda_0_gp => calc_thc_gp_1phase
        procedure, pass(self) :: calc_lambda_dispersivity_gp => calc_thc_dispersivity_gp_1phase
    end type type_thc_1phase

    interface
        module subroutine calc_thc_gp_1phase(self, state, lambda, tangent)
            implicit none
            class(type_thc_1phase), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: lambda
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_thc_gp_1phase

        module subroutine calc_thc_dispersivity_gp_1phase(self, state, lambda)
            implicit none
            class(type_thc_1phase), intent(in) :: self
            type(type_state), intent(in) :: state
            type(type_state_thc), intent(inout) :: lambda

        end subroutine calc_thc_dispersivity_gp_1phase
    end interface

    type, extends(abst_thc) :: type_thc_2phase
    contains
        procedure, pass(self) :: calc_lambda_0_gp => calc_thc_gp_2phase
        procedure, pass(self) :: calc_lambda_dispersivity_gp => calc_thc_dispersivity_gp_2phase
    end type type_thc_2phase

    interface
        module subroutine calc_thc_gp_2phase(self, state, lambda, tangent)
            implicit none
            class(type_thc_2phase), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: lambda
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_thc_gp_2phase

        module subroutine calc_thc_dispersivity_gp_2phase(self, state, lambda)
            implicit none
            class(type_thc_2phase), intent(in) :: self
            type(type_state), intent(in) :: state
            type(type_state_thc), intent(inout) :: lambda

        end subroutine calc_thc_dispersivity_gp_2phase
    end interface

    type, extends(abst_thc) :: type_thc_3phase
    contains
        procedure, pass(self) :: calc_lambda_0_gp => calc_thc_gp_3phase
        procedure, pass(self) :: calc_lambda_dispersivity_gp => calc_thc_dispersivity_gp_3phase
    end type type_thc_3phase

    interface
        module subroutine calc_thc_gp_3phase(self, state, lambda, tangent)
            implicit none
            class(type_thc_3phase), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: lambda
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_thc_gp_3phase

        module subroutine calc_thc_dispersivity_gp_3phase(self, state, lambda)
            implicit none
            class(type_thc_3phase), intent(in) :: self
            type(type_state), intent(in) :: state
            type(type_state_thc), intent(inout) :: lambda

        end subroutine calc_thc_dispersivity_gp_3phase
    end interface

    type, extends(abst_thc) :: type_thc_4phase
    contains
        procedure, pass(self) :: calc_lambda_0_gp => calc_thc_gp_4phase
        procedure, pass(self) :: calc_lambda_dispersivity_gp => calc_thc_dispersivity_gp_4phase
    end type type_thc_4phase

    interface
        module subroutine calc_thc_gp_4phase(self, state, lambda, tangent)
            implicit none
            class(type_thc_4phase), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: lambda
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_thc_gp_4phase

        module subroutine calc_thc_dispersivity_gp_4phase(self, state, lambda)
            implicit none
            class(type_thc_4phase), intent(in) :: self
            type(type_state), intent(in) :: state
            type(type_state_thc), intent(inout) :: lambda

        end subroutine calc_thc_dispersivity_gp_4phase
    end interface

    interface
        module subroutine calc_lambda_dispersivity_abst_thc(self, lambda_0, htc_water, q_x, q_y, q_z, lambda)
            implicit none
            class(abst_thc), intent(in) :: self
            real(real64), intent(in) :: lambda_0
            real(real64), intent(in) :: htc_water
            real(real64), intent(in) :: q_x
            real(real64), intent(in) :: q_y
            real(real64), intent(in) :: q_z
            type(type_state_thc), intent(inout) :: lambda

        end subroutine calc_lambda_dispersivity_abst_thc
    end interface

end module materials_thermal_conductivity
