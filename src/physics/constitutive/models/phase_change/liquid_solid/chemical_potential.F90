!> @brief Chemical potential utilities for the freezing fringe.
!>
!> Provides thermodynamic quantities derived from the generalized Clapeyron equation
!> (ice-pressure-constant approximation):
!>
!> \[ \Psi_{ice}(T) = L_f \ln\!\left(\frac{T_K}{T_{r,K}}\right) \quad [J/kg] \]
!>
!> The **freezing interface** is the isotherm \(T = T_{high}(p_w)\) at which the
!> capillary and cryogenic suctions coincide:
!>
!> \[ T_{high} = T_r \exp\!\left(\frac{\Psi_w}{L_f}\right), \qquad
!>    \Psi_w = \frac{P_w}{\rho_w} \]
!>
!> These utilities support the freezing-point relations that
!> integrates the discontinuous suction-switch coefficients in cut elements.
module models_phase_change_chemical_potential
    use, intrinsic :: iso_fortran_env
    use :: constitutive_constants, only: &
        lf => latent_heat_fusion_water_0c, &
        Tf0_K => water_freezing_point_at_standard_atmospheric_pressure_k, &
        Tf0 => water_freezing_point_at_standard_atmospheric_pressure
    implicit none
    private

    public :: calc_psi_ice
    public :: calc_dpsi_ice_dT
    public :: calc_T_high_celsius

contains

    !> @brief Specific chemical potential of ice relative to liquid at \(T_r\): \(\Psi_{ice} = L_f \ln(T_K/T_{r,K})\)
    !>
    !> \[ \Psi_{ice}(T) = L_f \ln\!\left(\frac{T_K}{T_{r,K}}\right) \quad [J/kg] \]
    !>
    !> Returns zero for T >= 0 degC (no ice contribution above freezing point).
    !> Assumptions: T in Celsius; L_f constant (evaluated at 0 degC).
    pure subroutine calc_psi_ice(T_C, psi_ice)
        implicit none
        real(real64), intent(in)  :: T_C      !< Temperature [degC]
        real(real64), intent(inout) :: psi_ice !< \(\Psi_{ice}\) [J/kg]

        real(real64) :: T_K

        if (T_C >= Tf0) then
            psi_ice = 0.0d0
            return
        end if
        T_K = T_C + Tf0_K
        if (T_K <= 0.0d0) then
            psi_ice = 0.0d0
            return
        end if
        psi_ice = lf * log(T_K / Tf0_K)
    end subroutine calc_psi_ice

    !> @brief Derivative \( d\Psi_{ice}/dT = L_f / T_K \) [J/(kg K)]
    !>
    !> Returns zero for T >= 0 degC.
    pure subroutine calc_dpsi_ice_dT(T_C, dpsi_dT)
        implicit none
        real(real64), intent(in)  :: T_C      !< Temperature [degC]
        real(real64), intent(inout) :: dpsi_dT !< \(d\Psi_{ice}/dT\) [J/(kg K)]

        real(real64) :: T_K

        if (T_C >= Tf0) then
            dpsi_dT = 0.0d0
            return
        end if
        T_K = T_C + Tf0_K
        if (T_K <= tiny(1.0d0)) then
            dpsi_dT = 0.0d0
            return
        end if
        dpsi_dT = lf / T_K
    end subroutine calc_dpsi_ice_dT

    !> @brief Nucleation temperature \(T_{high}\) from pore-water pressure.
    !>
    !> \[ T_{high} = T_r \exp\!\left(\frac{P_w}{\rho_w L_f}\right) - 273.15 \quad [^\circ\text{C}] \]
    !>
    !> For unsaturated soil with \(P_w < 0\): \(T_{high} < 0\).
    !> For \(P_w = 0\): \(T_{high} = 0\) degC (atmospheric).
    !>
    !> \param P_w    Pore-water pressure [Pa] (negative for suction)
    !> \param rho_w  Liquid water density [kg/m^3]
    !> \param T_high Nucleation temperature [degC]
    pure subroutine calc_T_high_celsius(P_w, rho_w, T_high)
        implicit none
        real(real64), intent(in)  :: P_w    !< Pore-water pressure [Pa]
        real(real64), intent(in)  :: rho_w  !< Water density [kg/m^3]
        real(real64), intent(inout) :: T_high !< Nucleation temperature [degC]

        real(real64) :: psi_w, exponent

        if (rho_w <= tiny(1.0d0) .or. lf <= tiny(1.0d0)) then
            T_high = Tf0
            return
        end if

        psi_w    = P_w / rho_w
        exponent = psi_w / lf
        ! Clamp exponent to avoid exp overflow / underflow
        exponent = max(-50.0d0, min(exponent, 0.0d0))
        T_high = Tf0_K * exp(exponent) - Tf0_K
    end subroutine calc_T_high_celsius

end module models_phase_change_chemical_potential
