module constitutive_constants
    use, intrinsic :: iso_fortran_env
    use :: core_physical_constants, only: reference_water_density, gravity_acceleration
    implicit none
    private

    public :: reference_water_density
    public :: gravity_acceleration

    !> Mathematical constant pi
    real(real64), parameter, public :: circle_ratio = 3.141592653589793d0

    !> Conversion from Celsius to Kelvin [K]
    real(real64), parameter, public :: celsius_to_kelvin = 273.15d0

    !> Standard atmospheric pressure [Pa]
    real(real64), parameter, public :: standard_atmospheric_pressure = 101325.0d0

    !> Universal gas constant [J/K/mol]
    real(real64), parameter, public :: universal_gas_constant = 8.31446261815324d0

    !------------------------------------------------------------------------------------------
    ! For water properties
    !------------------------------------------------------------------------------------------
    !> The specific gas constant of ordinary water [J/kg/K]
    real(real64), parameter, public :: specific_gas_constant_water = 461.526d0
    !> Water triple point temperature [K]
    real(real64), parameter, public :: water_triple_point_temperature = 273.16d0
    !> Water triple point pressure [Pa]
    real(real64), parameter, public :: water_triple_point_pressure = 611.657d0
    !> Water triple point density [kg/m^3]
    real(real64), parameter, public :: water_triple_point_density = 999.793d0
    !> Water critical point temperature [K]
    real(real64), parameter, public :: water_critical_point_temperature = 647.096d0
    !> Water critical point pressure [Pa]
    real(real64), parameter, public :: water_critical_point_pressure = 22.06464d6
    !> Water critical point density [kg/m^3]
    real(real64), parameter, public :: water_critical_point_density = 322.0d0

    !> Latent heat of fusion for water at 0 degC [J/kg]
    real(real64), parameter, public :: latent_heat_fusion_water_0C = 334560.0d0
    !> Water freezing point at standard atmospheric pressure [degC]
    real(real64), parameter, public :: water_freezing_point_at_standard_atmospheric_pressure = 0.0d0
    !>
    real(real64), parameter, public :: water_freezing_point_at_standard_atmospheric_pressure_K = &
                                       water_freezing_point_at_standard_atmospheric_pressure + celsius_to_kelvin

    !> Molar mass of water [kg/mol]
    real(real64), parameter, public :: molar_mass_water = 0.01801528d0

    !> Stefan-Boltzman constant [W/m^2/K^4]
    real(real64), parameter, public :: stefan_boltzmann_constant = 5.670374419d-8

    real(real64), parameter, public :: min_vapor_density = 1.0d-8

    !> Differentiability guard [Pa] for the generalized suction
    !> \(p_c^* = \max(\psi_{cap}, \psi_{cryo})\) and its subgradient weights.
    !> Single source of truth for every smooth-max evaluation of the freezing
    !> switch.  It is NOT a physical smoothing scale: elements cut by the
    !> so quadrature points never straddle the switch; this epsilon only keeps
    !> the derivative weights defined exactly at \(\psi_{cap} = \psi_{cryo}\).
    real(real64), parameter, public :: suction_blend_epsilon = 1.0d2

end module constitutive_constants
