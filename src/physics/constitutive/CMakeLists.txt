# ============================================================
# Constitutive Library
# ============================================================
add_library(SOLVER_CONSTITUTIVE STATIC)

# ------------------------------------------------------------
# Closed-form laws: SymPy states them, fypp expands them
# ------------------------------------------------------------
fypp_require_python_module(sympy)
fypp_F90("-m;sympy.printing.fortran;-m;textwrap" "laws/constitutive_laws.fypp" CONSTITUTIVE_LAWS_GENERATED)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_CONSTITUTIVE PRIVATE
    # ====================================================
    # 1. Base
    # ====================================================
    constitutive_constants.F90
    ${CONSTITUTIVE_LAWS_GENERATED}
    constitutive_base.F90

    # ====================================================
    # 2. Materials
    # ====================================================
    materials/materials_base.F90

    materials/density/density_interface.F90
    materials/density/density_base.F90
    materials/density/density_multiphase.F90

    materials/specific_heat/specific_heat_interface.F90
    materials/specific_heat/specific_heat_base.F90
    materials/specific_heat/specific_heat_multiphase.F90

    materials/thermal_conductivity/thermal_conductivity_interface.F90
    materials/thermal_conductivity/thermal_conductivity_base.F90
    materials/thermal_conductivity/thermal_conductivity_multiphase.F90

    materials/heat_capacity/heat_capacity_interface.F90
    materials/heat_capacity/heat_capacity_base.F90
    materials/heat_capacity/heat_capacity_multiphase.F90

    materials/materials_manager.F90

    # ====================================================
    # 3. Models
    # ====================================================
    models/permeability/hcf_interface.F90
    models/permeability/hcf_base.F90
    models/permeability/hcf_bc.F90
    models/permeability/hcf_vg.F90
    models/permeability/hcf_ko.F90
    models/permeability/hcf_mvg.F90
    models/permeability/hcf_durner.F90
    models/permeability/hcf_dvgch.F90
    models/permeability/hcf_impedance.F90
    models/permeability/hcf_viscosity.F90
    models/permeability/hcf_vapor.F90

    models/water_retention/wrf_interface.F90
    models/water_retention/wrf_base.F90
    models/water_retention/wrf_bc.F90
    models/water_retention/wrf_vg.F90
    models/water_retention/wrf_ko.F90
    models/water_retention/wrf_mvg.F90
    models/water_retention/wrf_durner.F90
    models/water_retention/wrf_dvgch.F90

    models/phase_change/liquid_solid/gcc_interface.F90
    models/phase_change/liquid_solid/gcc_base.F90
    models/phase_change/liquid_solid/fusion.F90
    models/phase_change/liquid_solid/chemical_potential.F90
    models/phase_change/liquid_vapor/vaporization.F90
    models/phase_change/segregation/segregation_potential.F90
    models/phase_change/phase_systems.F90

    models/models_manager.F90

    # ====================================================
    # 4. Service Layer
    # ====================================================
    constitutive_manager.F90
    constitutive.F90
)

# ------------------------------------------------------------
# Include directory
# ------------------------------------------------------------
target_include_directories(SOLVER_CONSTITUTIVE PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
target_link_libraries(SOLVER_CONSTITUTIVE PUBLIC
    SOLVER::CORE
    SOLVER::NUMERICAL_MATH
)

enable_build_flags(SOLVER_CONSTITUTIVE)
enable_thirdparty(SOLVER_CONSTITUTIVE)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::CONSTITUTIVE ALIAS SOLVER_CONSTITUTIVE)
