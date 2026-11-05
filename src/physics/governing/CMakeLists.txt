# ============================================================
# Governing Equation Library
# ============================================================
add_library(SOLVER_GOVERNING STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_GOVERNING PRIVATE
    # ====================================================
    # 1. Boundary Conditions
    # ====================================================
    condition/boundary/data_provider/data_provider_interface.F90
    condition/boundary/data_provider/data_provider_constant.F90
    condition/boundary/data_provider/data_provider_table.F90
    condition/boundary/data_provider/data_provider_dynamic.F90
    condition/boundary/data_provider/data_provider_factory.F90
    condition/boundary/strategy/strategy_interface.F90
    condition/boundary/strategy/strategy_base.F90
    condition/boundary/strategy/strategy_dirichlet.F90
    condition/boundary/strategy/strategy_flux.F90
    condition/boundary/strategy/strategy_robin.F90
    condition/boundary/strategy/strategy_switching.F90
    condition/boundary/strategy/strategy_factory.F90
    condition/boundary/boundary_manager.F90
    condition/boundary/boundary.F90
    # ====================================================
    # 4. Thermal
    # ====================================================
    thermal/thermal_interface.F90
    thermal/thermal_base.F90
    thermal/thermal_coefficients.F90
    # ====================================================
    # 5. Hydraulic
    # ====================================================
    hydraulic/hydraulic_interface.F90
    hydraulic/hydraulic_base.F90
    hydraulic/hydraulic_coefficients.F90
    # ====================================================
    # 6. Atmosphere / Data Assimilation (ETKF)
    # ====================================================
    atmosphere/ensemble_manager.F90
    atmosphere/observation_manager.F90
    atmosphere/etkf_solver.F90
    atmosphere/assimilation_controller.F90
    atmosphere/atmosphere.F90
    # ====================================================
    # \
    # ====================================================
    governing.F90


)

# ------------------------------------------------------------
# Include directory
# ------------------------------------------------------------
target_include_directories(SOLVER_GOVERNING PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
target_link_libraries(SOLVER_GOVERNING PUBLIC
    SOLVER::CORE
    SOLVER::INPUT
    SOLVER::CONSTITUTIVE
)

enable_build_flags(SOLVER_GOVERNING)
enable_thirdparty(SOLVER_GOVERNING)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::GOVERNING ALIAS SOLVER_GOVERNING)
