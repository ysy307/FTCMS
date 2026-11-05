# ============================================================
# Input Library
# ============================================================
add_library(SOLVER_INPUT STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_INPUT PRIVATE
    # ----------------------------
    # Reader layer
    # ----------------------------
    reader/input_base.F90

    # --- Basic ---
    reader/basic/basic_base.F90
    reader/basic/basic_interface.F90
    reader/basic/basic_simulation_setting.F90
    reader/basic/basic_analysis_controls.F90
    reader/basic/basic_geometry_settings.F90
    reader/basic/basic_materials.F90
    reader/basic/basic_solver_settings.F90
    reader/basic/basic_discretization_settings.F90

    # --- Conditions ---
    reader/conditions/conditions_base.F90
    reader/conditions/conditions_interface.F90
    reader/conditions/conditions_time_controls.F90
    reader/conditions/conditions_boundary.F90
    reader/conditions/conditions_initial.F90

    # --- Output ---
    reader/output/output_interface.F90
    reader/output/output_fields.F90
    reader/output/output_history.F90
    reader/output/output_standard.F90

    # --- Geometry ---
    reader/geometry/geometry_base.F90
    reader/geometry/geometry_interface.F90

    # Input manager
    reader/assimilation/assimilation_settings.F90
    reader/input_manager.F90

    # ----------------------------
    # Translator layer
    # ----------------------------
    translator/translator_interface.F90
    translator/translator_basic.F90
    translator/translator_conditions.F90
    translator/translator_output.F90

    # ----------------------------
    # Entry point
    # ----------------------------
    input.F90
)

# ------------------------------------------------------------
# Fortran module output directory
# ------------------------------------------------------------
target_include_directories(SOLVER_INPUT PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
target_link_libraries(SOLVER_INPUT PUBLIC
    SOLVER::CORE
    SOLVER::DOMAIN
)

# ------------------------------------------------------------
# Project-specific flags (kept as-is)
# ------------------------------------------------------------
enable_build_flags(SOLVER_INPUT)
enable_thirdparty(SOLVER_INPUT)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::INPUT ALIAS SOLVER_INPUT)
