# ============================================================
# Mixed hydraulic (RT0 flux / P0 pressure) solver library
# ============================================================
add_library(SOLVER_MIXED_HYDRAULIC STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_MIXED_HYDRAULIC PRIVATE
    mixed_hydraulic_physics.F90
    mixed_hydraulic_closure.F90
    mixed_hydraulic_frozen.F90
    mixed_hydraulic_pointwise.F90
    mixed_hydraulic_balance.F90
    mixed_hydraulic_gauge.F90
    mixed_hydraulic_interface.F90
    mixed_hydraulic_setup.F90
    mixed_hydraulic_solve.F90
    mixed_hydraulic_column_pointwise.F90
    mixed_hydraulic_column_primal.F90
    mixed_hydraulic_column_nodal.F90
    mixed_hydraulic_column.F90
)

# ------------------------------------------------------------
# Include directory (Fortran module path)
# ------------------------------------------------------------
target_include_directories(SOLVER_MIXED_HYDRAULIC PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
# The interop bridges, PETSc, and the constitutive models whose retention curve
# and hydraulic conductivity function the cell-local closure evaluates. Not the
# application layer: the mixed solver calls the material laws directly rather
# than through the monolithic path.
#
# The governing layer is linked for the frozen closure only. Its contract is
# (T, p_w, material) -> (W, enthalpy, liquid pressure, mobility), and every one
# of those is a procedure of type_hydraulic or type_thermal: evaluating them
# through the repository's own governing objects is what keeps the mixed solver
# from carrying a second statement of the freezing physics.
target_link_libraries(SOLVER_MIXED_HYDRAULIC PUBLIC
    SOLVER::CORE
    SOLVER::CONSTITUTIVE
    SOLVER::GOVERNING
)

find_package(PkgConfig REQUIRED)

set(ENV{PKG_CONFIG_PATH}
    "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(PETSC_MIXED_HYDRAULIC REQUIRED IMPORTED_TARGET PETSc)
target_link_libraries(SOLVER_MIXED_HYDRAULIC PUBLIC PkgConfig::PETSC_MIXED_HYDRAULIC)

enable_build_flags(SOLVER_MIXED_HYDRAULIC)
enable_thirdparty(SOLVER_MIXED_HYDRAULIC)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::MIXED_HYDRAULIC ALIAS SOLVER_MIXED_HYDRAULIC)
