# ============================================================
# Mixed hydraulic (RT0 flux / P0 pressure) solver library
# ============================================================
add_library(SOLVER_MIXED_HYDRAULIC STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_MIXED_HYDRAULIC PRIVATE
    mixed_hydraulic_physics.F90
    mixed_hydraulic_pointwise.F90
    mixed_hydraulic_balance.F90
    mixed_hydraulic_interface.F90
    mixed_hydraulic_setup.F90
    mixed_hydraulic_solve.F90
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
# Only the interop bridges and the PETSc library: this solver carries its own
# physics and does not reach into the application layer.
target_link_libraries(SOLVER_MIXED_HYDRAULIC PUBLIC
    SOLVER::CORE
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
