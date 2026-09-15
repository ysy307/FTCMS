# ============================================================
# PetscDS / TS Application Library
# ============================================================
add_library(SOLVER_PLEXFEM STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_PLEXFEM PRIVATE
    plexfem_physics.F90
    plexfem_pointwise.F90
    plexfem_interface.F90
    plexfem_setup.F90
    plexfem_run.F90
)

# ------------------------------------------------------------
# Include directory (Fortran module path)
# ------------------------------------------------------------
target_include_directories(SOLVER_PLEXFEM PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
target_link_libraries(SOLVER_PLEXFEM PUBLIC
    SOLVER::GOVERNING
    SOLVER::CONTROL
    SOLVER::CORE
    SOLVER::INPUT
)

find_package(PkgConfig REQUIRED)

# PETSc: the discrete system, the boundary conditions and the time stepper all
# come from it, through the shim in src/core/interop/petscds.
set(ENV{PKG_CONFIG_PATH}
    "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(PETSC_PLEXFEM REQUIRED IMPORTED_TARGET PETSc)
target_link_libraries(SOLVER_PLEXFEM PUBLIC PkgConfig::PETSC_PLEXFEM)

enable_build_flags(SOLVER_PLEXFEM)
enable_thirdparty(SOLVER_PLEXFEM)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::PLEXFEM ALIAS SOLVER_PLEXFEM)
