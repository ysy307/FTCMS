# ============================================================
# FTCMS Application Library
# ============================================================
add_library(SOLVER_FTCMS STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_FTCMS PRIVATE
    ftcms_darcy.F90
    ftcms_closure.F90
    ftcms_kernels.F90
    ftcms_kernels_lagrange.F90
    ftcms_kernels_advection.F90
    ftcms_nodal_flux.F90
    ftcms_boundary_registry.F90
    ftcms_interface.F90
    ftcms_setup.F90
    ftcms_step.F90
    ftcms_output.F90
    ftcms_assimilation.F90
    ftcms.F90
)

# ------------------------------------------------------------
# Include directory (Fortran module path)
# ------------------------------------------------------------
target_include_directories(SOLVER_FTCMS PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
# The governing layer supplies the constitutive models the closure evaluates
# and the boundary strategies; the input layer the project definition.
target_link_libraries(SOLVER_FTCMS PUBLIC
    SOLVER::GOVERNING
    SOLVER::CORE
    SOLVER::INPUT
)

enable_build_flags(SOLVER_FTCMS)
enable_thirdparty(SOLVER_FTCMS)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::FTCMS ALIAS SOLVER_FTCMS)
