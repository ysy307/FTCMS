# ============================================================
# Domain Library
# ============================================================
add_library(SOLVER_DOMAIN STATIC)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
target_sources(SOLVER_DOMAIN PRIVATE
    mesh/mesh_plex.F90
    # ====================================================
    # 1. Mesh (Finite Element Core)
    # ====================================================
    # Abstract finite element definitions
    fe/base/fe_interface.F90
    fe/base/fe_base.F90
    fe/base/fe_interpolation.F90
    fe/base/fe_compute.F90
    fe/base/fe_coefficients.F90
    # FE utilities
    fe/fe_integration.F90
    fe/fe_petsc.F90
    fe/fe_factory.F90
    fe/fe_manager.F90
    fe/fe.F90
    # ====================================================
    # 2. Adjacency
    # ====================================================
    adjacency/node.F90
    adjacency/node_element.F90
    adjacency/adjacency.F90
    # ====================================================
    # 3. Multicoloring
    # ====================================================
    multicoloring/multicoloring.F90
    # ====================================================
    # 4. Domain Components
    # ====================================================
    components/domain_nodes.F90
    components/domain_elements.F90
    components/domain_boundaries.F90
    # ====================================================
    # 5. Domain Service Layer
    # ====================================================
    domain_manager.F90
    domain.F90
)

# ------------------------------------------------------------
# Fortran module output directory
# ------------------------------------------------------------
target_include_directories(SOLVER_DOMAIN PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
target_link_libraries(SOLVER_DOMAIN PUBLIC
    SOLVER::CORE
    SOLVER::LINALG
)

enable_build_flags(SOLVER_DOMAIN)
enable_thirdparty(SOLVER_DOMAIN)

# ------------------------------------------------------------
# Namespaced alias
# ------------------------------------------------------------
add_library(SOLVER::DOMAIN ALIAS SOLVER_DOMAIN)

find_package(PkgConfig REQUIRED)
set(ENV{PKG_CONFIG_PATH} "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(PETSC_DOMAIN REQUIRED IMPORTED_TARGET PETSc)
target_link_libraries(SOLVER_DOMAIN PUBLIC PkgConfig::PETSC_DOMAIN)
