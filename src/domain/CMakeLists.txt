# ============================================================
# Domain Library: the mesh as a PETSc DMPlex
# ============================================================
add_library(SOLVER_DOMAIN STATIC)

target_sources(SOLVER_DOMAIN PRIVATE
    mesh/mesh_plex.F90
)

target_include_directories(SOLVER_DOMAIN PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(SOLVER_DOMAIN PUBLIC
    SOLVER::CORE
)

enable_build_flags(SOLVER_DOMAIN)
enable_thirdparty(SOLVER_DOMAIN)

add_library(SOLVER::DOMAIN ALIAS SOLVER_DOMAIN)
