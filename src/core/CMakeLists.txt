# ============================================================
# Core Library
# ============================================================

add_library(SOLVER_CORE STATIC)

target_sources(SOLVER_CORE PRIVATE
    parallel_mpi.F90
    parallel_reduce.F90
)

# ============================================================
# 1. Foundational Sources
# ============================================================

# ------------------------------------------------------------
# Constants
# ------------------------------------------------------------
set(CORE_CONSTANTS_FYPP
    constants/constants_ids.fypp
    constants/constants_values.fypp
    constants/constants_error.fypp
)
fypp_F90("" "${CORE_CONSTANTS_FYPP}" CORE_CONSTANTS_GENERATED)

target_sources(SOLVER_CORE PRIVATE
    constants/physical_constants.F90
    constants/constants_base.F90
    constants/constants.F90
    ${CORE_CONSTANTS_GENERATED}
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------
set(CORE_VALIDATION_FYPP
    validation/check_length.fypp
    validation/check_nan.fypp
    validation/check_range.fypp
)
fypp_F90("" "${CORE_VALIDATION_FYPP}" CORE_VALIDATION_GENERATED)

target_sources(SOLVER_CORE PRIVATE
    validation/error.F90
    validation/check_io.F90
    validation/validation.F90
    ${CORE_VALIDATION_GENERATED}
)

# ------------------------------------------------------------
# Memory
# ------------------------------------------------------------
set(CORE_MEMORY_FYPP
    memory/allocate.fypp
    memory/deallocate.fypp
)
fypp_F90("-D MAX_RANK=3" "${CORE_MEMORY_FYPP}" CORE_MEMORY_GENERATED)

target_sources(SOLVER_CORE PRIVATE
    memory/memory.F90
    ${CORE_MEMORY_GENERATED}
)

# ------------------------------------------------------------
# Utils
# ------------------------------------------------------------
set(CORE_UTILS_FYPP
    utils/unique.fypp
    utils/findings.fypp
)
fypp_F90("" "${CORE_UTILS_FYPP}" CORE_UTILS_GENERATED)

target_sources(SOLVER_CORE PRIVATE
    utils/string_utils.F90
    utils/system_env.F90
    utils/utils.F90
    ${CORE_UTILS_GENERATED}
)

# ------------------------------------------------------------
# Types
# ------------------------------------------------------------
target_sources(SOLVER_CORE PRIVATE
    # Config
    types/config/config_base.F90
    types/config/conditions/config_boundary.F90
    types/config/conditions/config_initial.F90
    types/config/conditions/config_conditions.F90
    types/config/physics/config_physics_base.F90
    types/config/physics/config_constitutive.F90
    types/config/physics/config_swcc.F90
    types/config/physics/config_gcc.F90
    types/config/physics/config_physics.F90
    types/config/control/config_output_manager.F90
    types/config/control/config_control.F90
    types/config/control/config_time.F90
    # geometry
    types/geometry/coordinate_array.F90
    types/geometry/coordinate.F90
    # Physics
    types/physics/meteorology.F90
    types/physics/solar.F90
    types/physics/state.F90
    types/physics/state_thc.F90
    # Utils
    types/utils/datetime.F90
    # Core
    types/types.F90

)

# ============================================================
# 2. Interop (C / C++ / Fortran wrappers)
# ============================================================

target_sources(SOLVER_CORE PRIVATE

    # --- Utils ---
    interop/interop_utils.F90

    # --- PETSc ---
    interop/petscfe/petscfe_shim.c
    interop/petscfe/petscfe_wrapper.F90
    interop/petscds/petscds_shim.c
    interop/petscds/petscds_wrapper.F90

    interop/interop.F90
)

# ============================================================
# 3. Core Entry
# ============================================================

target_sources(SOLVER_CORE PRIVATE
    core.F90
)

# ============================================================
# 4. Target Properties
# ============================================================

# 生成モジュール用
target_include_directories(SOLVER_CORE
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)


# ============================================================
# 6. Project-specific flags (existing functions kept)
# ============================================================

enable_build_flags(SOLVER_CORE)
enable_thirdparty(SOLVER_CORE)

# ============================================================
# 7. Alias
# ============================================================

add_library(SOLVER::CORE ALIAS SOLVER_CORE)


target_link_libraries(SOLVER_CORE PUBLIC PETSc::PETSc)
