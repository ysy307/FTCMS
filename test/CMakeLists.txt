# Ensure log output directory exists
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/log)

# PETSc smoke test dependency
find_package(PkgConfig REQUIRED)

set(PETSC_PREFIX "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}")
set(ENV{PKG_CONFIG_PATH} "${PETSC_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

add_executable(petsc_smoke test_petsc_smoke.F90)
set_target_properties(petsc_smoke PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
target_link_libraries(petsc_smoke PRIVATE
    PkgConfig::PETSC
    MPI::MPI_Fortran
)

add_test(NAME petsc_smoke_serial COMMAND $<TARGET_FILE:petsc_smoke>)

add_executable(vtu_probe vtu_probe.F90)
target_link_libraries(vtu_probe PRIVATE PkgConfig::PETSC MPI::MPI_Fortran)

add_executable(dmplex_probe test_dmplex_probe.F90)
set_target_properties(dmplex_probe PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
target_link_libraries(dmplex_probe PRIVATE PkgConfig::PETSC MPI::MPI_Fortran)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME petsc_smoke_mpi_2
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:petsc_smoke>
    )
endif()

# ====================================================================
# Function: add_solver_test
# ====================================================================
function(add_solver_test test_name)
    set(options "")
    set(oneValueArgs "")
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(target_name "test_${test_name}")

    # ソースファイルが指定されていない場合はデフォルトパスを使用
    if(ARG_SOURCES)
        set(test_sources ${ARG_SOURCES})
    else()
        set(test_sources "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.F90")
        if(NOT EXISTS "${test_sources}")
            message(FATAL_ERROR "Test source file not found: ${test_sources}")
        endif()
    endif()

    add_executable(${target_name} ${test_sources})
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE Fortran)

    # Apply general build flags and third-party libraries by default
    enable_build_flags(${target_name})
    enable_thirdparty(${target_name})

    # Link internal solver libraries passed via LIBS
    if(ARG_LIBS)
        target_link_libraries(${target_name} PRIVATE ${ARG_LIBS})
    endif()

    # Link test-drive framework
    target_link_libraries(${target_name} PRIVATE test-drive::test-drive)

    # Register to CTest
    add_test(NAME ${target_name} COMMAND $<TARGET_FILE:${target_name}>)

    # Custom run target
    set(log_dir "${PROJECT_SOURCE_DIR}/log")
    set(log_file "${log_dir}/results_${test_name}.md")

    if(ENABLE_MPI AND MPI_Fortran_FOUND)
        add_custom_target(run_${test_name}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${log_dir}
            COMMAND sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${target_name}> > ${log_file}"
            DEPENDS ${target_name}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            VERBATIM
        )
    else()
        add_custom_target(run_${test_name}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${log_dir}
            COMMAND sh -c "env OMP_NUM_THREADS=4 $<TARGET_FILE:${target_name}> > ${log_file}"
            DEPENDS ${target_name}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            VERBATIM
        )
    endif()
endfunction()

# --------------------------------------------------------------------
# Register each test with its specific dependencies
# --------------------------------------------------------------------

# test_matrix: Fypp-generated lifecycle and assembly tests
fypp_F90("" "test_matrix.fypp" TEST_MATRIX_GENERATED)
add_solver_test(matrix
    SOURCES
        ${TEST_MATRIX_GENERATED}
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
)

# test_vector: Fypp-generated tests for type_vector_dp and type_vector_int
fypp_F90("" "test_vector.fypp" TEST_VECTOR_GENERATED)
add_solver_test(vector
    SOURCES
        ${TEST_VECTOR_GENERATED}
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
)

# test_field: domain-free system layer (Jacobian/residual via type_system_topology)
add_solver_test(field
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
        SOLVER::SYSTEM
)

add_solver_test(constitutive
    LIBS
        SOLVER::CORE
        SOLVER::NUMERICAL_MATH
        SOLVER::INPUT
        SOLVER::CONSTITUTIVE
)

add_solver_test(ats
    LIBS
        SOLVER::CORE
        SOLVER::CONTROL
)

add_solver_test(numerical
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
        SOLVER::NUMERICAL_MATH
        SOLVER::SOLVER
)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME test_numerical_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_numerical>)
    add_test(NAME test_numerical_mpi_4
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_numerical>)
endif()

# test_main requires all solver libraries
add_solver_test(main
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
        SOLVER::SOLVER
        SOLVER::SYSTEM
        SOLVER::DOMAIN
        SOLVER::INPUT
        SOLVER::OUTPUT
        SOLVER::CONSTITUTIVE
        SOLVER::GOVERNING
        SOLVER::FTCMS
)

add_solver_test(coupling_assembly
    LIBS
        SOLVER::CORE
        SOLVER::LINALG
        SOLVER::SOLVER
        SOLVER::SYSTEM
        SOLVER::DOMAIN
        SOLVER::INPUT
        SOLVER::OUTPUT
        SOLVER::CONSTITUTIVE
        SOLVER::GOVERNING
        SOLVER::FTCMS
)
