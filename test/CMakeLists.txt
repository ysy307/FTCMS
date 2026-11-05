# Ensure log output directory exists
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/log)

add_executable(petsc_smoke test_petsc_smoke.F90)
set_target_properties(petsc_smoke PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
target_link_libraries(petsc_smoke PRIVATE
    PETSc::PETSc
    MPI::MPI_Fortran
)

add_test(NAME petsc_smoke_serial COMMAND $<TARGET_FILE:petsc_smoke>)

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
            COMMAND sh -c "$<TARGET_FILE:${target_name}> > ${log_file}"
            DEPENDS ${target_name}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            VERBATIM
        )
    endif()
endfunction()

# --------------------------------------------------------------------
# Register each test with its specific dependencies
# --------------------------------------------------------------------

add_solver_test(constitutive
    LIBS
        SOLVER::CORE
        SOLVER::NUMERICAL_MATH
        SOLVER::INPUT
        SOLVER::CONSTITUTIVE
)

# --------------------------------------------------------------------
# Verification cases: FTCMS runs on small projects under test/data/projects,
# each checked against one stated value (test/verification.py).
# --------------------------------------------------------------------
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(FTCMS_VERIFICATION_CASES
    conduction
    radiation
    hydrostatic_mixed
    hydrostatic_lagrange
    inflow_mixed_lagrange
    layers_mixed
    layers_lagrange
    water_table_mixed
    water_table_lagrange
    jacobian_mixed
    jacobian_lagrange
    jacobian_nodal
    freezing_mixed
    freezing_nodal
    atmosphere_surface
    refuse_mixed_heat
    refuse_nodal_isothermal
    refuse_vapour
    refuse_segregation
)
foreach(case IN LISTS FTCMS_VERIFICATION_CASES)
    add_test(NAME ftcms_${case}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verification.py ${case}
            --ftcms $<TARGET_FILE:FTCMS>
            --work ${CMAKE_BINARY_DIR}/verification
            --mpirun ${MPIEXEC_EXECUTABLE})
    set_tests_properties(ftcms_${case} PROPERTIES TIMEOUT 1800 LABELS verification)
endforeach()

# The same statements on two ranks: distribution, reductions, ghost exchange.
if(ENABLE_MPI)
    foreach(case IN ITEMS freezing_mixed freezing_nodal inflow_mixed_lagrange layers_mixed atmosphere_surface)
        add_test(NAME ftcms_${case}_mpi_2
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verification.py ${case}
                --ftcms $<TARGET_FILE:FTCMS>
                --work ${CMAKE_BINARY_DIR}/verification
                --mpirun ${MPIEXEC_EXECUTABLE} --ranks 2)
        set_tests_properties(ftcms_${case}_mpi_2 PROPERTIES TIMEOUT 1800 LABELS verification)
    endforeach()
endif()
