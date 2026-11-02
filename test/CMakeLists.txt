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

add_executable(test_mixed_hydraulic test_mixed_hydraulic.F90)
set_target_properties(test_mixed_hydraulic PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
target_link_libraries(test_mixed_hydraulic PRIVATE SOLVER::MIXED_HYDRAULIC MPI::MPI_Fortran)
enable_build_flags(test_mixed_hydraulic)
add_test(NAME mixed_space_serial COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case space)
set_tests_properties(mixed_space_serial PROPERTIES TIMEOUT 60)

# The same construction question for the rectangle: a manufactured state that
# each of the three spaces carries exactly, read back through the audit's own
# reconstruction. Every comparison it makes is against zero, so it separates a
# wrong reference cell, a wrong pushforward and a simplex element on a tensor
# cell from a discretisation error.
add_test(NAME mixed_rectangle
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case rectangle)
set_tests_properties(mixed_rectangle PROPERTIES TIMEOUT 60)

# Task 2 steady references. The driver installs the direct LU/MUMPS
# configuration itself and asserts the solver it actually ran with, so these
# entries carry only the case selection.
add_test(NAME mixed_darcy_constant
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case darcy_constant)
add_test(NAME mixed_hydrostatic
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case hydrostatic)
add_test(NAME mixed_layers
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case layers -mixed_contrast 1e7)
# The same code path at contrast one: the control that separates a coefficient
# jump failure from a failure the layered geometry has on its own.
add_test(NAME mixed_layers_equal
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case layers -mixed_contrast 1)
add_test(NAME mixed_flux_boundary
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case flux_boundary)
add_test(NAME mixed_jacobian
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case jacobian)
set_tests_properties(
    mixed_darcy_constant
    mixed_hydrostatic
    mixed_layers
    mixed_layers_equal
    mixed_flux_boundary
    mixed_jacobian
    PROPERTIES TIMEOUT 120
)

# The refusal is the result. WILL_FAIL alone cannot tell the expected refusal
# from a run that accepted the zero conductivity and then failed its own
# assertion, so the diagnostic itself is what this test matches.
add_test(NAME mixed_reject_mobility
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case reject_mobility)
set_tests_properties(mixed_reject_mobility PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "An impermeable region is a different model"
)

# The same shape of test for the geometry gate: a quadrilateral whose map from
# the reference square is bilinear rather than affine is refused, because that
# is the assumption the audit's reconstruction and the element loop both make.
add_test(NAME mixed_reject_warped
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case reject_warped)
set_tests_properties(mixed_reject_warped PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "triangles and axis-aligned rectangles only"
)

add_executable(vtu_probe vtu_probe.F90)
target_link_libraries(vtu_probe PRIVATE PkgConfig::PETSC MPI::MPI_Fortran)

add_executable(dmplex_probe test_dmplex_probe.F90)
set_target_properties(dmplex_probe PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
target_link_libraries(dmplex_probe PRIVATE PkgConfig::PETSC MPI::MPI_Fortran)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME mixed_space_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS} -mixed_case space
    )
    set_tests_properties(mixed_space_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 60)

    # The strip split across two ranks, so the shared facets' normal traces are
    # computed on different ranks and compared as two numbers. Its two layers
    # put interior facets across the manufactured field's own direction of
    # transport, where the trace is nonzero, and not only along it.
    add_test(NAME mixed_rectangle_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS} -mixed_case rectangle
    )
    set_tests_properties(mixed_rectangle_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 60)

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

# Task 2's gauge cases and Task 3's conservative nonlinear storage, failed-step
# rollback, pressure-dependent conductivity and per-cell material routing. Each
# runs serially and on two ranks, because the accepted history and the storage
# audit are both per-rank state.
foreach(gauge_case IN ITEMS all_flux all_flux_incompatible storage_pressure storage rollback
        richards material_interface freeze_closure)
    add_test(NAME mixed_${gauge_case}
        COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case ${gauge_case})
    set_tests_properties(mixed_${gauge_case} PROPERTIES TIMEOUT 60)
    if(ENABLE_MPI AND MPI_Fortran_FOUND)
        add_test(NAME mixed_${gauge_case}_mpi_2
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS} -mixed_case ${gauge_case})
        set_tests_properties(mixed_${gauge_case}_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 60)
    endif()
endforeach()

# Task 3's mesh and step ladders. The serial case carries the whole evidence:
# the two studies are separate on purpose, so running them on several rank
# counts would repeat the same discretisation argument rather than extend it.
# Only the window projection and its reductions are new sampling logic, so that
# is what the two-rank entry re-runs.
add_test(NAME mixed_refinement
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case refinement)
set_tests_properties(mixed_refinement PROPERTIES TIMEOUT 900)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME mixed_refinement_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS}
            -mixed_case refinement -mixed_refine_scope sample)
    set_tests_properties(mixed_refinement_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 120)
endif()

# The coupled column. `frozen_column` is the bounded regression -- a warm
# control that may not move water, then two hours of cooling that has to form
# ice, draw water up against gravity and conserve it -- and it is the case that
# guards the adaptive step and its rollback, so it runs on two ranks as well.
add_test(NAME mixed_frozen_column
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case frozen_column)
set_tests_properties(mixed_frozen_column PROPERTIES TIMEOUT 300)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME mixed_frozen_column_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS}
            -mixed_case frozen_column)
    set_tests_properties(mixed_frozen_column_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 300)
endif()

# The same column on rectangles, on the same vertex grid. The assertions are
# the triangular case's own, unchanged: what this entry establishes is that the
# tensor element conserves, freezes and draws water upward as well, which is
# the precondition for reading anything into a comparison of the two.
add_test(NAME mixed_frozen_column_quad
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case frozen_column
        -column_cell_type quadrilateral)
set_tests_properties(mixed_frozen_column_quad PROPERTIES TIMEOUT 300)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_test(NAME mixed_frozen_column_quad_mpi_2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS}
            -mixed_case frozen_column -column_cell_type quadrilateral)
    set_tests_properties(mixed_frozen_column_quad_mpi_2 PROPERTIES PROCESSORS 2 TIMEOUT 300)
endif()

# The experiment itself, to its first reported horizon. The registered entry
# stops at 12 h because that is what a test may cost; the 24 h and 50 h
# horizons are the same case run with `-mizo_final_hours`, and the profiles it
# is compared on come from `-mizo_profile_prefix`. Serial: the profile export
# is one rank's by construction, and the coupled solve's parallel path is what
# `mixed_frozen_column_mpi_2` covers.
add_test(NAME mixed_mizo_column
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case mizo_column -mizo_final_hours 12)
set_tests_properties(mixed_mizo_column PROPERTIES TIMEOUT 900)

# The cell profile export. The serial entry writes the reference and the two
# parallel ones are compared against it, so the three are ordered by a fixture
# rather than left to CTest's scheduling. Eight ranks is the point of the 4x16
# grid: it is the smallest of the experiment's own strips on which eight ranks
# each own real cells, so the gather is exercised and not short-circuited.
set(profile_export_dir "${CMAKE_CURRENT_BINARY_DIR}/profile_export")
file(MAKE_DIRECTORY ${profile_export_dir})

add_test(NAME mixed_profile_export
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case profile_export
        -export_prefix ${profile_export_dir}/serial)
set_tests_properties(mixed_profile_export PROPERTIES
    TIMEOUT 300
    FIXTURES_SETUP profile_export_reference)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    foreach(export_ranks IN ITEMS 2 8)
        add_test(NAME mixed_profile_export_mpi_${export_ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${export_ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_mixed_hydraulic> ${MPIEXEC_POSTFLAGS} -mixed_case profile_export
                -export_prefix ${profile_export_dir}/mpi_${export_ranks}
                -export_reference ${profile_export_dir}/serial)
        set_tests_properties(mixed_profile_export_mpi_${export_ranks} PROPERTIES
            PROCESSORS ${export_ranks}
            TIMEOUT 300
            FIXTURES_REQUIRED profile_export_reference)
    endforeach()
endif()

# The experiment on rectangles, at the same `-mizo_nx` and `-mizo_nz`, so the
# two runs share their 195 temperature vertices and differ in 128 quadrilateral
# cells against 256 triangular ones. This entry asserts what the triangular one
# asserts and nothing about agreement between the two: a comparison of their
# profiles is a run of both with `-mizo_profile_prefix`, not a test.
add_test(NAME mixed_mizo_column_quad
    COMMAND $<TARGET_FILE:test_mixed_hydraulic> -mixed_case mizo_column -mizo_final_hours 12
        -mizo_cell_type quadrilateral)
set_tests_properties(mixed_mizo_column_quad PROPERTIES TIMEOUT 900)
