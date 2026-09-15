cmake_minimum_required(VERSION 3.20)

# ====================================================================
# 1. Project Setup
# ====================================================================
project(FTCMS LANGUAGES Fortran C CXX)

option(BUILD_TESTING "Build the testing tree" OFF)
option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(ENABLE_MKL "Enable MKL support" ON)
option(ENABLE_MPI "Enable MPI support" ON)
option(ENABLE_PETSC "Enable the PETSc KSP linear solver back-end" ON)
option(FTCMS_PORTABLE_BUILD "Avoid host-specific CPU instructions" OFF)
option(FTCMS_BUNDLE_RUNTIME "Bundle non-system runtime libraries when installing" OFF)

# 出力先ディレクトリを統一
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/.mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# カスタムモジュールの読み込み
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Toolchain OPTIONAL)
include(fypp OPTIONAL)

# ====================================================================
# 3. Source Tree (Libraries)
# ====================================================================
# 各サブディレクトリ内の CMakeLists.txt で target_link_libraries() を用い，
# 個別に依存関係を定義することで，CMakeが自動的にビルド順序を解決する．
add_subdirectory(src/core)
add_subdirectory(src/numerical)
add_subdirectory(src/control)
add_subdirectory(src/domain)
add_subdirectory(src/io)
add_subdirectory(src/physics)
add_subdirectory(src/app/ftcms)
add_subdirectory(src/app/plexfem)
add_subdirectory(src/app/mixed_hydraulic)
# ====================================================================
# 4. Executable
# ====================================================================
# プロジェクト直下の app ディレクトリにあるエントリーポイント
add_executable(FTCMS app/ftcms.f90)
set_target_properties(FTCMS PROPERTIES LINKER_LANGUAGE Fortran)

# 実行ファイルに必要な最上位ライブラリターゲットをリンクする（名称は仮）
# 必要なフラグ（インクルードパスやOpenMP等）もターゲット経由で自動伝搬する
target_link_libraries(FTCMS PRIVATE SOLVER::FTCMS)

# The PetscDS/TS path is a separate executable so that the existing one is not
# rebuilt against it: the two share the input layer and the mesh, nothing else.
add_executable(FTCMS_PLEXFEM app/ftcms_plexfem.f90)
set_target_properties(FTCMS_PLEXFEM PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(FTCMS_PLEXFEM PRIVATE SOLVER::PLEXFEM)

include(GNUInstallDirs)

if(FTCMS_BUNDLE_RUNTIME)
    set(ftcms_runtime_directories
        ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
        ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
        ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}
    )
    if(WIN32)
        list(APPEND ftcms_runtime_directories "$ENV{PATH}")
    endif()
    list(REMOVE_DUPLICATES ftcms_runtime_directories)

    install(TARGETS FTCMS
        RUNTIME_DEPENDENCY_SET ftcms_runtime
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    install(RUNTIME_DEPENDENCY_SET ftcms_runtime
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        DIRECTORIES ${ftcms_runtime_directories}
        PRE_EXCLUDE_REGEXES
            "^api-ms-"
            "^ext-ms-"
            "^ld-linux"
            "^libc\\.so"
            "^libdl\\.so"
            "^libgcc_s\\.so"
            "^libm\\.so"
            "^libpthread\\.so"
            "^librt\\.so"
            "^libstdc\\+\\+\\.so"
        POST_EXCLUDE_REGEXES
            ".*[Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\][Ss][Yy][Ss][Tt][Ee][Mm]32[/\\\\].*"
    )

    if(UNIX)
        set_target_properties(FTCMS PROPERTIES
            BUILD_WITH_INSTALL_RPATH ON
            INSTALL_RPATH "$ORIGIN"
        )
    endif()
else()
    install(TARGETS FTCMS
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

install(FILES cmake/package/README.md DESTINATION .)
install(FILES LICENSE DESTINATION .)

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    target_link_libraries(FTCMS PRIVATE MPI::MPI_Fortran)
    add_custom_target(run
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:FTCMS>
        DEPENDS FTCMS
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
else()
    add_custom_target(run
        COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=4 $<TARGET_FILE:FTCMS>
        DEPENDS FTCMS
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()

if(FTCMS_BUNDLE_RUNTIME AND UNIX AND NOT APPLE)
    # Keep $ORIGIN effective for indirect Intel runtime dependencies.
    target_link_libraries(FTCMS PRIVATE "-Wl,--disable-new-dtags")
endif()

# ====================================================================
# 5. Testing
# ====================================================================
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

set(CPACK_PACKAGE_NAME "FTCMS")
set(CPACK_PACKAGE_FILE_NAME "FTCMS-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
if(WIN32)
    set(CPACK_GENERATOR "ZIP")
else()
    set(CPACK_GENERATOR "TGZ")
endif()
include(CPack)
