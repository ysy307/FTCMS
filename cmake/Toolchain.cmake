# =========================================================================
# src/cmake/BuildSettings.cmake
# =========================================================================

# 1. Determine compiler-specific paths and package names
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(COMPILER_DIR "intel")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set(JSON_FORTRAN_PKG "jsonfortran-intelllvm")
    else()
        set(JSON_FORTRAN_PKG "jsonfortran-intel")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(COMPILER_DIR "gnu")
    set(JSON_FORTRAN_PKG "jsonfortran-gnu")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(COMPILER_DIR "nvidia")
    set(JSON_FORTRAN_PKG "jsonfortran-nvhpc")
else()
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# 2. Add local third-party directory to search path
list(PREPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}")

# 3. Find required packages
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

# 4. PETSc, found once for every library. Its BLAS and LAPACK (MKL in the
#    Intel build) are the ones fortran_stdlib links as well.
find_package(PkgConfig REQUIRED)
set(ENV{PKG_CONFIG_PATH} "${PROJECT_SOURCE_DIR}/third_party/.local/${COMPILER_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET GLOBAL PETSc)
add_library(PETSc::PETSc ALIAS PkgConfig::PETSC)
# The engine hands real64 and int32 values to PETSc without conversion.
file(READ "${PETSC_PREFIX}/include/petscconf.h" PETSC_CONF)
if(NOT PETSC_CONF MATCHES "#define PETSC_USE_REAL_DOUBLE 1" OR PETSC_CONF MATCHES "#define PETSC_USE_(COMPLEX|64BIT_INDICES) 1")
    message(FATAL_ERROR "FTCMS needs a PETSc with real double scalars and 32-bit indices: ${PETSC_PREFIX}")
endif()

# 5. BLAS/LAPACK: the ones PETSc was configured with, for fortran_stdlib and
#    the ETKF's dense eigensolver. No second BLAS is linked beside PETSc's.
file(STRINGS "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" PETSC_BLASLAPACK REGEX "^BLASLAPACK_LIB =")
string(REGEX REPLACE "^BLASLAPACK_LIB = *" "" PETSC_BLASLAPACK "${PETSC_BLASLAPACK}")
separate_arguments(PETSC_BLASLAPACK UNIX_COMMAND "${PETSC_BLASLAPACK}")

if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    target_link_libraries(BLAS::BLAS INTERFACE ${PETSC_BLASLAPACK})
endif()

if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    target_link_libraries(LAPACK::LAPACK INTERFACE BLAS::BLAS)
endif()

# 6. Find third-party libraries
find_package(fortran_stdlib REQUIRED)
find_package(${JSON_FORTRAN_PKG} REQUIRED)
find_package(IAPWS REQUIRED)
if(BUILD_TESTING)
    find_package(test-drive REQUIRED)
endif()

# =========================================================================
# Function: enable_build_flags
# =========================================================================
function(enable_build_flags target)
    get_target_property(TARGET_TYPE ${target} TYPE)
    set(KEYWORD PUBLIC)
    if(TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
        set(KEYWORD INTERFACE)
    endif()

    if(NOT KEYWORD STREQUAL "INTERFACE")
        target_compile_options(${target} ${KEYWORD}
            # Fortran: Intel / IntelLLVM on Windows
            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel,IntelLLVM>,$<PLATFORM_ID:Windows>>:
            /stand:f18 /fpp /fpscomp:logicals /extend-source /traceback
            $<$<CONFIG:Release>:/O2>
            $<$<CONFIG:Debug>:/Od /debug:full /check:all /init:snan /init:arrays /warn:all /warn:errors /implicitnone>
            >
            # Fortran: Intel / IntelLLVM on Unix
            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel,IntelLLVM>,$<NOT:$<PLATFORM_ID:Windows>>>:
            -stand f18 -fpp -fpscomp logicals -extend-source -g -traceback
            $<$<CONFIG:Release>:-O2>
            $<$<CONFIG:Debug>:-O0 -check all -init=snan -init=arrays -warn all -warn errors -implicitnone -fstack-protector-all>
            >
            # Fortran: GNU
            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:
            -std=f2018 -cpp
            $<$<CONFIG:Release>:-O3 -march=native>
            $<$<CONFIG:Debug>:-g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-9999999
            -Wall -Wextra
            -Wno-maybe-uninitialized -Wno-uninitialized -Wno-c-binding-type -Wno-surprising
            -Wno-unused-dummy-argument -Wno-compare-reals -Wno-unused-function -Wno-unused-value>
            >
            # Fortran: NVHPC
            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:NVHPC,PGI>>:
            -Mpreprocess
            $<$<CONFIG:Release>:-fast>
            $<$<CONFIG:Debug>:-g -O0 -Mbounds -Mchkptr -traceback -Ktrap=fp>
            >

            # ---------------------------------------------------------
            # C Compile Options
            # ---------------------------------------------------------
            # C: Intel / IntelLLVM
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:Intel,IntelLLVM>>:
            $<$<CONFIG:Release>:-O3 -g -traceback
            $<$<NOT:$<BOOL:${FTCMS_PORTABLE_BUILD}>>:-xHost>>
            $<$<CONFIG:Debug>:-O0 -g -traceback -Wall -Wextra -fstack-protector-all -ftrapv>
            >
            # C: GNU
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:
            $<$<CONFIG:Release>:-O3
            $<$<NOT:$<BOOL:${FTCMS_PORTABLE_BUILD}>>:-march=native>>
            $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra -pedantic -fstack-protector-all -ftrapv>
            >
            # C: NVHPC
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:NVHPC,PGI>>:
            $<$<CONFIG:Release>:-fast>
            $<$<CONFIG:Debug>:-g -O0 -Mbounds -Mchkptr -traceback -Ktrap=fp>
            >

            # ---------------------------------------------------------
            # CXX Compile Options
            # ---------------------------------------------------------
            # CXX: Intel / IntelLLVM
            $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Intel,IntelLLVM>>:
            $<$<CONFIG:Release>:-O3 -g -traceback
            $<$<NOT:$<BOOL:${FTCMS_PORTABLE_BUILD}>>:-xHost>>
            $<$<CONFIG:Debug>:-O0 -g -traceback -Wall -Wextra -fstack-protector-all -ftrapv>
            >
            # CXX: GNU
            $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:
            $<$<CONFIG:Release>:-O3
            $<$<NOT:$<BOOL:${FTCMS_PORTABLE_BUILD}>>:-march=native>>
            $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra -pedantic -fstack-protector-all -ftrapv>
            >
            # CXX: NVHPC
            $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:NVHPC,PGI>>:
            $<$<CONFIG:Release>:-fast>
            $<$<CONFIG:Debug>:-g -O0 -Mbounds -Mchkptr -traceback -Ktrap=fp>
            >
        )
    endif()

    target_compile_definitions(${target} ${KEYWORD}
        $<$<OR:$<CONFIG:Debug>,$<BOOL:${ENABLE_DEBUG}>>:USE_DEBUG>
    )

    if(ENABLE_MPI)
        target_compile_definitions(${target} ${KEYWORD} _MPI)
        target_link_libraries(${target} ${KEYWORD} MPI::MPI_Fortran)
    endif()
    target_link_libraries(${target} ${KEYWORD} PETSc::PETSc LAPACK::LAPACK)
endfunction()

# =========================================================================
# Function: enable_thirdparty
# =========================================================================
function(enable_thirdparty target)
    target_include_directories(${target} PUBLIC
        $<TARGET_PROPERTY:fortran_stdlib::fortran_stdlib,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:${JSON_FORTRAN_PKG}::jsonfortran-static,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:IAPWS::IAPWS,INTERFACE_INCLUDE_DIRECTORIES>
    )

    target_link_libraries(${target} PUBLIC
        fortran_stdlib::fortran_stdlib
        ${JSON_FORTRAN_PKG}::jsonfortran-static
        IAPWS::IAPWS
    )
endfunction()
