cmake_minimum_required(VERSION 3.25)

project(FTCMSDependencies LANGUAGES NONE)

include(ExternalProject)

set(FTCMS_DEPENDENCY_PREFIX
    "${CMAKE_CURRENT_BINARY_DIR}/install"
    CACHE PATH
    "Installation prefix for FTCMS dependencies"
)

foreach(language Fortran C CXX)
    if(NOT CMAKE_${language}_COMPILER)
        message(FATAL_ERROR "CMAKE_${language}_COMPILER must be specified")
    endif()
endforeach()

set(common_cmake_args
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_INSTALL_PREFIX:PATH=${FTCMS_DEPENDENCY_PREFIX}
    -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
)

ExternalProject_Add(fortran_stdlib
    GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
    GIT_TAG v0.7.0
    GIT_SHALLOW TRUE
    UPDATE_DISCONNECTED TRUE
    CMAKE_ARGS
        ${common_cmake_args}
        -DFIND_BLAS:BOOL=OFF
        -DCMAKE_MAXIMUM_RANK:STRING=4
)

ExternalProject_Add(json_fortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG 9.0.3
    GIT_SHALLOW TRUE
    UPDATE_DISCONNECTED TRUE
    CMAKE_ARGS
        ${common_cmake_args}
        -DCMAKE_Fortran_MODULE_DIRECTORY:PATH=<BINARY_DIR>/include
        -DJSONFORTRAN_ENABLE_TESTS:BOOL=OFF
        -DJSONFORTRAN_STATIC_LIBRARY_ONLY:BOOL=ON
)

ExternalProject_Add(iapws
    GIT_REPOSITORY https://github.com/ysy307/IAPWS.git
    GIT_TAG 9c75bae6243de938a5bd2bdd1b90b4e850832ddd
    UPDATE_DISCONNECTED TRUE
    CMAKE_ARGS ${common_cmake_args}
)
