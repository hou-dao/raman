cmake_minimum_required (VERSION 2.8)
project (HEOM)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

set (CMAKE_CXX_COMPILER g++)
set (CMAKE_CXX_FLAGS "-std=c++11 -g -Wall -O2 -ffast-math")

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/deom ${PROJECT_SOURCE_DIR}/src/resp)

link_directories (${PROJECT_SOURCE_DIR}/lib)

set (src_deom)
set (src_resp)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/deom)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/resp)

file (GLOB apps_src ${PROJECT_SOURCE_DIR}/src/apps/*.cpp)

foreach (p ${apps_src})
    get_filename_component(p_bin ${p} NAME_WE)
    add_executable(${p_bin} ${p} ${src_deom} ${src_resp})
    target_link_libraries(${p_bin} openblas lapack json11 armadillo)
endforeach ()

enable_language(Fortran)

include(ExternalProject)

find_package(LAPACK)

if(${LAPACK_FOUND})
else()
ExternalProject_Add (lapack
    PREFIX ${PROJECT_BINARY_DIR}
    DOWNLOAD_DIR download
    SOURCE_DIR source
    URL http://www.netlib.org/lapack/lapack-3.6.0.tgz
    INSTALL_DIR ${PROJECT_SOURCE_DIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)
endif()

if (EXISTS "${PROJECT_SOURCE_DIR}/include/json11.hpp")
else()
  ExternalProject_Add(json11
    GIT_REPOSITORY https://github.com/dropbox/json11.git
    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
    INSTALL_DIR ${PROJECT_SOURCE_DIR}
    INSTALL_COMMAND mkdir -p ${PROJECT_SOURCE_DIR}/include && cp <SOURCE_DIR>/json11.hpp <INSTALL_DIR>/include && cp <BINARY_DIR>/libjson11.a <INSTALL_DIR>/lib
  )
endif()

if (EXISTS "${PROJECT_SOURCE_DIR}/lib/libopenblas.a")
else()
  ExternalProject_Add (openblas
      GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
      PREFIX ${PROJECT_BINARY_DIR}/thirdparty
      CONFIGURE_COMMAND ""
      BUILD_IN_SOURCE 1
      BUILD_COMMAND cd <SOURCE_DIR> && $(MAKE)
      INSTALL_DIR ${PROJECT_SOURCE_DIR}
      INSTALL_COMMAND cd <SOURCE_DIR> && $(MAKE) install PREFIX=<INSTALL_DIR>
  )
endif()

if (EXISTS "${PROJECT_SOURCE_DIR}/include/armadillo")
else()
  ExternalProject_Add(armadillo
      URL http://sourceforge.net/projects/arma/files/armadillo-6.400.3.tar.gz
      PREFIX ${PROJECT_BINARY_DIR}/thirdparty
      INSTALL_DIR ${PROJECT_SOURCE_DIR}
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  )
endif()
