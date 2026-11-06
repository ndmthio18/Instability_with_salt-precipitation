cmake_minimum_required(VERSION 3.13)
project(p1compositional CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE ${PROJECT_NAME} dune-common dune-common_DIR ${PROJECT_BINARY_DIR})
endif()

find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})
include(DuneMacros)
dune_project()
dune_enable_all_packages()
add_subdirectory(appl)
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
