# P2-D isolated C++ crust module.
#
# Standalone on purpose: it has no DuMuX or Dune dependency, so it can be
# built and checked long before the pinned stack in
# model/p1_r2_conservative_richards/provenance/source_manifest.json is
# available. `pratcrust.hh` is header-only; this target exists only to run the
# reference-vector check.
cmake_minimum_required(VERSION 3.13)
project(pratcrust CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(verify_vectors verify_vectors.cc)
target_include_directories(verify_vectors PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
add_test(NAME reference_vectors
         COMMAND verify_vectors ${CMAKE_CURRENT_SOURCE_DIR}/reference_vectors.txt)
