cmake_minimum_required(VERSION 3.13)
project(pratcoupling CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(verify_vectors ../../p2_d_coupling/cpp/verify_vectors.cc)
target_include_directories(verify_vectors PRIVATE ../../p2_d_coupling/cpp)

add_executable(verify_lifecycle verify_lifecycle.cc)
target_include_directories(verify_lifecycle PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
add_test(NAME reference_vectors
         COMMAND verify_vectors ${CMAKE_CURRENT_SOURCE_DIR}/../../p2_d_coupling/cpp/reference_vectors.txt)
add_test(NAME lifecycle COMMAND verify_lifecycle)
