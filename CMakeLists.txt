cmake_minimum_required(VERSION 3.15)
project(Sturm LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(sturm STATIC sturm.cpp)
# Export the project root as an include directory so tests and consumers can find headers
target_include_directories(sturm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(sturm_app main.cpp)
target_link_libraries(sturm_app PRIVATE sturm)

enable_testing()
add_subdirectory(tests)
