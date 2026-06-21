cmake_minimum_required(VERSION 3.21)
project(Sturm LANGUAGES CXX VERSION 1.0.0)

option(STURM_BUILD_TESTS "Build the unit tests" ON)
option(STURM_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(STURM_WERROR "Treat compiler warnings as errors" OFF)

add_library(sturm STATIC src/sturm.cpp)
add_library(Sturm::sturm ALIAS sturm)

target_include_directories(sturm PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(sturm PUBLIC cxx_std_20)

set(STURM_WARNINGS
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic;-Wshadow>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>)
target_compile_options(sturm PRIVATE ${STURM_WARNINGS})
if(STURM_WERROR)
  target_compile_options(sturm PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Werror>
    $<$<CXX_COMPILER_ID:MSVC>:/WX>)
endif()

if(STURM_ENABLE_SANITIZERS)
  target_compile_options(sturm PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(sturm PUBLIC -fsanitize=address,undefined)
endif()

add_executable(sturm_app src/main.cpp)
target_link_libraries(sturm_app PRIVATE Sturm::sturm)
target_compile_options(sturm_app PRIVATE ${STURM_WARNINGS})

if(STURM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
