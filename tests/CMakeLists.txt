add_executable(test_sturm test_sturm.cpp)

target_link_libraries(test_sturm PRIVATE sturm)

include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.10
)
FetchContent_MakeAvailable(catch2)

target_include_directories(test_sturm PRIVATE ${catch2_SOURCE_DIR}/single_include)

enable_testing()
add_test(NAME sturm_tests COMMAND test_sturm)
