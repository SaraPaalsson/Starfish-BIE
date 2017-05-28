cmake_minimum_required(VERSION 2.8)
project(MEX)

if(CMAKE_VERSION VERSION_LESS "3.3") 
  message(WARNING "CMake version less than 3.3\n"
    "using local FindMatlab module extracted from v3.3")
  # Where we keep FindMatlab module
  set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake)
endif()

find_package(Matlab REQUIRED MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  message(FATAL_ERROR "MATLAB not found.")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")


# Assuming gcc for now
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -O3  -fopenmp -msse3 -falign-loops=16")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -fopenmp -falign-loops=16 -msse4.1")


include_directories(
  ${Matlab_INCLUDE_DIRS}
)

# Add modules with MEX to be built
add_subdirectory("${PROJECT_SOURCE_DIR}/mex")

# Enable running Matlab test suite through ctest
enable_testing()
matlab_add_unit_test(
  NAME matlab_test_suite
  UNITTEST_FILE "${PROJECT_SOURCE_DIR}/run_unit_tests.m"
  UNITTEST_PRECOMMAND "cd ${PROJECT_SOURCE_DIR}"
  NO_UNITTEST_FRAMEWORK
)
