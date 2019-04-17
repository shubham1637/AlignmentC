project(Alignment)
include(CTest)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set (CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES affinealignobj.cpp alignment.cpp affinealignment.cpp)

add_executable(${PROJECT_NAME} main.cpp ${SOURCE_FILES})

# Test executables
add_executable(FullTest test.cpp ${SOURCE_FILES})
add_test(RunTest FullTest) 

