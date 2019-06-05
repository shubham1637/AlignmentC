project(DIAlign)
include(CTest)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set(CMAKE_CXX_STANDARD 11) # enable C++11 standard
set(SOURCE_FILES affinealignobj.cpp alignment.cpp affinealignment.cpp chromSimMatrix.cpp)

add_executable(${PROJECT_NAME} main.cpp ${SOURCE_FILES})

# Test executables
add_executable(FullTest test.cpp ${SOURCE_FILES})
add_executable(testAlignment test_alignment.cpp ${SOURCE_FILES})
add_test(RunTest FullTest) 
add_test(RunTest2 FullTest)
add_test(RunTest3 FullTest --option 1)
add_test(RuntestAlignment testAlignment)

