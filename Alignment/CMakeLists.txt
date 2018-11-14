project(Alignment)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set(SOURCE_FILES main.cpp affinealignobj.cpp alignment.cpp affinealignment.cpp)
add_executable(${PROJECT_NAME} ${SRC_LIST})

