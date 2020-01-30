cmake_minimum_required(VERSION 3.1)

project(opinion_socp)

set(CMAKE_CXX_STANDARD 11)

if(UNIX)
	set(ENABLE_MULTITHREADING_COMPILE_FLAGS "-pthread")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ENABLE_MULTITHREADING_COMPILE_FLAGS}")
endif(UNIX)

find_package(SOCP REQUIRED)

add_subdirectory(src)
add_subdirectory(app)
