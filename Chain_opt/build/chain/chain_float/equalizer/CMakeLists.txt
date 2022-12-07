cmake_minimum_required(VERSION 3.12)

set(LIB_NAME "equalizer")

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-O3")


project(${LIB_NAME})

add_library(${LIB_NAME} 
		../../../../chain_float/${LIB_NAME}/src/equalizer_control.c
		../../../../chain_float/${LIB_NAME}/src/equalizer_process.c)

target_include_directories(${LIB_NAME} PRIVATE
		../../../../chain_float/${LIB_NAME}/inc
		../../../../lib/FixedPoint/inc
)

target_link_libraries(${LIB_NAME} PRIVATE
	FixedPoint)
