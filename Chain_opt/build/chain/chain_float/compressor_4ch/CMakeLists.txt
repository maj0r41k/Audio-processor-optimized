cmake_minimum_required(VERSION 3.12)

set(LIB_NAME "compressor_4ch")

project(${LIB_NAME})

add_library(${LIB_NAME} 
		../../../../chain_float/${LIB_NAME}/src/compressor_4ch_control.c
		../../../../chain_float/${LIB_NAME}/src/compressor_4ch_process.c)

target_include_directories(${LIB_NAME} PRIVATE
		../../../../chain_float/${LIB_NAME}/inc
		../../../../lib/FixedPoint/inc
)

target_link_libraries(${LIB_NAME} PRIVATE
	FixedPoint)
