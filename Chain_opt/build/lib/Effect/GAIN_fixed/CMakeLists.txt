cmake_minimum_required(VERSION 3.12)

set(LIB_NAME "GAIN_fixed")

project(${LIB_NAME})

add_library(${LIB_NAME} STATIC
	../../../../lib/Effect/src/${LIB_NAME}/effect_control.c
	../../../../lib/Effect/src/${LIB_NAME}/effect_process.c)

target_include_directories(${LIB_NAME} PRIVATE  ../../../../lib/Effect/inc
						../../../../lib/FixedPoint/inc)