cmake_minimum_required(VERSION 3.12)

set(LIB_NAME "FixedPoint")

project(${LIB_NAME})

add_library(${LIB_NAME} STATIC
	../../../lib/${LIB_NAME}/src/fixedpoint.c
	../../../lib/${LIB_NAME}/src/tables.c)

target_include_directories(${LIB_NAME} PRIVATE ../../../lib/${LIB_NAME}/inc)