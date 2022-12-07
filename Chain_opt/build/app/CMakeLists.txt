cmake_minimum_required(VERSION 3.12)

set(PR_NAME "Filters")

project(${PR_NAME})

add_executable(${PR_NAME}
	../../src/main.c
        ../../src/cJSON.c
	../../src/file_handler.c)

target_include_directories(${PR_NAME} PRIVATE ../../inc)
target_include_directories(${PR_NAME} PRIVATE ../../lib/FixedPoint/inc)
target_include_directories(${PR_NAME} PRIVATE ../../lib/Effect/inc)
target_include_directories(${PR_NAME} PRIVATE ../../chain_float/chain/inc)


target_link_libraries(${PR_NAME} PRIVATE FixedPoint)
target_link_libraries(${PR_NAME} PRIVATE chain)

SET_TARGET_PROPERTIES(${PR_NAME} PROPERTIES LINK_FLAGS "/PROFILE")

#target_link_libraries(${PR_NAME} PRIVATE EQUALIZER_float)
#target_link_libraries(${PR_NAME} PRIVATE CROSSOVER_float)

#target_link_libraries(${PR_NAME} PRIVATE EXPANDER_fixed)
#target_link_libraries(${PR_NAME} PRIVATE EXPANDER_float)

#target_link_libraries(${PR_NAME} PRIVATE COMPRESSOR_fixed)
#target_link_libraries(${PR_NAME} PRIVATE COMPRESSOR_float)

#target_link_libraries(${PR_NAME} PRIVATE GAIN_fixed)
#target_link_libraries(${PR_NAME} PRIVATE GAIN_float)

#target_link_libraries(${PR_NAME} PRIVATE IIR_fixed)
#target_link_libraries(${PR_NAME} PRIVATE IIR_float)

#target_link_libraries(${PR_NAME} PRIVATE FIR_fixed)
#target_link_libraries(${PR_NAME} PRIVATE FIR_float)

if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()