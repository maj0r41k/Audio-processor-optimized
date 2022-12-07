cmake_minimum_required(VERSION 3.12)

project(Solution)

add_subdirectory(lib/FixedPoint)

add_subdirectory(chain/chain_float/compressor_4ch)
add_subdirectory(chain/chain_float/crossover)
add_subdirectory(chain/chain_float/equalizer)
add_subdirectory(chain/chain_float/chain)

add_subdirectory(app)











#add_subdirectory(lib/Effect/CROSSOVER_float)
#add_subdirectory(lib/Effect/EQUALIZER_float)

#add_subdirectory(lib/Effect/EXPANDER_float)
#add_subdirectory(lib/Effect/EXPANDER_fixed)

#add_subdirectory(lib/Effect/COMPRESSOR_float)
#add_subdirectory(lib/Effect/COMPRESSOR_fixed)

#add_subdirectory(lib/Effect/GAIN_fixed)
#add_subdirectory(lib/Effect/GAIN_float)

#add_subdirectory(lib/Effect/FIR_float)
#add_subdirectory(lib/Effect/FIR_fixed)

#add_subdirectory(lib/Effect/IIR_float)
#add_subdirectory(lib/Effect/IIR_fixed)