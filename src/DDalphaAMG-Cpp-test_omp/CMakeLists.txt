# main CMake configuration file for DDalphaAMG

# ---------------------------------------------------------

# basic CMake settings

cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_sTANDARD_REQUIRED ON)

project(DDalphaAMG VERSION 1.0.0)

#Gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# sources files

file(GLOB_RECURSE SRC_FILES lib/*.cpp)
file(GLOB_RECURSE TST_FILES tests/*.cpp)
add_executable(ddalphaamg ${SRC_FILES} ${TST_FILES})

# compilation options

target_compile_options(ddalphaamg PRIVATE
			-Wall
			-pedantic
			-O3
			-ffast-math
			-fopenmp
			-DOPENMP
			-Wfatal-errors
)

add_compile_definitions(BLEN=2)

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  target_compile_options(ddalphaamg PRIVATE
			-DDEBUG
			-DPARAMOUTPUT
			-DTRACK_RES
			-DFGMRES_RESTEST
			#-DPROFILING
			-DCOARSE_RES
			-DSCHWARZ_RES
			-DTESTVECTOR_ANALYSIS
  )
elseif(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
  target_compile_options(ddalphaamg PRIVATE
			-DPARAMOUTPUT
			-DTRACK_RES
  )
else()
  message( FATAL_ERROR "Unknown (or unaccepted) build type, use -DCMAKE_BUILD_TYPE=[Debug|RelWithDebInfo]")
endif()

# include files

target_include_directories(ddalphaamg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# MPI and LAPACK linkage

find_package(MPI REQUIRED)
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
find_package(LAPACK REQUIRED)
find_package(OpenMP)

target_link_libraries(ddalphaamg PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX ${LAPACK_LIBRARIES} gtest_main)

include(GoogleTest)

# ------------------------------

# TODOs :

#	1. manually include LIME. Not only as a library, but automate its installation (qmp and qio)
# 	2. use of Google Tests
#	3. enable compilation of tests written by Simone ---> do this via Google Tests already?
#	4. enable <doxygen> target
#	5. compile DDalphaAMG as a library, and call this library from the tests. Do this as the tmLQCD people would
