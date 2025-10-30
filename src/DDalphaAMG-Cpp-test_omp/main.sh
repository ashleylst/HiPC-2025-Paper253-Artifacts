#!/bin/bash

# this script automates cmake-ing, compiling and running

# run for help : ./main --help

# run as (we assume this is done correctly) : ./main.sh --build=1 --compile=1 --run=1

# ------------------------------

# displaying help, if requested

for var in "$@"
do
    if [[ "$var" == "--help" ]]; then
        echo "The possible flags for this script are:"
        echo "	--build=0/1 : to build with CMake, it creates ../build/ and builds there. The build mode can be Debug or RelWithDebInfo, and to set this use --buildmode"
        echo "	--compile=0/1 : assuming ../build/, it compile with make"
        echo "	--run=0/1 : it runs DDalphaAMG. The default input file is sample.ini, but another can be specified via --infile=NAME"
        echo "	--clean=0/1 : applies <make clean> inside of ../build/"
        echo "	--purge=0/1 : forcefully removes ../build/"
        echo "	--infile=FILENAME : see --run above for more information"
        echo "	--buildmode=Debug/RelWithDebInfo : the mode we want to use when building with CMake"
        exit 0
    fi
done

# ------------------------------

# SOME PARAMETERS THAT CAN BE CHANGED

# possible values for build type :
#	Debug: Usually a classic debug build including debugging information, no optimization etc.
#	Release: Your typical release build with no debugging information and full optimization.
#	RelWithDebInfo: Same as Release, but with debugging information.
#	MinSizeRel: A special Release build optimized for size.
BUILD_TYPE=RelWithDebInfo

NR_CORES_TO_COMPILE=32

# ------------------------------

# SOME PARAMETERS THAT SHOULDN'T BE CHANGED

SRC_DIR=`pwd`

BUILD="0"
COMPILE="0"
RUN="0"
CLEAN="0"
PURGE="0"
INFILE="sample.ini"
for var in "$@"
do
    P="$var"
    arrP=(${P//=/ })
    if [[ "${arrP[0]}" == "--build" ]]; then
        BUILD="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--compile" ]]; then
        COMPILE="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--run" ]]; then
        RUN="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--clean" ]]; then
        CLEAN="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--purge" ]]; then
        PURGE="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--infile" ]]; then
        INFILE="${arrP[1]}"
    fi
    if [[ "${arrP[0]}" == "--buildmode" ]]; then
        BUILD_TYPE="${arrP[1]}"
    fi
done

# ------------------------------

# BUILDING, COMPILING AND/OR RUNNING

if [[ "$BUILD" == "1" ]]; then
    DIR="../build/"
    if ! [ -d "$DIR" ]; then
      # if the directory ../build/ doesn't exist, create it
      mkdir ../build/
    fi
    cd ../build/

    MPIC_LOCATION=`which mpicc`
    MPICXX_LOCATION=`which mpicxx`
    cmake -DMPI_C_COMPILER=$MPIC_LOCATION -DMPI_CXX_COMPILER=$MPICXX_LOCATION -DCMAKE_BUILD_TYPE=$BUILD_TYPE $SRC_DIR

    cd $SRC_DIR
fi

if [[ "$COMPILE" == "1" ]]; then
    cd ../build/
    make -j $NR_CORES_TO_COMPILE
    cd $SRC_DIR
fi

if [[ "$RUN" == "1" ]]; then
    . run -i $INFILE
fi

if [[ "$CLEAN" == "1" ]]; then
    cd ../build/
    make clean
    cd $SRC_DIR
fi

if [[ "$PURGE" == "1" ]]; then
    rm -R -f ../build/
fi
