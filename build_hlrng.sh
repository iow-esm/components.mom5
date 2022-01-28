#!/bin/bash

debug=${1:-"release"}
rebuild=${2:-"fast"}

# LOAD REQUIRED MODULES
module load intel/18.0.5
module load impi/2018.5
module load netcdf/intel/4.7.3

# GET IOW ESM ROOT PATH
export IOW_ESM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/../.."

# SET SYSTEM-SPECIFIC COMPILER OPTIONS AND PATHS
# compile mode: "PRODUCTION" or "DEBUG"
if [ $debug == "debug" ]; then
	export IOW_ESM_COMPILE_MODE="DEBUG"
elif [ $debug == "release" ]; then
	export IOW_ESM_COMPILE_MODE="PRODUCTION"
else
	echo "Compile mode is not specified correctly. Use debug or release"
	exit;
fi

# include paths
export IOW_ESM_NETCDF_INCLUDE="/sw/dataformats/netcdf/intel.18/4.7.3/skl/include"
export IOW_ESM_NETCDF_LIBRARY="/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib"
# executables
export IOW_ESM_MAKE="/usr/bin/make"
export IOW_ESM_FC="mpiifort"
export IOW_ESM_CC="mpiicc"
export IOW_ESM_LD="mpiifort"


if [ $debug == "debug" ]; then
	# compiler flags
	export IOW_ESM_CPPDEFS="-DOASIS_IOW_ESM -DIOW_ESM_DEBUG"

	export IOW_ESM_FFLAGS="-O0 -r8 -g -traceback -check all -fp-model fast=2 -xHost -I${IOW_ESM_NETCDF_INCLUDE}/"
	export IOW_ESM_CFLAGS="-O0 -r8 -g -traceback -fp-model fast=2 -xHost -I${IOW_ESM_NETCDF_INCLUDE}/"
	export IOW_ESM_LDFLAGS="-g -traceback -L${IOW_ESM_NETCDF_LIBRARY} -lnetcdf -lnetcdff -Wl,-rpath,${IOW_ESM_NETCDF_LIBRARY}"
else
	# compiler flags
	export IOW_ESM_CPPDEFS="-DOASIS_IOW_ESM"
	
	export IOW_ESM_FFLAGS="-O3 -r8 -no-prec-div -fp-model fast=2 -xHost -I${IOW_ESM_NETCDF_INCLUDE}/"
	export IOW_ESM_CFLAGS="-O3 -r8 -no-prec-div -fp-model fast=2 -xHost -I${IOW_ESM_NETCDF_INCLUDE}/"
	export IOW_ESM_LDFLAGS="-L${IOW_ESM_NETCDF_LIBRARY} -lnetcdf -lnetcdff -Wl,-rpath,${IOW_ESM_NETCDF_LIBRARY}"
fi

# MAKE CLEAN
if [ $rebuild == "rebuild" ]; then
	rm -r ${IOW_ESM_ROOT}/components/MOM5/exec/IOW_ESM_${IOW_ESM_COMPILE_MODE}
fi

# RUN BUILD COMMAND
cd ${IOW_ESM_ROOT}/components/MOM5/exp
./MOM_compile.csh
cd ${IOW_ESM_ROOT}/components/MOM5

