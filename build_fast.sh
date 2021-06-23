#!/bin/bash
# LOAD REQUIRED MODULES
module load intel/18.0.6
module load impi/2018.5
module load netcdf/intel/4.7.3

# GET IOW ESM ROOT PATH
export IOW_ESM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/../.."

# SET SYSTEM-SPECIFIC COMPILER OPTIONS AND PATHS
# compile mode: "PRODUCTION" or "DEBUG"
export IOW_ESM_COMPILE_MODE="PRODUCTION"
# include paths
export IOW_ESM_NETCDF_INCLUDE="/sw/dataformats/netcdf/intel.18/4.7.3/skl/include"
export IOW_ESM_NETCDF_LIBRARY="/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib"
# executables
export IOW_ESM_MAKE="/usr/bin/make"
export IOW_ESM_FC="mpiifort"
export IOW_ESM_CC="mpiicc"
export IOW_ESM_LD="mpiifort"
# compiler flags
export IOW_ESM_CPPDEFS="-DCOUP_OAS -DOASIS_IOW_ESM"
export IOW_ESM_FFLAGS="-O3 -xCORE-AVX512 -r8 -fp-model precise -align array64byte -I${IOW_ESM_NETCDF_INCLUDE}/"
export IOW_ESM_CFLAGS="-O3 -I${IOW_ESM_NETCDF_INCLUDE}/"
export IOW_ESM_LDFLAGS="-g -traceback -L${IOW_ESM_NETCDF_LIBRARY} -lnetcdf -lnetcdff -Wl,-rpath,${IOW_ESM_NETCDF_LIBRARY}"

# MAKE CLEAN
#rm -r ${IOW_ESM_ROOT}/components/MOM5/exec/IOW_ESM_${IOW_ESM_COMPILE_MODE}

# RUN BUILD COMMAND
cd ${IOW_ESM_ROOT}/components/MOM5/exp
./MOM_compile.csh
cd ${IOW_ESM_ROOT}/components/MOM5

