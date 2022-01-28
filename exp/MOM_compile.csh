#!/bin/csh -f
# Minimal compile script for fully coupled model CM2M experiments

#module unload intel.compiler mpt
#module load intel.compiler/12.0.4
#module load mpt/1.26

set echo
set platform      = IOW_ESM_${IOW_ESM_COMPILE_MODE}  # A unique identifier for your platform
                                  # This corresponds to the mkmf templates in $root/bin dir.
set type          = MOM_SIS       # Type of the experiment
set help = 0

set coupled       = 1             # 1: in case of COSMO_clm coupled to MOM_SIS

set argv = (`getopt -u -o h -l type: -l platform:  -l help  --  $*`)
while ("$argv[1]" != "--")
    switch ($argv[1])
        case --type:
                set type = $argv[2]; shift argv; breaksw
        case --platform:
                set platform = $argv[2]; shift argv; breaksw
        case --help:
                set help = 1;  breaksw
        case -h:
                set help = 1;  breaksw
    endsw
    shift argv
end
shift argv
if ( $help ) then
    echo "The optional arguments are:"
    echo "--type       followed by the type of the experiment, currently one of the following:"
    echo "             MOM_solo : solo ocean model"
    echo "             MOM_SIS  : ocean-seaice model"
    echo "             CM2M     : ocean-seaice-land-atmosphere coupled climate model"
    echo "             ESM2M    : ocean-seaice-land-atmosphere coupled climate model with biogeochemistry, EarthSystemModel"
    echo "             ICCM     : ocean-seaice-land-atmosphere coupled model"
    echo "             EBM      : ocean-seaice-land-atmosphere coupled model with energy balance atmosphere"
    echo
    echo "--platform   followed by the platform name that has a corresponding environ file in the ../bin dir, default is ncrc.intel"
    echo
    echo
    exit 0
endif

#
# User does not need to change anything below!
#
set root          = $cwd:h                            # The directory you created when you checkout
set code_dir      = $root/src                         # source code directory
set executable    = $root/exec/$platform/$type/fms_$type.x      # executable created after compilation
set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set mkmfTemplate  = $root/bin/mkmf.template.$platform # path to template for your platform
set mkmf          = $root/bin/mkmf                    # path to executable mkmf
set cppDefs  = ( "-Duse_netCDF -Duse_netCDF4 -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS ${IOW_ESM_CPPDEFS}")
#On Altrix systems you may include "-Duse_shared_pointers -Duse_SGI_GSM" in cppDefs for perfomance.
#These are included in the GFDL configuration of the model.
  
set static        = 0              # 1 if you want static memory allocation, 0 for dynamic
if($static) then
  set executable = $root/exec/$platform/${type}_static/fms_$type.x 
  set cppDefs = "$cppDefs -DMOM_STATIC_ARRAYS -DNI_=360 -DNJ_=200 -DNK_=50 -DNI_LOCAL_=60 -DNJ_LOCAL_=50"
endif

if ( $type == EBM ) then
  set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DLAND_BND_TRACERS -DOVERLOAD_C8 -DOVERLOAD_C4 -DOVERLOAD_R4" )
endif

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules


#
## compile mppnccombine.c, needed only if $npes > 1
#  if ( ! -f $mppnccombine ) then
#    cc -O -o $mppnccombine -I/usr/local/include -L/usr/local/lib $code_dir/postprocessing/mppnccombine/mppnccombine.c -lnetcdf
#  endif

set mkmf_lib = "$mkmf -f -m Makefile -a $code_dir -t $mkmfTemplate"
set lib_include_dirs = "$root/include $code_dir/shared/include $code_dir/shared/mpp/include -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib"

source ./FMS_compile.csh
if ( $status ) exit $status

cd $root/exp
source ./ocean_compile.csh
if ( $status ) exit $status

if( $type != MOM_solo) then
    cd $root/exp
    source ./ice_compile.csh
    if ( $status ) exit $status
endif
if( $type == MOM_SIS) then
    cd $root/exp
    source ./land_null_compile.csh
    if ( $status ) exit $status

    cd $root/exp
    source ./atmos_null_compile.csh
    if ( $status ) exit $status

    if ($coupled == 1) then
      cd $root/exp
      source ./oasis_compile.csh
      if ( $status ) exit $status
    endif
endif
if( $type == EBM) then
    cd $root/exp
    source ./atmos_ebm_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ESM2M | $type == ICCM ) then
    cd $root/exp
    source ./atmos_phys_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ESM2M ) then
    cd $root/exp
    source ./atmos_fv_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ICCM | $type == EBM ) then
    cd $root/exp
    source ./land_lad_compile.csh
    if ( $status ) exit $status
endif
if( $type == ESM2M ) then
    cd $root/exp
    source ./land_lad2_compile.csh
    if ( $status ) exit $status
endif
if( $type == ICCM ) then
    cd $root/exp
    source ./atmos_bg_compile.csh
    if ( $status ) exit $status
endif

# Build the executable
set mkmf_exec = "$mkmf -f -m Makefile -a $code_dir -t $mkmfTemplate -p $executable:t"
mkdir -p $executable:h
cd $executable:h
if( $type == MOM_solo ) then
    set srcList = ( mom5/drivers )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == MOM_SIS ) then
    set srcList = ( coupler )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice  -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/mct -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/psmile.MPI1 -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/scrip -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/mctdir -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib -I$executable:h:h/lib_oasis -I$executable:h:h/lib_atmos_null -I$executable:h:h/lib_land_null"
    set libs = "$executable:h:h/lib_oasis/lib_oasis.a $executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_null/lib_atmos_null.a ${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib/libpsmile.MPI1.a ${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib/libmct.a ${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib/libmpeu.a ${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/lib/libscrip.a $executable:h:h/lib_land_null/lib_land_null.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == EBM ) then
    set srcList = ( coupler )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_ebm  -I$executable:h:h/lib_land_lad"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_ebm/lib_atmos_ebm.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == CM2M ) then
    set srcList = ( coupler )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_fv -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_fv/lib_atmos_fv.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == ESM2M ) then
    set srcList = ( coupler )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_fv -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad2"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_fv/lib_atmos_fv.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad2/lib_land_lad2.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == ICCM ) then
    set srcList = ( coupler )
    set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_bg -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad" 
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_bg/lib_atmos_bg.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
endif
$mkmf_exec -o "$includes" -l "$libs" -c "$cppDefs" $srcList
${IOW_ESM_MAKE}
if( $status ) then
    echo "Make failed to create the $type executable"
    exit 1
endif    

exit
