# Build the Ice library
# The list of source files that should be compiled for this component.

set srcList = ( oasis_interface )

set lib_name = "lib_oasis"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name

$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean -I$executable:h:h/lib_ice  -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/mct -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/psmile.MPI1 -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/scrip -I${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${IOW_ESM_COMPILE_MODE}/build/lib/mctdir" $srcList $lib_include_dirs

${IOW_ESM_MAKE}

if( $status ) then
    echo "Make failed to create $lib_name"
    exit 1
endif
