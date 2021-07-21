target=$1
debug=${2:-"release"}
fast=${3:-"fast"}

source ../../local_scripts/identify_target.sh $target $debug $fast

# component to build
component="MOM5"

# create component's directory (if not existing)
echo ssh -t "${user_at_dest}" "mkdir -p ${dest_folder}/components/${component}"
ssh -t "${user_at_dest}" "mkdir -p ${dest_folder}/components/${component}"

# deploy the code from local source
echo rsync -r -i -u src ${dest}/components/${component}/.
echo rsync -r -i -u bin ${dest}/components/${component}/.
echo rsync -r -i -u exp ${dest}/components/${component}/.
echo rsync -i -u start_build_${target}.sh ${dest}/components/${component}/
echo rsync -i -u build_${target}.sh ${dest}/components/${component}/
rsync -r -i -u src ${dest}/components/${component}/.
rsync -r -i -u bin ${dest}/components/${component}/.
rsync -r -i -u exp ${dest}/components/${component}/.
rsync -i -u start_build_${target}.sh ${dest}/components/${component}/
rsync -i -u build_${target}.sh ${dest}/components/${component}/

# start the build process
echo ssh -t "${user_at_dest}" "cd ${dest_folder}/components/${component}/; bash -c \"source ~/.bash_profile; source start_build_${target}.sh $debug $fast\""
ssh -t "${user_at_dest}" "cd ${dest_folder}/components/${component}/; bash -c \"source ~/.bash_profile; source start_build_${target}.sh $debug $fast\""

../../local_scripts/tag_build.sh "$target" "$debug" "$fast"
