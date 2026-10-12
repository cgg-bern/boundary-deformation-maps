#!/usr/bin/bsh

if [ $# -lt 2 ]; then
    echo "Usage: run_experiment_on_original_meshes [executable_path] [original_meshes_mapped_to_a_boundary_type_path]"
    echo "For the dataset path this should be for instance path/to/original/mesh/tet_mapped/"
    exit
fi

parallel_tasks=20
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
dataset_path=$2


find $dataset_path -mindepth 1 -maxdepth 1 -type d | parallel -j$parallel_tasks --eta timeout $timeout_t \
 $bin_path -tet -v2 --skip-pre-processing \
 --domain-mesh {}/remeshed_input.ovm \
 --boundary-txt {}/boundary_conditions.txt \
 --json-stats-filepath {}/boundary_deformation_map_output.json \
 --output-domain-mesh-filepath {}/boundary_deformation_map_domain.ovm \
 --output-codomain-mesh-filepath {}/boundary_deformation_map_codomain.ovm \
 --mesh-name {/.} #> /dev/null 2>&1
