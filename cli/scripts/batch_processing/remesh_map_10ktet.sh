#!/usr/bin/bsh

if [ $# -lt 3 ]; then
    echo "Usage: remesh_map_10ktet [executable_path] [10ktet_dataset_path] [codomain_type]"
    echo "with [10ktet_dataset_path] pointing to a directory with all .ovm meshes and .txt boundary conditions files"
    echo "and [codomain_type] in [tet	stiff_tet	ball random_star_shape] (Note: not enforced or checked, this is just a reminder)."
    exit
fi

parallel_tasks=120
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
data_path=$2
codomain_type=$3

#find $data_path/*tetmesh.ovm -maxdepth 0 | parallel echo {}

ls -S -r $data_path/*tetmesh.ovm  | parallel -j${parallel_tasks} --eta \
 $bin_path -v2 -tet \
 --domain-mesh {} \
 --boundary-txt                  {.}_${codomain_type}_boundary.txt \
 --json-stats-filepath           {.}_${codomain_type}_boundary_boundary_deformation_map_output.json \
 --output-domain-mesh-filepath   {.}_${codomain_type}_boundary_boundary_deformation_map_domain.ovm \
 --output-codomain-mesh-filepath {.}_${codomain_type}_boundary_deformation_map_codomain.ovm \
 --mesh-name {/.} #> /dev/null 2>&1
