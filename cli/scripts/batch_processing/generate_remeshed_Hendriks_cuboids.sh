#!/usr/bin/bsh

if [ $# -lt 2 ]; then
    echo "Usage: run_experiment_on_Hendriks_cuboids [executable_path] [Hendrik\'s_cuboids_path]"
    exit
fi

parallel_tasks=100
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
dataset_path=$2


find $dataset_path -mindepth 1 -maxdepth 1 -type d | parallel -j$parallel_tasks --eta timeout $timeout_t \
 $bin_path -tet -v2 --invert-input \
 --domain-mesh {}/object.ovm \
 --boundary-txt {}/boundary_conditions.txt \
 --remeshed-domain-mesh-filepath {}/remeshed_object.ovm \
 --remeshed-boundary-conditions-filepath {}/remeshed_boundary_conditions.txt \
 --mesh-name {/} #> /dev/null 2>&1
