#!/usr/bin/bsh

if [ $# -lt 4 ]; then
    echo "Usage: run_experiment_on_VOLMAP_dataset [executable_path] [VOLMAP_dataset_path] [VOLMAP_group_index] [codomain_type]"
    echo "with [codomain_type] in [cube	tet	pyramid	octa star	sphere	polycube] (Note: not enforced or checked, this is just a reminder)."
    exit
fi

parallel_tasks=200
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
volmap_path=$2
volmap_idx=$3
codomain_type=$4

find $volmap_path/G$volmap_idx/*.ovm -maxdepth 0 | parallel -j200 --eta timeout 6h \
 $bin_path -tet --invert-input \
 --domain-mesh {} \
 --boundary-txt $volmap_path/G${volmap_idx}_${codomain_type}/{/.}_${codomain_type}.txt \
 --remeshed-domain-mesh-filepath $volmap_path/G${volmap_idx}/{/.}_remeshed.ovm \
 --remeshed-boundary-conditions-filepath $volmap_path/G${volmap_idx}_${codomain_type}/{/.}_${codomain_type}_remeshed.txt \
 --mesh-name {/.} #> /dev/null 2>&1
