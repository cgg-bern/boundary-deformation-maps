#!/usr/bin/bsh

if [ $# -lt 3 ]; then
    echo "Usage: remesh_map_Armadillo [executable_path] [Armadillo path]"
    exit
fi

parallel_tasks=60
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
dataset_path=$2

parallel -j$parallel_tasks --eta \
$bin_path -tet -v2 \
--domain-mesh $dataset_path/rest.ovm\
--codomain-mesh $dataset_path/{}/init.ovm\
--json-stats-filepath $dataset_path/{}/boundary_deformation_map_output.json\
--mesh-name armadillo_{} #> /dev/null 2>&1

exit

parallel -j60 --eta ./cli/boundary_deformation_maps -tet -v2 --invert-input --domain-mesh $arma_dir/rest.ovm --codomain-mesh $arma_dir/{}/init.ovm --json-stats-filepath $arma_dir/{}/boundary_deformation_map_output.json --mesh-name armadillo_{} ::: {0..600..20}

for i in {20..600..20}; do ./cli/boundary_deformation_maps -tet -v2 --invert-input --domain-mesh $arma_dir/$(($i - 20))/seq_boundary_deformation_map_output_codomain.ovm --codomain-mesh $arma_dir/$i/init.ovm --output-codomain-mesh-filepath $arma_dir/$i/seq_boundary_deformation_map_output_codomain.ovm --json-stats-filepath $arma_dir/$i/seq_boundary_deformation_map_output.json --mesh-name armadillo_seq_$i; done

