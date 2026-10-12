#!/usr/bin/bsh

if [ $# -lt 2 ]; then
    echo "Usage: remesh_map_HC_stars [executable_path] [Hendrik\'s_cuboids_path] (optional)[output_suffix]"
    exit
fi

parallel_tasks=100
timeout_t="0h"

echo "running $parallel_tasks in parallel, with a timeout limit of $timeout_t"

bin_path=$1
dataset_path=$2
hc_star_suffix=$3

find ${dataset_path} -name "star_codomain_*.ovm" | parallel -j$parallel_tasks --eta \
'i=$(basename {} | grep -oP "\d+(?=\.ovm)"); $bin_path -tet -v2 \
--domain-mesh $(dirname {})/star_domain_${i}.ovm \
--codomain-mesh $(dirname {})/star_codomain_${i}.ovm \
--json-stats-filepath $(dirname {})/star_${i}_boundary_deformation_map_output${hc_star_suffix}.json \
 --output-domain-mesh-filepath $(dirname {})/star_${i}_boundary_deformation_map_domain${hc_star_suffix}.ovm \
 --output-codomain-mesh-filepath $(dirname {})/star_${i}_boundary_deformation_map_codomain${hc_star_suffix}.ovm \
 --mesh-name {/} #> /dev/null 2>&1'

export hc_star_suffix=$hc_star_suffix && find ${hc_dir} -name "star_codomain_*.ovm" | parallel -j100 --eta 'i=$(basename {} | grep -oP "\d+(?=\.ovm)"); ./cli/boundary_deformation_maps -tet -v2 --domain-mesh $(dirname {})/star_domain_${i}.ovm --codomain-mesh $(dirname {})/star_codomain_${i}.ovm --json-stats-filepath $(dirname {})/star_${i}_boundary_deformation_map_output${hc_star_suffix}.json  --output-domain-mesh-filepath $(dirname {})/star_${i}_boundary_deformation_map_domain${hc_star_suffix}.ovm --output-codomain-mesh-filepath $(dirname {})/star_${i}_boundary_deformation_map_codomain${hc_star_suffix}.ovm --mesh-name {/} > /dev/null 2>&1'
