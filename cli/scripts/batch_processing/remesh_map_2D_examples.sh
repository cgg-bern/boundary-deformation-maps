#!/usr/bin/bsh

START=$(date +%s.%N)
if [ $# -lt 3 ]; then
    echo "Usage: run_remeshing_on_2D_examples [executable_path] [experiments_results_path] [delete_results_path_contents_first (y or n)]"
    exit
fi

bin_path=$1
results_path=$2
erase_res_path=$3

if [ ! -d "$results_path" ]; then
  echo "results path $results_path does not exist."
  exit
fi


if [ $erase_res_path = "y" ]; then
  rm -rf $results_path/*
fi

mesh_types=("minimal-non-star-shaped" "sine" "spiral" "clover" "paraspiral" "compression")

sine_type_args=(
"--mesh-len 10 --sine-x-scale 10.0 --sine-y-scale 1.0 --sine-period 0.5"
"--mesh-len 10 --sine-x-scale 10.0 --sine-y-scale 1.0 --sine-period 1.0"
"--mesh-len 10 --sine-x-scale 10.0 --sine-y-scale 1.0 --sine-period 1.5"
"--mesh-len 10 --sine-x-scale 10.0 --sine-y-scale 1.0 --sine-period 2.0"

"--mesh-len 20 --sine-x-scale 10.0 --sine-y-scale 4.0 --sine-period 0.5"
"--mesh-len 20 --sine-x-scale 10.0 --sine-y-scale 4.0 --sine-period 1.0"
"--mesh-len 20 --sine-x-scale 10.0 --sine-y-scale 4.0 --sine-period 1.5"
"--mesh-len 20 --sine-x-scale 10.0 --sine-y-scale 4.0 --sine-period 2.0"

"--mesh-len 20 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 0.5"
"--mesh-len 20 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 1.0"
"--mesh-len 20 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 1.5"
"--mesh-len 20 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 2.0"

"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 2.5"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 3.0"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 3.5"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 4.0"
)

spiral_type_args=(
"--mesh-len 5 --spiral-turns 0.5"
"--mesh-len 10 --spiral-turns 0.5"
"--mesh-len 15 --spiral-turns 0.5"
"--mesh-len 20 --spiral-turns 0.5"

"--mesh-len 5 --spiral-turns 0.75"
"--mesh-len 10 --spiral-turns 0.75"
"--mesh-len 15 --spiral-turns 0.75"
"--mesh-len 20 --spiral-turns 0.75"

"--mesh-len 5 --spiral-turns 1.0"
"--mesh-len 10 --spiral-turns 1.0"
"--mesh-len 15 --spiral-turns 1.0"
"--mesh-len 20 --spiral-turns 1.0"

"--mesh-len 10 --spiral-turns 1.5"
"--mesh-len 15 --spiral-turns 1.5"
"--mesh-len 20 --spiral-turns 1.5"
"--mesh-len 25 --spiral-turns 1.5"

"--mesh-len 10 --spiral-turns 2.0"
"--mesh-len 15 --spiral-turns 2.0"
"--mesh-len 20 --spiral-turns 2.0"
"--mesh-len 25 --spiral-turns 2.0"

"--mesh-len 10 --spiral-turns 2.5"
"--mesh-len 15 --spiral-turns 2.5"
"--mesh-len 20 --spiral-turns 2.5"
"--mesh-len 25 --spiral-turns 2.5"

"--mesh-len 10 --spiral-turns 3.0"
"--mesh-len 15 --spiral-turns 3.0"
"--mesh-len 20 --spiral-turns 3.0"
"--mesh-len 25 --spiral-turns 3.0"

"--mesh-len 15 --spiral-turns 3.5"
"--mesh-len 20 --spiral-turns 3.5"
"--mesh-len 25 --spiral-turns 3.5"
"--mesh-len 30 --spiral-turns 3.5"
)


clover_type_args=(
"--mesh-len 20 --clover-center-radius 0.1 --clover-leaves 4 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 0.2 --clover-leaves 4 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 0.5 --clover-leaves 4 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 1.0 --clover-leaves 4 --clover-shift 4.0"

"--mesh-len 20 --clover-center-radius 0.1 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 0.2 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 0.5 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 20 --clover-center-radius 1.0 --clover-leaves 8 --clover-shift 4.0"

"--mesh-len 40 --clover-center-radius 0.1 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 40 --clover-center-radius 0.2 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 40 --clover-center-radius 0.5 --clover-leaves 8 --clover-shift 4.0"
"--mesh-len 40 --clover-center-radius 1.0 --clover-leaves 8 --clover-shift 4.0"
)


mouth_type_args=(
"--mesh-len 40 --mouth-teeth-count 5 --mouth-opening 0.1"
"--mesh-len 40 --mouth-teeth-count 10 --mouth-opening 0.1"
"--mesh-len 40 --mouth-teeth-count 15 --mouth-opening 0.1"

"--mesh-len 60 --mouth-teeth-count 5 --mouth-opening 0.1"
"--mesh-len 60 --mouth-teeth-count 10 --mouth-opening 0.1"
"--mesh-len 60 --mouth-teeth-count 15 --mouth-opening 0.1"

"--mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.1"
"--mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.5"
"--mesh-len 100 --mouth-teeth-count 15 --mouth-opening 1.0"
)


paraspiral_type_args=(
"--mesh-len 30 --paraspiral-center-radius 1.0 --paraspiral-branches 3 --paraspiral-branch-len 0.5"
"--mesh-len 50 --paraspiral-center-radius 1.0 --paraspiral-branches 5 --paraspiral-branch-len 0.5"
"--mesh-len 100 --paraspiral-center-radius 1.0 --paraspiral-branches 10 --paraspiral-branch-len 0.5"

"--mesh-len 30 --paraspiral-center-radius 1.0 --paraspiral-branches 3 --paraspiral-branch-len 1.0"
"--mesh-len 50 --paraspiral-center-radius 1.0 --paraspiral-branches 5 --paraspiral-branch-len 1.0"
"--mesh-len 100 --paraspiral-center-radius 1.0 --paraspiral-branches 10 --paraspiral-branch-len 1.0"

"--mesh-len 60 --paraspiral-center-radius 1.0 --paraspiral-branches 3 --paraspiral-branch-len 2.0"
"--mesh-len 100 --paraspiral-center-radius 1.0 --paraspiral-branches 5 --paraspiral-branch-len 2.0"
"--mesh-len 200 --paraspiral-center-radius 1.0 --paraspiral-branches 10 --paraspiral-branch-len 2.0"

"--mesh-len 30 --paraspiral-center-radius 2.0 --paraspiral-branches 3 --paraspiral-branch-len 0.5"
"--mesh-len 50 --paraspiral-center-radius 2.0 --paraspiral-branches 5 --paraspiral-branch-len 0.5"
"--mesh-len 100 --paraspiral-center-radius 2.0 --paraspiral-branches 10 --paraspiral-branch-len 0.5"
)


compression_type_args=(
"--test-mesh \"sine\" --mesh-len 30 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 2.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"sine\" --mesh-len 30 --sine-x-scale 20.0 --sine-y-scale 2.0 --sine-period 2.0 --compression-mode 1 --i-to-b-v-count-ratio 100"
"--test-mesh \"sine\" --mesh-len 30 --sine-x-scale 20.0 --sine-y-scale 5.0 --sine-period 1.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"sine\" --mesh-len 30 --sine-x-scale 20.0 --sine-y-scale 5.0 --sine-period 1.0 --compression-mode 1 --i-to-b-v-count-ratio 100"

"--test-mesh \"spiral\" --mesh-len 20 --spiral-turns 1.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"spiral\" --mesh-len 20 --spiral-turns 1.0 --compression-mode 1 --i-to-b-v-count-ratio 100"
"--test-mesh \"spiral\" --mesh-len 40 --spiral-turns 4.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"spiral\" --mesh-len 40 --spiral-turns 4.0 --compression-mode 1 --i-to-b-v-count-ratio 100"

"--test-mesh \"clover\" --mesh-len 60 --clover-center-radius 0.1 --clover-leaves 8 --clover-shift 4.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"clover\" --mesh-len 60 --clover-center-radius 0.1 --clover-leaves 8 --clover-shift 4.0 --compression-mode 1 --i-to-b-v-count-ratio 100"
"--test-mesh \"clover\" --mesh-len 60 --clover-center-radius 1.0 --clover-leaves 4 --clover-shift 4.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"clover\" --mesh-len 60 --clover-center-radius 1.0 --clover-leaves 4 --clover-shift 4.0 --compression-mode 1 --i-to-b-v-count-ratio 100"

"--test-mesh \"mouth\" --mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.5 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"mouth\" --mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.5 --compression-mode 1 --i-to-b-v-count-ratio 100"
"--test-mesh \"mouth\" --mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.1 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"mouth\" --mesh-len 100 --mouth-teeth-count 15 --mouth-opening 0.1 --compression-mode 1 --i-to-b-v-count-ratio 100"

"--test-mesh \"paraspiral\" --mesh-len 30 --paraspiral-center-radius 1.0 --paraspiral-branches 3 --paraspiral-branch-len 0.5 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"paraspiral\" --mesh-len 30 --paraspiral-center-radius 1.0 --paraspiral-branches 3 --paraspiral-branch-len 0.5 --compression-mode 1 --i-to-b-v-count-ratio 50"
"--test-mesh \"paraspiral\" --mesh-len 50 --paraspiral-center-radius 0.1 --paraspiral-branches 5 --paraspiral-branch-len 1.0 --compression-mode 1 --i-to-b-v-count-ratio 20"
"--test-mesh \"paraspiral\" --mesh-len 50 --paraspiral-center-radius 0.1 --paraspiral-branches 5 --paraspiral-branch-len 1.0 --compression-mode 1 --i-to-b-v-count-ratio 100"
)

tab="     "

echo "================================================================================================================"
echo "================================================================================================================"

echo " binary path: $bin_path"
echo " experiments results path: $results_path"
echo "---------------------------"
echo " testing method on test mesh types: ${mesh_types[@]}"
echo "---------------------------"
#echo "  * sine mesh arguments: "
#for ((i = 0; i < ${#sine_type_args[@]}; i++))
#do
#    echo "$tab ${sine_type_args[$i]}"
#done
#echo "---------------------------"
#echo "  * spiral mesh arguments: "
#for ((i = 0; i < ${#spiral_type_args[@]}; i++))
#do
#    echo "$tab ${spiral_type_args[$i]}"
#done
#echo "---------------------------"
echo "---------------------------"

echo " ------------------------------------------------------ "
echo " running method on minimal non-star-shaped mesh..."
eval $bin_path --test-mesh ${mesh_types[0]} --json-stats-filepath min_non_starshaped_output.json > /dev/null 2>&1
NSS_END=$(date +%s.%N)
DIFF=$(echo "$NSS_END - $START" | bc)
echo "minimal non-star-shaped test ran in $DIFF seconds"



echo " ------------------------------------------------------ "
SINE_START=$(date +%s.%N)
echo " running method on sine meshes..."
parallel -j8 --eta eval '{1} --test-mesh "sine" {3} --json-stats-filepath "{2}/sine_{#}_output.json" --output-domain-mesh-filepath "{2}/sine_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/sine_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${sine_type_args[@]}"
SINE_END=$(date +%s.%N)
DIFF=$(echo "$SINE_END - $SINE_START" | bc)
echo " -> sine mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "



echo " ------------------------------------------------------ "
SPIRAL_START=$(date +%s.%N)
echo " running method on spiral meshes..."
parallel -j8 --eta eval '{1} --test-mesh "spiral" {3} --json-stats-filepath "{2}/spiral_{#}_output.json" --output-domain-mesh-filepath "{2}/spiral_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/spiral_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${spiral_type_args[@]}"
SPIRAL_END=$(date +%s.%N)
DIFF=$(echo "$SPIRAL_END - $SPIRAL_START" | bc)
echo " -> spiral mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "


echo " ------------------------------------------------------ "
CLOVER_START=$(date +%s.%N)
echo " running method on clover meshes..."
parallel -j8 --eta eval '{1} --test-mesh "clover" {3} --json-stats-filepath "{2}/clover_{#}_output.json" --output-domain-mesh-filepath "{2}/clover_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/clover_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${clover_type_args[@]}"
CLOVER_END=$(date +%s.%N)
DIFF=$(echo "$CLOVER_END - $CLOVER_START" | bc)
echo " -> clover mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "


echo " ------------------------------------------------------ "
MOUTH_START=$(date +%s.%N)
echo " running method on mouth meshes..."
parallel -j8 --eta eval '{1} --test-mesh "mouth" {3} --json-stats-filepath "{2}/mouth_{#}_output.json" --output-domain-mesh-filepath "{2}/mouth_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/mouth_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${mouth_type_args[@]}"
MOUTH_END=$(date +%s.%N)
DIFF=$(echo "$MOUTH_END - $MOUTH_START" | bc)
echo " -> mouth mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "


echo " ------------------------------------------------------ "
PARASPIRAL_START=$(date +%s.%N)
echo " running method on parametric spiral meshes..."
parallel -j8 --eta eval '{1} --test-mesh "paraspiral" {3} --json-stats-filepath "{2}/paraspiral_{#}_output.json" --output-domain-mesh-filepath "{2}/paraspiral_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/paraspiral_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${paraspiral_type_args[@]}"
PARASPIRAL_END=$(date +%s.%N)
DIFF=$(echo "$PARASPIRAL_END - $PARASPIRAL_START" | bc)
echo " -> parametric spiral mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "

echo " ------------------------------------------------------ "
COMPRESSION_START=$(date +%s.%N)
echo " running method on compression meshes..."
parallel -j8 --eta eval '{1} {3} --json-stats-filepath "{2}/paraspiral_{#}_output.json" --output-domain-mesh-filepath "{2}/paraspiral_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/paraspiral_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${compression_type_args[@]}"
COMPRESSION_END=$(date +%s.%N)
DIFF=$(echo "$COMPRESSION_END - $COMPRESSION_START" | bc)
echo " -> compression mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "

echo " ------------------------------------------------------ "
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "full experiment ran in $DIFF seconds"
echo "================================================================================================================"
echo "================================================================================================================"
