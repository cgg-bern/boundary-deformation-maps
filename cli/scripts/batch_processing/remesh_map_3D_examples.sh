#!/usr/bin/bsh

START=$(date +%s.%N)
if [ $# -lt 3 ]; then
    echo "Usage: run_remeshing_on_3D_examples [executable_path] [experiments_results_path] [delete_results_path_contents_first (y or n)]"
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

mesh_types=("minimal-non-star-shaped" "sine" "spiral" "trefoil-knot" "complex-spiral")

sine_type_args=(
"--mesh-len 10 --sine-z-scale 10.0 --sine-y-scale 1.0 --sine-period 0.5 --sine-x-scale 0.0 --torsion 0.0"
"--mesh-len 10 --sine-z-scale 10.0 --sine-y-scale 1.0 --sine-period 1.0 --sine-x-scale 0.0 --torsion 90.0"
"--mesh-len 10 --sine-z-scale 10.0 --sine-y-scale 1.0 --sine-period 1.5 --sine-x-scale 0.0 --torsion 180.0"
"--mesh-len 10 --sine-z-scale 10.0 --sine-y-scale 1.0 --sine-period 2.0 --sine-x-scale 0.0 --torsion 360.0"

"--mesh-len 20 --sine-z-scale 10.0 --sine-y-scale 4.0 --sine-period 0.5 --sine-x-scale 0.0 --torsion 0.0"
"--mesh-len 20 --sine-z-scale 10.0 --sine-y-scale 4.0 --sine-period 1.0 --sine-x-scale 0.0 --torsion 90.0"
"--mesh-len 20 --sine-z-scale 10.0 --sine-y-scale 4.0 --sine-period 1.5 --sine-x-scale 0.0 --torsion 180.0"
"--mesh-len 20 --sine-z-scale 10.0 --sine-y-scale 4.0 --sine-period 2.0 --sine-x-scale 0.0 --torsion 360.0"

"--mesh-len 20 --sine-z-scale 20.0 --sine-y-scale 2.0 --sine-period 0.5 --sine-x-scale 0.0 --torsion 0.0"
"--mesh-len 20 --sine-z-scale 20.0 --sine-y-scale 2.0 --sine-period 1.0 --sine-x-scale 0.0 --torsion 90.0"
"--mesh-len 20 --sine-z-scale 20.0 --sine-y-scale 2.0 --sine-period 1.5 --sine-x-scale 0.0 --torsion 180.0"
"--mesh-len 20 --sine-z-scale 20.0 --sine-y-scale 2.0 --sine-period 2.0 --sine-x-scale 0.0 --torsion 360.0"

"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 2.5 --sine-x-scale 0.0 --torsion 90.0"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 3.0 --sine-x-scale 0.0 --torsion 180.0"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 3.5 --sine-x-scale 0.0 --torsion 270.0"
"--mesh-len 30 --sine-z-scale 30.0 --sine-y-scale 2.0 --sine-period 4.0 --sine-x-scale 0.0 --torsion 360.0"
)

spiral_type_args=(
"--mesh-len 5 --spiral-turns 0.5 --spiral-x-scale 1.0 --torsion 0.0"
"--mesh-len 10 --spiral-turns 0.5 --spiral-x-scale 2.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 0.5 --spiral-x-scale 5.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 0.5 --spiral-x-scale 10.0 --torsion 360.0"

"--mesh-len 5 --spiral-turns 1.0 --spiral-x-scale 1.0 --torsion 0.0"
"--mesh-len 10 --spiral-turns 1.0 --spiral-x-scale 2.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 1.0 --spiral-x-scale 5.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 1.0 --spiral-x-scale 10.0 --torsion 360.0"

"--mesh-len 10 --spiral-turns 1.5 --spiral-x-scale 1.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 1.5 --spiral-x-scale 2.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 1.5 --spiral-x-scale 5.0 --torsion 360.0"
"--mesh-len 25 --spiral-turns 1.5 --spiral-x-scale 10.0 --torsion 450.0"

"--mesh-len 10 --spiral-turns 2.0 --spiral-x-scale 1.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 2.0 --spiral-x-scale 2.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 2.0 --spiral-x-scale 5.0 --torsion 360.0"
"--mesh-len 25 --spiral-turns 2.0 --spiral-x-scale 10.0 --torsion 450.0"

"--mesh-len 10 --spiral-turns 2.5 --spiral-x-scale 1.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 2.5 --spiral-x-scale 2.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 2.5 --spiral-x-scale 5.0 --torsion 360.0"
"--mesh-len 25 --spiral-turns 2.5 --spiral-x-scale 10.0 --torsion 450.0"

"--mesh-len 10 --spiral-turns 3.0 --spiral-x-scale 1.0 --torsion 90.0"
"--mesh-len 15 --spiral-turns 3.0 --spiral-x-scale 2.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 3.0 --spiral-x-scale 5.0 --torsion 360.0"
"--mesh-len 25 --spiral-turns 3.0 --spiral-x-scale 10.0 --torsion 450.0"

"--mesh-len 15 --spiral-turns 3.5 --spiral-x-scale 1.0 --torsion 180.0"
"--mesh-len 20 --spiral-turns 3.5 --spiral-x-scale 2.0 --torsion 360.0"
"--mesh-len 25 --spiral-turns 3.5 --spiral-x-scale 5.0 --torsion 450.0"
"--mesh-len 30 --spiral-turns 3.5 --spiral-x-scale 10.0 --torsion 540.0"
)



trefoil_knot_type_args=(
"--mesh-len 20 --parametric-range 0.2 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 20 --parametric-range 0.5 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 20 --parametric-range 0.9 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 20 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"

"--mesh-len 20 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 30 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 40 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"

"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.1 --torsion 0.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.2 --torsion 0.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 0.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.7 --torsion 0.0"

"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 90.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 180.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 360.0"
"--mesh-len 50 --parametric-range 0.99 --parametric-width 0.5 --torsion 720.0"
)


complex_spiral_type_args=(

"--mesh-len 10 --parametric-range 0.01 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 20 --parametric-range 0.05 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 80 --parametric-range 0.2 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 160 --parametric-range 0.4 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"

"--mesh-len 30 --parametric-range 0.05 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 90.0"
"--mesh-len 30 --parametric-range 0.05 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 180.0"
"--mesh-len 30 --parametric-range 0.05 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 360.0"
"--mesh-len 30 --parametric-range 0.05 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 720.0"

"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 60 --parametric-range 0.1 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 80 --parametric-range 0.1 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 100 --parametric-range 0.1 --parametric-width 0.3 --parametric-z-scale 5.0 --torsion 0.0"

"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.1 --parametric-z-scale 1.0 --torsion 0.0"
"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.2 --parametric-z-scale 3.0 --torsion 0.0"
"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.4 --parametric-z-scale 5.0 --torsion 0.0"
"--mesh-len 40 --parametric-range 0.1 --parametric-width 0.5 --parametric-z-scale 7.0 --torsion 0.0"
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
parallel -j8 --eta eval '{1} -tet --test-mesh "sine" {3} --json-stats-filepath "{2}/sine_{#}_output.json" --output-domain-mesh-filepath "{2}/sine_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/sine_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${sine_type_args[@]}"
SINE_END=$(date +%s.%N)
DIFF=$(echo "$SINE_END - $SINE_START" | bc)
echo " -> sine mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "



echo " ------------------------------------------------------ "
SPIRAL_START=$(date +%s.%N)
echo " running method on spiral meshes..."
parallel -j8 --eta eval '{1} -tet --test-mesh "spiral" {3} --json-stats-filepath "{2}/spiral_{#}_output.json" --output-domain-mesh-filepath "{2}/spiral_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/spiral_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${spiral_type_args[@]}"
SPIRAL_END=$(date +%s.%N)
DIFF=$(echo "$SPIRAL_END - $SPIRAL_START" | bc)
echo " -> spiral mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "


echo " ------------------------------------------------------ "
TREFOIL_START=$(date +%s.%N)
echo " running method on trefoil knot meshes..."
parallel -j8 --eta eval '{1} -tet --test-mesh "trefoil-knot" {3} --json-stats-filepath "{2}/trefoil_knot_{#}_output.json" --output-domain-mesh-filepath "{2}/trefoil_knot_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/trefoil_knot_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${trefoil_knot_type_args[@]}"
TREFOIL_END=$(date +%s.%N)
DIFF=$(echo "$TREFOIL_END - $TREFOIL_START" | bc)
echo " -> trefoil knot mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "


echo " ------------------------------------------------------ "
COMPLEX_SPIRAL_START=$(date +%s.%N)
echo " running method on complex spiral meshes..."
parallel -j8 --eta eval '{1} -tet --test-mesh "spiral" {3} --json-stats-filepath "{2}/complex_spiral_{#}_output.json" --output-domain-mesh-filepath "{2}/complex_spiral_{#}_domain.ovm" --output-codomain-mesh-filepath "{2}/complex_spiral_{#}_codomain.ovm" > /dev/null 2>&1' ::: $bin_path ::: $results_path ::: "${complex_spiral_type_args[@]}"
COMPLEX_SPIRAL_END=$(date +%s.%N)
DIFF=$(echo "$COMPLEX_SPIRAL_END - $COMPLEX_SPIRAL_START" | bc)
echo " -> complex spiral mesh tests ran in $DIFF seconds"
echo " ------------------------------------------------------ "

echo " ------------------------------------------------------ "
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "full experiment ran in $DIFF seconds"
echo "================================================================================================================"
echo "================================================================================================================"
