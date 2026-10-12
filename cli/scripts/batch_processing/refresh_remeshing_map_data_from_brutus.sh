
echo "---------------------------------------------"
echo "reminder: this script uses fish syntax, so run it with fish"

set base_dir /Users/val/Desktop/

set hc_dir $base_dir/hendriks_cuboids/
set hs_dir $base_dir/stars_from_hendriks_cuboids/
set kh_dir $base_dir/knotted_hole_meshes/
set tenktet_dir $base_dir/10ktetovmmeshes/
set volmap_dir $base_dir/VOLMAP/star/

echo "HC dir: " $hc_dir
echo "HS dir: " $hs_dir
echo "KH dir: " $kh_dir
echo "10ktet dir: " $tenktet_dir
echo "VOLMAP dir: " $volmap_dir

echo "cleaning up remeshing-maps data and re-downloading from brutus..."

rm $hc_dir/remeshing-maps-data/*.json
rm $hs_dir/remeshing-maps-data/*.json
rm $kh_dir/remeshing-maps-data/*.json
rm $tenktet_dir/remeshing-maps-data/*.json
rm $volmap_dir/remeshing-maps-data/*.json

echo "cleaned up, downloading json files..."

#Hendrik's cuboids
echo "---------------------------------------------"
echo "downloading Hendrik's cuboids json files..."
for file in $(ls $hc_dir/); scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/hendriks_cuboids/$(basename $file)/\*block\*boundary_deformation_map_output.json $hc_dir/remeshing-maps-data/ > /dev/null 2>&1; end
echo "downloaded $(ls $hc_dir/remeshing-maps-data/*.json | wc -l) json files"

#Hendrik's stars (download, then move and rename)
echo "---------------------------------------------"
echo "downloading Steffen's stars json files..."
for file in $(ls $hc_dir/); scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/hendriks_cuboids/$(basename $file)/star\*boundary_deformation_map_output.json $hs_dir/$(basename $file)/ > /dev/null 2>&1 ; end
echo "moving (and renaming) Steffen's stars json files..."
for file in $(ls $hs_dir/); for star in $(ls $hs_dir/$file/star*_boundary_deformation_map_output.json); mv $star $hs_dir/remeshing-maps-data/{$file}_$(basename $star) ;end   ; end
echo "downloaded $(ls $hs_dir/remeshing-maps-data/*.json | wc -l) json files"


echo "---------------------------------------------"
echo "downloading Knotted Hole json files..."
scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/knotted_hole_clean/\*boundary_deformation_map_output.json $kh_dir/remeshing-maps-data/  > /dev/null 2>&1
echo "downloaded $(ls $kh_dir/remeshing-maps-data/*.json | wc -l) json files"


# 10ktet
echo "---------------------------------------------"
echo "downloading 10ktet json files..."
scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/10ktetovmmeshes/SGA24/\*tet\*json $tenktet_dir/remeshing-maps-data/  > /dev/null 2>&1
echo "downloaded $(ls $tenktet_dir/remeshing-maps-data/*.json | wc -l) json files"


echo "---------------------------------------------"
echo "downloading VOLMAP star (G1 and G3) json files..."
scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/VOLMAP/G1_star/\*json $volmap_dir/remeshing-maps-data/  > /dev/null 2>&1
scp -r -oPort=2003 root@brutus.inf.unibe.ch:/vol/valentin/root/VOLMAP/G3_star/\*json $volmap_dir/remeshing-maps-data/  > /dev/null 2>&1
echo "downloaded $(ls $volmap_dir/remeshing-maps-data/*.json | wc -l) json files"

echo "---------------------------------------------"
