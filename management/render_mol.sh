printf "\nThis script renders a xyz file with vmd and produces \n"
printf " views from different perspectives.\n\n"

# Set the path to the external Taychon renderer!
tachyon_path=/home/jsteffen/Software/libs_local/vmd/tachyon_LINUXAMD64

# Set the resolutions if needed 
resolution_x=1500
resolution_y=1500

if [ -z "$1" ]
then
   printf "Please give a xyz file name as command line argument!\n\n"
   exit
fi 
filename=$1
rootname="${filename%.*}"
printf "Render the front perspective ...\n"
# Render the front view of the molecule
cat > tmp.vmd <<'EOF'
# load molecule
set filename [lindex $argv 0]
mol new $filename type xyz
# set display style
light 0 on
light 1 off
light 2 off
light 3 off
# position the stage and axes
axes location off
stage location off
display projection orthographic
display shadows on
display ambientocclusion on
display aoambient 0.8
display aodirect 0.7
display depthcue   off
# set representation of molecule
# colors are defined automatically by vmdrc file
mol representation CPK 1.400000 0.900000 70.000000 70.000000
mol color Element
mol material HardPlastic
mol addrep top
# render the tga file with Tachyon, transparent backround
set renderer Tachyon
# FRONT view (default orientation)
display resetview
render $renderer [file rootname ${filename}]_front
quit
EOF
vmd -dispdev rext -e tmp.vmd -args $filename
$tachyon_path -aasamples 12 -res $resolution_x $resolution_y  ${rootname}_front -format TARGA -o ${rootname}_front.tga
convert ${rootname}_front.tga -trim +repage -fuzz 3% -transparent white ${rootname}_front.png
rm ${rootname}_front.tga
rm ${rootname}_front
printf "...done! \n\n"

printf "Render the side perspective ...\n"
# Render the side view of the molecule
cat > tmp.vmd <<'EOF'
# load molecule
set filename [lindex $argv 0]
mol new $filename type xyz
# set display style
light 0 on
light 1 off
light 2 off
light 3 off
# position the stage and axes
axes location off
stage location off
display projection orthographic
display shadows on
display ambientocclusion on
display aoambient 0.8
display aodirect 0.7
display depthcue   off
# set representation of molecule
# colors are defined automatically by vmdrc file
mol representation CPK 1.400000 0.900000 70.000000 70.000000
mol color Element
mol material HardPlastic
mol addrep top
# render the tga file with Tachyon, transparent backround
set renderer Tachyon
# FRONT view (default orientation)
display resetview
rotate y by 90
render $renderer [file rootname ${filename}]_side
quit
EOF
vmd -dispdev rext -e tmp.vmd -args $filename
$tachyon_path -aasamples 12 -res $resolution_x $resolution_y  ${rootname}_side -format TARGA -o ${rootname}_side.tga
convert ${rootname}_side.tga -trim +repage -fuzz 3% -transparent white ${rootname}_side.png
rm ${rootname}_side.tga
rm ${rootname}_side
printf "...done! \n\n"

printf "Render the top perspective ...\n"
# Render the top view of the molecule
cat > tmp.vmd <<'EOF'
# load molecule
set filename [lindex $argv 0]
mol new $filename type xyz
# set display style
light 0 on
light 1 off
light 2 off
light 3 off
# position the stage and axes
axes location off
stage location off
display projection orthographic
display shadows on
display ambientocclusion on
display aoambient 0.8
display aodirect 0.7
display depthcue   off
# set representation of molecule
# colors are defined automatically by vmdrc file
mol representation CPK 1.400000 0.900000 70.000000 70.000000
mol color Element
mol material HardPlastic
mol addrep top
# render the tga file with Tachyon, transparent backround
set renderer Tachyon
# FRONT view (default orientation)
display resetview
rotate x by 90
render $renderer [file rootname ${filename}]_top
quit
EOF
vmd -dispdev rext -e tmp.vmd -args $filename
$tachyon_path -aasamples 12 -res $resolution_x $resolution_y  ${rootname}_top -format TARGA -o ${rootname}_top.tga
convert ${rootname}_top.tga -trim +repage -fuzz 3% -transparent white ${rootname}_top.png
rm ${rootname}_top.tga
rm ${rootname}_top
printf "...done! \n\n"

printf "Render the tilted perspective ...\n"
# Render the tilted view of the molecule
cat > tmp.vmd <<'EOF'
# load molecule
set filename [lindex $argv 0]
mol new $filename type xyz
# set display style
light 0 on
light 1 off
light 2 off
light 3 off
# position the stage and axes
axes location off
stage location off
display projection orthographic
display shadows on
display ambientocclusion on
display aoambient 0.8
display aodirect 0.7
display depthcue   off
# set representation of molecule
# colors are defined automatically by vmdrc file
mol representation CPK 1.400000 0.900000 70.000000 70.000000
mol color Element
mol material HardPlastic
mol addrep top
# render the tga file with Tachyon, transparent backround
set renderer Tachyon
# FRONT view (default orientation)
display resetview
rotate x by 25
rotate y by -5
render $renderer [file rootname ${filename}]_tilted
quit
EOF
vmd -dispdev rext -e tmp.vmd -args $filename
$tachyon_path -aasamples 12 -res $resolution_x $resolution_y  ${rootname}_tilted -format TARGA -o ${rootname}_tilted.tga
convert ${rootname}_tilted.tga -trim +repage -fuzz 3% -transparent white ${rootname}_tilted.png
rm ${rootname}_tilted.tga
rm ${rootname}_tilted
printf "...done! \n\n"

printf " Produced png files of different perspectives:\n"
printf " * Front view: ${rootname}_font.png \n"
printf " * Front view: ${rootname}_side.png \n"
printf " * Front view: ${rootname}_top.png \n"
printf " * Front view: ${rootname}_tilted.png \n"
printf " render_mol.sh exited normally.\n\n" 
