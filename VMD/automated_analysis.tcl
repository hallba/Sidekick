set location "~/Work/Sidekick/VMD/"
source ${location}TM_analysis.tcl
set numframes [expr {[molinfo top get numframes]}]
comz_distances "name \"B.*\"" "resname DPPC DPPG DPPE" bilayer_position.dat
#lipid_deformation "name \"B.*\"" "resname DPPC DPPG DPPE" up.dat down.dat
lipid_deformation_refined "name \"B.*\"" "resname DPPC DPPG DPPE" up.dat down.dat [expr $numframes/5]
helix_rotation "not resname W ION DPPC" "resname DPPC DPPG DPPE" rot.dat

rotate x by 90
scale by 2

mol delrep 0 top
mol material Diffuse
mol selection {name W}
mol representation VDW 2.4 8.0
mol addrep top
mol showperiodic top 0 xX

mol selection {name PO4}
mol addrep top
mol showperiodic top 1 xX 

mol selection {not resname DPPC DPPE DPPG ION W}
mol addrep top

mol material Ghost
mol selection {resname DPPC DPPE DPPG}
mol representation DynamicBonds 4.600000 0.300000 6.000000
mol addrep top
mol showperiodic top 3 xX

mol drawframes top 0 {now}

set numframes [expr {[molinfo top get numframes]}]
display resize 200 200
for {set i 0} {$i < $numframes} {incr i 50} {
	animate goto $i
	#Linux
	#render Tachyon $i.dat /usr/local/lib/vmd/tachyon_LINUX  -auto_skylight 0.8 %s -o $i.tga
	#Mac
	render Tachyon $i.dat '/Applications/VMD\ 1.8.6.app/Contents/vmd/tachyon_MACOSXX86'  -auto_skylight 0.8 %s -o $i.tga
	}
	
exit
