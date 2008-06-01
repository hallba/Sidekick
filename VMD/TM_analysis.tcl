proc comz_distances {proteinselection lipidselection outfile_name} {
	
	proc main {proteinselection lipidselection outfile_name} {
		set numframes [expr {[molinfo top get numframes]}]
		puts "total number of frames: $numframes"
		set outfile [open $outfile_name w]
		
		set protein [atomselect top $proteinselection ]
		set lipid [atomselect top $lipidselection ] 
		set all [atomselect top "all" ]
		
		set increment 1
		
		for {set i 0} {$i < $numframes} {incr i $increment} {
			$all frame $i
			$all update
			$protein frame $i
			$protein update
			$lipid frame $i
			$lipid update
			
			
			set pcom [measure center $protein]
			set lcom [measure center $lipid]
			set topcom [measure center [atomselect top "resid 1" frame $i]]
			
			set pz [lindex $pcom 2]
			set lz [lindex $lcom 2]
			set topz [lindex $topcom 2]
			
			if { $pz > $lz } {
			set com_distance [expr $pz - $lz]
			} else {
			set com_distance [expr $lz - $pz]
			}
			
			if { $pz < $topz && $pz > $lz } {
				set com_distance [expr $com_distance]
			} elseif { $pz > $topz && $pz < $lz } {
				set com_distance [expr $com_distance]
			} else {
				set com_distance [expr 0 - $com_distance]
			}

			
			puts $outfile "$i $com_distance"
			flush $outfile
			if { $i == 0 } {continue}
			if { [expr $i % [expr 10 * $increment]] == 0 } { puts -nonewline "." }
      			if { [expr $i % [expr 500 * $increment]] == 0 } { puts " " }
			flush stdout
    			}
		}
		
	main $proteinselection $lipidselection $outfile_name
}

#comz_distances "name \"B.*\"" "resname DPPC" bilayer_position.dat

proc lipid_deformation {proteinselection lipidselection upper_outfile_name lower_outfile_name} {
	
	proc fill_boxes {upper_radial_boxes lower_radial_boxes lipid_list px py lz frame topz} {

			foreach resid $lipid_list {
				set resx [[atomselect top "resid $resid and name PO4" frame $frame] get x]
				set resy [[atomselect top "resid $resid and name PO4" frame $frame] get y]
				set resz [[atomselect top "resid $resid and name PO4" frame $frame] get z]
				
				#which bin/box? calculate distance, make positive, then floor
				set d2 [expr ($resx - $px)*($resx - $px) + ($resy - $py)*($resy - $py)]
				#puts "Distance2 resx px resy py $d2 $resx $px $resy $py"
				#if {$d2 < 0} { set d2 [expr 0 - $d2] }
				
				set distance [expr sqrt($d2) ]
				set target_box [expr $distance / 5]
				set target_box [expr int($target_box)]
				
				set bin_no [llength $upper_radial_boxes]
				set thickness [expr $resz - $lz]
				if {$thickness < 0} { set thickness [expr 0 - $thickness] }
				#puts "T $target_box L $bin_no"
				if {$target_box >= $bin_no} { continue }
				if {$resz > $lz && $topz > $lz} {
					#lappend upper_leaflet $resid
					
					#lappend [lindex $upper_radial_boxes $target_box] [$resz - $lz]
					#set upper_radial_boxes [lreplace $upper_radial_boxes $target_box $target_box [concat [lindex $upper_radial_boxes $target_box] [list $thickness ] ] ]
					lset upper_radial_boxes $target_box [concat [lindex $upper_radial_boxes $target_box] [list $thickness ] ] 
				} elseif {$resz < $lz && $topz < $lz} {
					#lappend upper_leaflet $resid

					#lappend [lindex $upper_radial_boxes $target_box] [$resz - $lz]
					#set upper_radial_boxes [lreplace $upper_radial_boxes $target_box $target_box [concat [lindex $upper_radial_boxes $target_box] [list $thickness ] ] ]
					lset upper_radial_boxes $target_box [concat [lindex $upper_radial_boxes $target_box] [list $thickness ] ] 
				} else {
					#lappend lower_leaflet $resid
					#lappend [lindex $lower_radial_boxes $target_box] [$lz - $resz]
					#set lower_radial_boxes [lreplace $lower_radial_boxes $target_box $target_box [concat [lindex $lower_radial_boxes $target_box] [list $thickness ] ] ]
					lset lower_radial_boxes $target_box [concat [lindex $lower_radial_boxes $target_box] [list $thickness ] ] 
				}
			}
		#puts "lz resz $lz $resz"
		return [list $upper_radial_boxes $lower_radial_boxes]
	}
	
	proc main {proteinselection lipidselection upper_outfile_name lower_outfile_name} {
		set numframes [expr {[molinfo top get numframes]}]
		puts "total number of frames: $numframes"
		set upper_outfile [open $upper_outfile_name w]
		set lower_outfile [open $lower_outfile_name w]
		
		set protein [atomselect top $proteinselection ]
		set lipid [atomselect top $lipidselection ] 
		set all [atomselect top "all" ]
		
		set increment 1
		
		if { [molinfo top get a] > [molinfo top get b] } {
			set box_rad [molinfo top get b]
		} else {
			set box_rad [molinfo top get a]
		}
		

		
		set lipid_list [lsort -integer -unique [$lipid get resid] ]
		
		for {set i 0} {$i < $numframes} {incr i $increment} {
			$all frame $i
			$all update
			$protein frame $i
			$protein update
			$lipid frame $i
			$lipid update
			
			
			set pcom [measure center $protein]
			set lcom [measure center $lipid]
			set topcom [measure center [atomselect top "resid 1" frame $i]]
			
			set px [lindex $pcom 0]
			set py [lindex $pcom 1]
			set pz [lindex $pcom 2]
			set lz [lindex $lcom 2]
			set topz [lindex $topcom 2]
			
			if { $pz > $lz } {
			set com_distance [expr $pz - $lz]
			} else {
			set com_distance [expr $lz - $pz]
			}
			
			#set upper_leaflet [list]
			#set lower_leaflet [list]
			
			set upper_radial_boxes [list]
			set lower_radial_boxes [list]
			#5A banding for concentric circles
			for {set j 0 } {$j < ($box_rad/2) } { incr j 5} {
				lappend upper_radial_boxes [list]
				lappend lower_radial_boxes [list]
			}
		
			#set number_of_bins [llength $radial_boxes]
			set full_boxes [fill_boxes $upper_radial_boxes $lower_radial_boxes $lipid_list $px $py $lz $i $topz]
			set upper_radial_boxes [lindex $full_boxes 0]
			set lower_radial_boxes [lindex $full_boxes 1]
			
			#commented out and put in a proc due to a massive memory leak
			if 0 {
			#The following code attempts to append items to lists in the list.
			#This is at least one route to the dark side
			foreach resid $lipid_list {
				set resx [[atomselect top "resid $resid and name PO4"] get x]
				set resy [[atomselect top "resid $resid and name PO4"] get y]
				set resz [[atomselect top "resid $resid and name PO4"] get z]
				#which bin/box?
				set d2 [expr ($resx - $px)*($resx - $px) + ($resy - $py)*($resy - $py)]
				if {$d2 < 0} { set d2 [expr 0 - $d2] }
				set distance [expr sqrt($d2) ]
				set target_box [expr $distance / 5]
				set target_box [expr int($target_box)]
				set bin_no [llength $upper_radial_boxes]
				#puts "T $target_box L $bin_no"
				if {$target_box >= $bin_no} { continue }
				if {$resz > $lz} {
					#lappend upper_leaflet $resid

					#lappend [lindex $upper_radial_boxes $target_box] [$resz - $lz]
					set upper_radial_boxes [lreplace $upper_radial_boxes $target_box $target_box [concat [lindex $upper_radial_boxes $target_box] [list [expr $resz - $lz] ] ] ]
				} else {
					#lappend lower_leaflet $resid
					#lappend [lindex $lower_radial_boxes $target_box] [$lz - $resz]
					set lower_radial_boxes [lreplace $lower_radial_boxes $target_box $target_box [concat [lindex $lower_radial_boxes $target_box] [list [expr $lz - $resz] ] ] ]
				}
			}
			#Let the hate flow through you
			}
			
			#average the bins and append to lists for printing
			
			set upper_leaflet [list]
			set lower_leaflet [list]

			foreach bin $upper_radial_boxes {
				set len [llength $bin]
				if {$len == 0} {lappend upper_leaflet 0 ; continue}
				set sum [expr [join $bin +]]
				set average [expr $sum / [llength $bin] ]
				lappend upper_leaflet $average
			}
			
			foreach bin $lower_radial_boxes {
				set len [llength $bin]
				if {$len == 0} {lappend lower_leaflet 0 ; continue}
				set sum [expr [join $bin +]]
				set average [expr $sum / [llength $bin] ]
				lappend lower_leaflet $average
			}
			puts $upper_outfile "$upper_leaflet"
			puts $lower_outfile "$lower_leaflet"
			flush $upper_outfile
			flush $lower_outfile
						
			if { $i == 0 } {continue}
			if { [expr $i % [expr 10 * $increment]] == 0 } { puts -nonewline "." }
      			if { [expr $i % [expr 500 * $increment]] == 0 } { puts " " }
			
			flush stdout
    			}
		}
		
	main $proteinselection $lipidselection $upper_outfile_name $lower_outfile_name
}

proc lipid_deformation_test {proteinselection lipidselection upper_outfile_name lower_outfile_name} {
	
	proc fill_boxes {upper_radial_boxes lower_radial_boxes lipid_list px py lz} {

			foreach resid $lipid_list {
				set resx [[atomselect top "resid $resid and name PO4"] get x]
				set resy [[atomselect top "resid $resid and name PO4"] get y]
				set resz [[atomselect top "resid $resid and name PO4"] get z]
				
				#which bin/box? calculate distance, make positive, then floor
				set d2 [expr ($resx - $px)*($resx - $px) + ($resy - $py)*($resy - $py)]
				#No reason to force positivity on it...
				#puts "Distance2 resx px resy py $d2 $resx $px $resy $py"
				#if {$d2 < 0} { set d2 [expr 0 - $d2] }
				set distance [expr sqrt($d2) ]
				set target_box [expr $distance / 5]
				set target_box [expr int($target_box)]
				
				set bin_no [llength $upper_radial_boxes]
				#puts "T $target_box L $bin_no"
				if {$target_box >= $bin_no} { continue }
				
				puts "lz resz $lz $resz"
				if {$resz > $lz} {
					#lappend upper_leaflet $resid

					#lappend [lindex $upper_radial_boxes $target_box] [$resz - $lz]
					set upper_radial_boxes [lreplace $upper_radial_boxes $target_box $target_box [concat [lindex $upper_radial_boxes $target_box] [list [expr $resz - $lz] ] ] ]
				} else {
					#lappend lower_leaflet $resid
					#lappend [lindex $lower_radial_boxes $target_box] [$lz - $resz]
					set lower_radial_boxes [lreplace $lower_radial_boxes $target_box $target_box [concat [lindex $lower_radial_boxes $target_box] [list [expr $lz - $resz] ] ] ]
				}
				
			}
		return [list $upper_radial_boxes $lower_radial_boxes]
	}
	
	proc main {proteinselection lipidselection upper_outfile_name lower_outfile_name} {
		set numframes [expr {[molinfo top get numframes]}]
		puts "total number of frames: $numframes"
		set upper_outfile [open $upper_outfile_name w]
		set lower_outfile [open $lower_outfile_name w]
		
		set protein [atomselect top $proteinselection ]
		set lipid [atomselect top $lipidselection ] 
		set all [atomselect top "all" ]
		
		set increment 1
		
		if { [molinfo top get a] > [molinfo top get b] } {
			set box_rad [molinfo top get b]
		} else {
			set box_rad [molinfo top get a]
		}
		

		
		set lipid_list [lsort -integer -unique [$lipid get resid] ]
		
		for {set i 0} {$i < 1} {incr i $increment} {
			$all frame $i
			$all update
			$protein frame $i
			$protein update
			$lipid frame $i
			$lipid update
			
			
			set pcom [measure center $protein]
			set lcom [measure center $lipid]
			set topcom [measure center [atomselect top "resid 1"]]
			
			set px [lindex $pcom 0]
			set py [lindex $pcom 1]
			set pz [lindex $pcom 2]
			set lz [lindex $lcom 2]
			set topz [lindex $topcom 2]
			
			if { $pz > $lz } {
			set com_distance [expr $pz - $lz]
			} else {
			set com_distance [expr $lz - $pz]
			}
			
			#set upper_leaflet [list]
			#set lower_leaflet [list]
			
			set upper_radial_boxes [list]
			set lower_radial_boxes [list]
			#5A banding for concentric circles
			for {set j 0 } {$j < ($box_rad/2) } { incr j 5} {
				lappend upper_radial_boxes [list]
				lappend lower_radial_boxes [list]
			}
		
			#set number_of_bins [llength $radial_boxes]
			set full_boxes [fill_boxes $upper_radial_boxes $lower_radial_boxes $lipid_list $px $py $lz]
			set upper_radial_boxes [lindex $full_boxes 0]
			set lower_radial_boxes [lindex $full_boxes 1]

			if 0 {
			#The following code attempts to append items to lists in the list.
			#This is at least one route to the dark side
			foreach resid $lipid_list {
				set resx [[atomselect top "resid $resid and name PO4"] get x]
				set resy [[atomselect top "resid $resid and name PO4"] get y]
				set resz [[atomselect top "resid $resid and name PO4"] get z]
				#which bin/box?
				set d2 [expr ($resx - $px)*($resx - $px) + ($resy - $py)*($resy - $py)]
				if {$d2 < 0} { set d2 [expr 0 - $d2] }
				set distance [expr sqrt($d2) ]
				set target_box [expr $distance / 5]
				set target_box [expr int($target_box)]
				set bin_no [llength $upper_radial_boxes]
				#puts "T $target_box L $bin_no"
				if {$target_box >= $bin_no} { continue }
				if {$resz > $lz} {
					#lappend upper_leaflet $resid

					#lappend [lindex $upper_radial_boxes $target_box] [$resz - $lz]
					set upper_radial_boxes [lreplace $upper_radial_boxes $target_box $target_box [concat [lindex $upper_radial_boxes $target_box] [list [expr $resz - $lz] ] ] ]
				} else {
					#lappend lower_leaflet $resid
					#lappend [lindex $lower_radial_boxes $target_box] [$lz - $resz]
					set lower_radial_boxes [lreplace $lower_radial_boxes $target_box $target_box [concat [lindex $lower_radial_boxes $target_box] [list [expr $lz - $resz] ] ] ]
				}
			}
			#Let the hate flow through you
			}
			
			#average the bins and append to lists for printing
			
			set upper_leaflet [list]
			set lower_leaflet [list]

			foreach bin $upper_radial_boxes {
				set len [llength $bin]
				if {$len == 0} {lappend upper_leaflet 0 ; continue}
				set sum [expr [join $bin +]]
				set average [expr $sum / [llength $bin] ]
				lappend upper_leaflet $average
			}
			
			foreach bin $lower_radial_boxes {
				set len [llength $bin]
				if {$len == 0} {lappend lower_leaflet 0 ; continue}
				set sum [expr [join $bin +]]
				set average [expr $sum / [llength $bin] ]
				lappend lower_leaflet $average
			}
			puts $upper_outfile "$upper_leaflet"
			puts $lower_outfile "$lower_leaflet"
			flush $upper_outfile
			flush $lower_outfile
						
			if { $i == 0 } {continue}
			if { [expr $i % [expr 10 * $increment]] == 0 } { puts -nonewline "." }
      			if { [expr $i % [expr 500 * $increment]] == 0 } { puts " " }
			
			flush stdout
    			}
		}
		
	main $proteinselection $lipidselection $upper_outfile_name $lower_outfile_name
}

proc helix_rotation {proteinselection lipidselection outfile_name} {
	if 0 {
	proc calculate_unit_vector_to_axis_z {start_com end_com ref_atom} {
		set Rxyz [lindex [$ref_atom get {x y z}] 0]
		set SR [vecsub $Rxyz $start_com]
		set SE [vecsub $end_com $start_com]
		set SN [vecscale $SE [expr [veclength $SR]/[veclength $SE] * [vecdot [vecnorm $SE] [vecnorm $SR ] ] ] ]
		set vector_from_axis [ vecsub $SR $SN ]
		set unit_vector_from_axis [vecnorm $vector_from_axis]
		return [ lindex $unit_vector_from_axis 2 ]
	}
	
	proc calculate_unit_vector_in_xy {start_com end_com ref_atom} {
		set Rxyz [lindex [$ref_atom get {x y z}] 0]
		set SR [vecsub $Rxyz $start_com]
		set SE [vecsub $end_com $start_com]
		set SN [vecscale $SE [expr [veclength $SR]/[veclength $SE] * [vecdot [vecnorm $SE] [vecnorm $SR ] ] ] ]
		set vector_from_axis [ vecsub $SR $SN ]
		set unit_vector_from_axis [vecnorm [list [lindex $vector_from_axis 0] [lindex $vector_from_axis 1] 0 ]]
		return $unit_vector_from_axis
	}
	}
	
	proc calculate_unit_vector_xyz {start_com end_com ref_atom} {
		set Rxyz [lindex [$ref_atom get {x y z}] 0]
		set SR [vecsub $Rxyz $start_com]
		set SE [vecsub $end_com $start_com]
		set SN [vecscale $SE [expr [veclength $SR]/[veclength $SE] * [vecdot [vecnorm $SE] [vecnorm $SR ] ] ] ]
		#set vector_from_axis [ vecsub $SR $SN ]
		set Nxyz [vecadd $SN $start_com]
		set vector_from_axis [vecsub $Rxyz $Nxyz]
		set unit_vector_from_axis [vecnorm $vector_from_axis ]
		return $unit_vector_from_axis
	}
	
	proc uv_start {start_com end_com ref_atom} {
		set Rxyz [lindex [$ref_atom get {x y z}] 0]
		set SR [vecsub $Rxyz $start_com]
		set SE [vecsub $end_com $start_com]
		set SN [vecscale $SE [expr [veclength $SR]/[veclength $SE] * [vecdot [vecnorm $SE] [vecnorm $SR ] ] ] ]
		return [vecadd $SN $start_com]
	}

	
	proc angle_wrap {angle} {
		global M_PI
		if {$angle > $M_PI } {set angle [expr $angle - 2*$M_PI]}
		if {$angle < -$M_PI } {set angle [expr $angle + 2*$M_PI]}
		return $angle
	}
	
	proc make_positive {input_number} {
		if {$input_number < 0 } { return [expr 0 - $input_number] 
		} else { return $input_number } 
	}
	
	proc min_element {numbers} {
		set minimum [lindex $numbers 0]
		set smallest_element 0
		for {set i 1} {$i < [llength $numbers]} {incr i 1} {
			if { [lindex $numbers $i] < $minimum } {set smallest_element $i; set minimum [lindex $numbers $i]}
		}
		return $smallest_element
	}
	
	proc main {proteinselection lipidselection outfile_name} {
		global M_PI 
		set numframes [expr {[molinfo top get numframes]}]
		puts "total number of frames: $numframes"
		set outfile [open $outfile_name w]
		
		set protein [atomselect top $proteinselection ]
		set lipid [atomselect top $lipidselection]
		set all [atomselect top "all" ]
		
		set increment 1
		
		#define helix axis based on com of first/last 3 residues helix backbone
		
		set protein_residues [lsort -integer -unique [[atomselect top $proteinselection] get resid]]
			
		set start [list [lindex $protein_residues 1] [lindex $protein_residues 2] [lindex $protein_residues 3]]
		set end [list [lindex $protein_residues end-3] [lindex $protein_residues end-2] [lindex $protein_residues end-1] ]
		
		set start_sel [atomselect top "name \"B.*\" and resid $start"]
		set end_sel [atomselect top "name \"B.*\" and resid $end"]
		
		#Select residues 1/2 for examining rotation. Residues end-1,end for examining mismatch
		#The lists read 1,2 and end-1,end
		
		set start_res [lindex $protein_residues 2]		
		set end_res [lindex $protein_residues end-2]
		
		set rotation_start [list [atomselect top "name \"B.*\" and resid $start_res"] ]
		set rotation_end [list [atomselect top "name \"B.*\" and resid $end_res"] ]
		
		for {set i 0} {$i < $numframes} {incr i $increment} {
			$all frame $i
			$all update
			$protein frame $i
			$protein update
			$lipid frame $i
			$lipid update
			$start_sel frame $i
			$start_sel update
			$end_sel frame $i
			$end_sel update
			
			foreach atom_r [concat $rotation_start $rotation_end] {
				$atom_r frame $i
				$atom_r update
			}
			
			#I need 2 points and the helical axis
			#if i take the points in the middle of the axis then they might be inaccurate.
			#perhaps a future version should take the start and end residues for ref, to look at mismatch
			
			set start_com [measure center $start_sel]
			set end_com [measure center $end_sel]
			
			set helix_axis [vecsub $end_com $start_com ]
			set helix_angle [expr acos([vecdot [vecnorm $helix_axis] {0 0 1}] ) ]
			
			#Find the matrix to bring the helix_axis to x
			set rotate_helix_to_x [transvecinv $helix_axis]
			
			#Calculate the unit vector of the residue with the helix
			set unit_direction_vector_xyz [calculate_unit_vector_xyz $start_com $end_com [lindex $rotation_start 0] ]
			#puts "Early $unit_direction_vector_xyz"
			
			#Rotate the unit vector so that it is moving around the x axis
			set vector_around_x [vectrans $rotate_helix_to_x $unit_direction_vector_xyz]
			#puts "Length vx [veclength $vector_around_x]"
			
			#The reference vector is a unit vector in z
			set unit_vector_z {0 0 1}
			
			#The angle of these two is the angle of interest- without directionality
			set rad_rotation [expr acos([vecdot $unit_vector_z $vector_around_x])]
			
			#puts "Rad $rad_rotation"
			
			#Calculate which side the angle is using the cross product & modify rad_rotation accordingly
			#set side_modifier [lindex [veccross $vector_around_x $unit_vector_z] 1]
			#puts "Mod [veccross $vector_around_x $unit_vector_z]"
			if {[lindex $vector_around_x 1] < 0} { set side_modifier -1 } else { set side_modifier 1}
			set rad_rotation [expr $side_modifier * $rad_rotation]
			
			#Determine if the start is higher or lower than the lipid center- if lower, rotate so that 0 always points to the upper leaflet	
			set start_comz [lindex $start_com 2]
			set lipid_comz [lindex [measure center $lipid] 2]
			if { $start_comz < $lipid_comz} {set rad_rotation [angle_wrap [expr $rad_rotation + $M_PI]]}
			
			#Get rotation in degrees
			#puts "Side mod $side_modifier"
			#puts "New Rad $rad_rotation"
			#puts "Pi $M_PI"
			#puts "Rad/Pi [expr ($rad_rotation/$M_PI)]"
			puts "Frame $i"
			set rotation [expr ($rad_rotation/$M_PI)*180]
			
			#print output to terminal and file, then continue
			puts $outfile "$i $rotation"
			#if {$noideality < 30} {puts $outfile "$i $rotation"}
			
			flush $outfile
			
			#Debug
			#
			if 0 {
			#graphics top line $start_com [vecadd $start_com $vector_around_x]
			graphics top color green
			#graphics top line $start_com $end_com
			graphics top color blue
			graphics top line [uv_start $start_com $end_com [lindex $rotation_start 0] ] [vecadd [uv_start $start_com $end_com [lindex $rotation_start 0]] $unit_direction_vector_xyz]
			#graphics top point $start_com
			graphics top color red
			#graphics top point [lindex [[lindex $rotation_start 0] get {x y z}] 0 ]
			#graphics top line [uv_start $start_com $end_com [lindex $rotation_start 0] ] [vecadd [uv_start $start_com $end_com [lindex $rotation_start 0]] $unit_direction_vector_xyz]
			graphics top line [uv_start $start_com $end_com [lindex $rotation_start 0] ] [lindex [[lindex $rotation_start 0] get {x y z}] 0 ]
			
			graphics top line {-31.2 22.5 21.2} [vecadd {-31.2 22.5 21.2} $vector_around_x]
			graphics top color green
			graphics top line {-31.2 22.5 21.2} {-31.2 22.5 22.2}
			
			puts "Rotated Vector $vector_around_x"
			puts "Rotation $rotation"

			}
			if 0 {
			if { $i == 0 } {continue}
			if { [expr $i % [expr 10 * $increment]] == 0 } { puts -nonewline "." }
      			if { [expr $i % [expr 500 * $increment]] == 0 } { puts " " }
			}
			
			flush stdout
    			}
		}
	
	main $proteinselection $lipidselection $outfile_name
}
