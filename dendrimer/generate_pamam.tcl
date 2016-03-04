proc veccut {coor_c coor_n} {
		global PI bondlength
		set x [lindex $coor_c 0]
		set y [lindex $coor_c 1]
		set z [lindex $coor_c 2]
		set a [lindex $coor_n 0]
		set b [lindex $coor_n 1]
		set c [lindex $coor_n 2]
		if {$a>=0 && $b>=0} {
			set a [expr $a+$bondlength]
			set b [expr $b+$bondlength]
		} elseif {$a<0 && $b>0} {
			set a [expr $a-$bondlength]
			set b [expr $b+$bondlength]
		} elseif {$a<0 && $b<0} {
			set a [expr $a-$bondlength]
			set b [expr $b-$bondlength]
		} elseif {$a>0 && $b<0} {
			set a [expr $a+$bondlength]
			set b [expr $b-$bondlength]
		}
		#set a [expr $a+$bondlength*cos((45*$PI)/180)]
		#set b [expr $b+$bondlength*sin((45*$PI)/180)]
		set dx [expr $x-$a]
		set dy [expr $y-$b]
		set dz [expr $z-$c]
		set list {}
		lappend list $dx
		lappend list $dy
		lappend list $dz
		return $list
}

mol delete all
set unitdir [ file join [pwd] units ]
file mkdir $unitdir
set PI 3.1415926535897931
set bondlength 1.5
set TOPOLOGY  ./top_pam_f_amber.top
set newfile   ./unit
set core_atom_tail_name(1) "name CC1"
set core_atom_tail_name(2) "name CC2"
set core_atom_tail_name(3) "name CC3"
set core_atom_tail_name(4) "name CC4"
set bru_atom_head_name "name C2"
set bru_atom_mid_name  "name C3"
set bru_atom_tail_name(0) "name CB1"
set bru_atom_tail_name(1) "name CB2"
set ter_atom_head_name "name CT2"
set nGen        2           ; ###  gen number
set core_id [mol load pdb core.pdb]
set repeat_id [mol load pdb repeat.pdb]
set terminal_id [mol load pdb f-ter.pdb]
set selCore [atomselect $core_id all]
set selRep  [atomselect $repeat_id  all]
set selTerm [ atomselect $terminal_id all]
	
	####### creat core#######
	set coreveclist [vecinvert [measure center $selCore]]
	$selCore moveby $coreveclist
	set resc 1
	$selCore set resid $resc
	$selCore set resname COR
	$selCore set segname PAM
	$selCore writepdb ./units/resid_1.pdb
	#######   get core head coordinate #####
	
	########coordinate the repeat units ###########
	for {set i 0} { $i <= [expr $nGen-1] } { incr i } {
		if {$i==0} {
			mol load pdb ./units/resid_1.pdb
			set molid [mol load pdb ./units/resid_1.pdb]
			for { set d 2 } { $d<=5 } { incr d } {
				set j [expr $d-1]
				set coor_core_atom [atomselect $molid $core_atom_tail_name($j)]
				set coor_branch_atom [atomselect $repeat_id $bru_atom_head_name]
				set mid [vecinvert [lindex [$coor_branch_atom get {x y z}] 0 ]]
				set pol [vecinvert [lindex [$coor_core_atom get {x y z}] 0 ]]
				set movelist [veccut $mid $pol]
				$selRep moveby $movelist
				#######rotate the molecular############
				if {1} {
				if {$d<4} {
					set pos_mid [lindex [$coor_branch_atom get {x y z}] 0 ]
					set pos_pol [lindex [$coor_core_atom get {x y z}] 0 ]
					$selRep move [trans center $pos_mid bond {0 0 0} $pos_pol 180]
				} else {
					set pos_mid [lindex [$coor_branch_atom get {x y z}] 0 ]
					set pos_pol [lindex [$coor_core_atom get {x y z}] 0 ]
					$selRep move [trans center $pos_mid bond  {0 0 0} $pos_pol 180]
				}
				}
				$selRep set resid $d  
				$selRep set resname BRU
				$selRep set segname PAM
				$selRep writepdb ./units/resid_$d.pdb
			}
		} else {
			set star_resid [expr pow(2,$i+2)-2]
			set star_resid [expr int($star_resid)]
			set end_resid [expr pow(2,$i+3)-3]
			set end_resid [expr int($end_resid)]
			for {set loop $star_resid} { $loop <= $end_resid } { incr loop 2 } {
				set prevresid [expr $loop/2-1]
				set prevresid [expr int($prevresid)]
				#mol load pdb ./units/resid_$prevresid.pdb
				set tempresid [mol load pdb ./units/resid_$prevresid.pdb]
				set coor_tail_atom(0) [atomselect $tempresid $bru_atom_tail_name(0)]
				set coor_tail_atom(1) [atomselect $tempresid $bru_atom_tail_name(1)]
				for {set temp 0} {$temp<=1} { incr temp } {
					set coor_branch_atom [atomselect $repeat_id $bru_atom_head_name]
					set mid [vecinvert [lindex [$coor_branch_atom get {x y z}] 0 ]]
					set pol [vecinvert [lindex [$coor_tail_atom($temp) get {x y z}] 0 ]]
					set movelist [veccut $mid $pol]
					$selRep moveby $movelist
					if {1} {
					set pos_mid [lindex [$coor_branch_atom get {x y z}] 0 ]
					set pos_pol [lindex [$coor_tail_atom($temp) get {x y z}] 0 ]
					$selRep move [trans center $pos_mid bond  {0 0 0} $pos_pol 270];#$coor_atom_cb $coor_atom_n2 $pos_pol 10]
					}
					set resc [expr $loop+$temp]
					$selRep set resid $resc  
					$selRep set resname BRU
					$selRep set segname PAM
					$selRep writepdb ./units/resid_$resc.pdb
				}
			}
		}

	}
	####################################################################
	
	##############coordinate the terminal units ##################
	set star_resid [expr pow(2,$nGen+2)-2]
	set star_resid [expr int($star_resid)]
	set end_resid [expr pow(2,$nGen+3)-3]
	set end_resid [expr int($end_resid)]
	for {set loop $star_resid} { $loop <= $end_resid } { incr loop 2 } {
		set prevresid [expr $loop/2-1]
		set prevresid [expr int($prevresid)]
		#mol load pdb ./units/resid_$prevresid.pdb
		set tempresid [mol load pdb ./units/resid_$prevresid.pdb]
		set coor_tail_atom(0) [atomselect $tempresid $bru_atom_tail_name(0)]
		set coor_tail_atom(1) [atomselect $tempresid $bru_atom_tail_name(1)]
		for {set temp 0} {$temp<=1} { incr temp } {
			set branch_atom [atomselect $terminal_id $ter_atom_head_name]
			set mid [vecinvert [lindex [$branch_atom get {x y z}] 0 ]]
			set pol [vecinvert [lindex [$coor_tail_atom($temp) get {x y z}] 0 ]]
			set movelist [veccut $mid $pol]
			$selTerm moveby $movelist
			if {1} {
			set pos_mid [lindex [$branch_atom get {x y z}] 0 ]
			set pos_pol [lindex [$coor_tail_atom($temp) get {x y z}] 0 ]
			$selTerm move [trans center $pos_mid bond  {0 0 0} $pos_pol 180];#$coor_atom_cb $coor_atom_n2 $pos_pol 10]
			}
			set resc [expr $loop+$temp]
			#set pdb [format "copy%-0.4i" $resc]
			$selTerm set resid $resc  
			$selTerm set resname FER
			$selTerm set segname PAM
			$selTerm writepdb ./units/resid_$resc.pdb
		}
	}

	###################creat psf pdb file  #####
	mol delete all
	resetpsf
	package require psfgen
	topology $TOPOLOGY
	for { set id 1} { $id<=$resc} {incr id 1} {
		segment PA$id {
			pdb ./units/resid_$id.pdb
		}
		coordpdb ./units/resid_$id.pdb PA$id
	}
	guesscoord
	writepsf dendrimer.psf
	writepdb dendrimer.pdb
	mol delete all
	mol load psf dendrimer.psf pdb dendrimer.pdb
	set all [atomselect top all]
	$all set segname PAM
	$all writepdb dendrimer.pdb
	$all writepsf dendrimer.psf
	mol delete all
	
	resetpsf
	package require psfgen
	#topology $TOPOLOGY
	segment PAM {
		pdb dendrimer.pdb
	}
	puts "fuck !!!!!"
	patch LIN1 PAM:2 PAM:1
	patch LIN2 PAM:3 PAM:1
	patch LIN3 PAM:4 PAM:1
	patch LIN4 PAM:5 PAM:1
	for { set id 6} { $id<=[expr pow(2, ($nGen-1)+3 )-3]} {incr id 2} {
		set id [expr int($id)]
		set reid [expr $id/2-1]
		set secid [expr $id+1]
	    #set reid [expr int($reid)]
		#set secid [expr int($secid)]
		patch CON1 PAM:$id PAM:$reid
		patch CON2 PAM:$secid PAM:$reid
	}
	##########patch ternimal#####
	for { set terid [expr pow(2,$nGen+2)-2] } { $terid<=$resc} {incr terid 2} {
		set terid [expr int($terid)]
		set reid [expr $terid/2-1]
		set secid [expr $terid+1]
	    set reid [expr int($reid)]
		set secid [expr int($secid)]
		patch TEN1 PAM:$terid PAM:$reid
		patch TEN2 PAM:$secid PAM:$reid
	}
	coordpdb dendrimer.pdb PAM
	puts "fuck !!!!!"
	guesscoord
	writepsf g.psf
	writepdb g.pdb
	mol delete all
	mol new g.psf
	mol addfile g.pdb
	set all [atomselect top all]
	$all writemol2 g.mol2
	
