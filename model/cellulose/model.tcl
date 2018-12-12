set degree  20        ;#degree of polymerization
set branch_degree "random" ;#half all none random
file mkdir  temp
set bondlength 1.5


####需要一个整体的base pdb生成文件####
####################################
proc veccut {coor_c coor_n} {
        global bondlength
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

proc RandomRange { min max } { 
    # 获得[0.0,1.0)之间的随机数 
    set rd [expr rand()] 
    # 将$rd放大到[$min, $max) 
    set result [expr $rd * ($max - $min) + $min] 
    return $result 
} 

# #FUNC:获取[min, max)区间是随机整数 # 
proc RandomRangeInt { min max } { 
    return [expr int([RandomRange $min $max])] 
}


proc CelluloseGenbranch { degree string } {
    file mkdir branch_temp
    mol delete all
    mol load psf chain.psf pdb chain.pdb
    set all [atomselect top all]
    $all writepdb branch_temp/chain_0.pdb
    $all writepsf branch_temp/chain_0.psf
    $all delete
    mol delete all

    if {$string == "half"} {

        set numlist {}
        for {set i 2} { $i <= [expr $degree*2] } {set i [expr $i+2]} {lappend numlist $i}

    } elseif {$string == "all"} {

        set numlist {}
        for {set i 1} { $i <= [expr $degree*2] } {incr i} {lappend numlist $i}

    } elseif {$string == "random"} {

        set numlist {}
        for {set i 1} { $i <= $degree } {incr i} {lappend numlist [ RandomRangeInt 1 [expr $degree*2] ] }
        set numlist [lsort -integer $numlist]
        set numlist [lsort -unique $numlist]

    }

    foreach i $numlist {
        set chain   [mol load psf branch_temp/chain_0.psf pdb branch_temp/chain_0.pdb]
        set ethanol [mol load psf Ethylene_glycol_use.psf pdb Ethylene_glycol_use.pdb]
        set resid_atom    [atomselect $chain "resid $i and name O6"]
        set ethanol_all   [atomselect $ethanol "all"]
        set ethanol_atom  [atomselect $ethanol "name C1"]
        set mid [lindex [$resid_atom   get {x y z}] 0 ]
        set pol [lindex [$ethanol_atom get {x y z}] 0 ]
        set movelist    [veccut $mid $pol]
        if {[expr $i%2] == 0 } {
            $ethanol_all move [trans center $mid bond  {0 0 0} $pol 180];#$coor_atom_cb $coor_atom_n2 $pos_pol 10]
        }
        $ethanol_all moveby    $movelist
        $ethanol_all set segname A1
        $ethanol_all set resid  [expr $degree*10+$i]
        $ethanol_all writepdb  branch_temp/ethanol_$i.pdb
        $ethanol_all writepsf  branch_temp/ethanol_$i.psf
        $resid_atom     delete
        $ethanol_all    delete
        $ethanol_atom   delete
    }

    package require psfgen
    resetpsf
    readpsf  chain.psf
    coordpdb chain.pdb
    foreach i $numlist {
        readpsf  branch_temp/ethanol_$i.psf
        coordpdb branch_temp/ethanol_$i.pdb
    }
    writepsf branch_temp/simulation.psf
    writepdb branch_temp/simulation.pdb

    mol delete all
    package require psfgen
    resetpsf
    topology top_all36_branch.rtf
    segment A1 {
        pdb branch_temp/simulation.pdb
    }
    foreach i $numlist {
        patch 20cc A1:$i A1:[expr $degree*10+$i]
    }
    set patch_line [expr $degree*2-1]
    for {set i 1} {$i <= $patch_line} {incr i} {
        set j [expr $i+1]
        patch 14bb A1:$j A1:$i
    }
    regenerate angles dihedrals
    coordpdb branch_temp/simulation.pdb A1
    guesscoord
    writepsf chain_branch.psf
    writepdb chain_branch.pdb

    #file delete -force branch_temp

}


######生成单条纤维素链坐标#######################################
mol new base.pdb
set selCore [ atomselect top all]
set i 0
### BUILD CORE #####
set pdb [format "temp/copy%-0.4i" $i]
#set pdb copy$i
set coreveclist [vecinvert [measure center $selCore]]
$selCore moveby $coreveclist
set clist {}
    lappend clist [expr $i*10]
    lappend clist 0.0
    lappend clist 0.0
$selCore moveby $clist
$selCore set chain A
$selCore set segname A1
# check - puts $pdb
$selCore writepdb $pdb.pdb
unset pdb
set all [atomselect top all]
set veclist [list [measure minmax $all]]
set xmin [lindex $veclist 0 0 0]
set xmax [lindex $veclist 0 1 0]
set ymin [lindex $veclist 0 0 1]
set ymax [lindex $veclist 0 1 1]
set zmin [lindex $veclist 0 0 2]
set zmax [lindex $veclist 0 1 2]
set delta_z [expr abs($zmin)+abs($zmax)]
set delta_x [expr abs($xmin)+abs($xmax)]
set delta_y [expr abs($ymin)+abs($ymax)]
$selCore delete
$all     delete
mol delete all

package require topotools
for {set i 1} {$i <= $degree } {incr i} {
        set pdb     [format "temp/copy%-0.4i" $i]
        set j       [expr $i-1]
        set oldpdb  [format "temp/copy%-0.4i" $j]

        mol new $oldpdb.pdb
        set selCore [ atomselect top "resid 1 2"]
        ### BUILD CORE #####
        
        #set pdb copy$i
        set coreveclist [vecinvert [measure center $selCore]]
        $selCore moveby $coreveclist
        set clist {}
            lappend clist 0.0
            lappend clist 0.0
            lappend clist [expr $i*($delta_z+1.0)]
        $selCore moveby $clist
        set resid1 [atomselect top "resid 1"]
        set resid2 [atomselect top "resid 2"]
        $resid1 set resid [expr $i*2+1]
        $resid2 set resid [expr $i*2+2]
        set selCore [ atomselect top "resid [expr $i*2+1] [expr $i*2+2]"]
        $selCore set chain A
        #set SNAME [format "A%-0.1i" $i]
        $selCore set segname A1
        # check - puts $pdb
        $selCore writepdb $pdb.pdb
        $resid1 delete
        $resid2 delete
        $selCore delete
        mol delete all

        set first_sel  [mol new $pdb.pdb]
        set second_sel [mol new $oldpdb.pdb]
        set sel_first  [atomselect $first_sel  "all"]
        set sel_second [atomselect $second_sel "all"]
        set merge_sel  [::TopoTools::selections2mol "$sel_first $sel_second"]
        animate write pdb $pdb.pdb $merge_sel
        $sel_second delete
        $sel_first  delete
        mol delete all

        unset pdb
        unset oldpdb
}
mol delete all

set pdb [format "temp/copy%-0.4i" $degree]
mol new $pdb.pdb
set all [atomselect top "not resid 1 2"]
$all writepdb tmp.pdb
mol delete all
mol new tmp.pdb
for {set i 3} {$i <= [expr $degree*2+2]} {incr i} {
    set resid_sel [atomselect top "resid $i"]
    $resid_sel set resid [expr $i-2]
    $resid_sel delete
}
set all [atomselect top all]
$all writepdb tmp.pdb
mol delete all
file delete -force temp
#############纤维素坐标代码末端#############################
#exit

############生成单条纤维素结构文件##############################
package require psfgen
resetpsf
topology top_all36_carb.rtf
segment A1 {
    pdb tmp.pdb
}
set patch_line [expr $degree*2-1]
for {set i 1} {$i <= $patch_line} {incr i} {
    set j [expr $i+1]
    patch 14bb A1:$j A1:$i
}
    
regenerate angles dihedrals
coordpdb tmp.pdb A1
guesscoord

writepsf chain.psf
writepdb chain.pdb
#################################################################

######################  branch 接枝过程   ########################################
mol load psf chain.psf pdb chain.pdb
set selCore [atomselect top all]
set coreveclist [vecinvert [measure center $selCore]]
$selCore moveby $coreveclist
$selCore writepdb chain.pdb
$selCore delete
mol delete all

if {$branch_degree == "none"} {
    puts "Great!!!"
} elseif {$branch_degree == "half"} {
    CelluloseGenbranch $degree "half"
} elseif {$branch_degree == "all"} {
    CelluloseGenbranch $degree "all"
} elseif {$branch_degree == "random"} {
    CelluloseGenbranch $degree "random"
}

#################genrate the gromacs file ###########
mol delete all
mol load psf chain_branch.psf pdb chain_branch.pdb
set all [atomselect top all]
set sizelist [measure minmax $all]
set size_x   [expr abs([lindex $sizelist 0 0])+abs([lindex $sizelist 1 0])+200.0]
set size_y   [expr abs([lindex $sizelist 0 1])+abs([lindex $sizelist 1 1])+200.0]
set size_z   [expr abs([lindex $sizelist 0 2])+abs([lindex $sizelist 1 2])+200.0]
set sizelist {}
    lappend sizelist $size_x
    lappend sizelist $size_y
    lappend sizelist $size_z
pbc set $sizelist -all -molid top
pbc box -center origin -shiftcenter {0 0 0}
set all [atomselect top all]
$all writepdb chain_branch.pdb
package require topotools
topo writegmxtop use.top [list par/par_all36_branch.prm]
set all [atomselect top all]
$all writegro use.gro

exit
