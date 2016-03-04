set runfile mem
set runfolder ./deal_data
set vmdpath  ~/VMD/vmd-1.9.1/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
set namdpath ~/NAMD_2.10b1_Source/Linux-x86_64-g++/namd2

mol delete all
mol new ./$runfile.psf
mol addfile ./result/$runfile.dcd first 0 last 8000 waitfor all
file mkdir $runfolder

molinfo top set frame 1200
set carbon [atomselect top "segname CT1"]
set mes [measure minmax $carbon]
set lowerEndZ [lindex $mes 0 2]
set upperEndZ [lindex $mes 1 2]
set lowerEndY [lindex $mes 0 1]
set upperEndY [lindex $mes 1 1]

###################################################################

set water_energy  [atomselect top "water and (z>$lowerEndZ and z<$upperEndZ)"]
set carbon_energy [atomselect top "segname CT1"] 
source $vmdpath

set output_flie [open $runfolder/namdenergy.dat w]
for {set i $lowerEndZ} {$i<=$upperEndZ} {incr i $precision} {
		set data_file [open $runfolder/namdenergy_{$i}.dat r]
		set j [expr $i+$precision]
		set sum 0
		set numfra 1
		set water_energy  [atomselect top "water and (z>$i and z<$j)"]
		namdenergy -vdw -sel $water_energy -ofile $data_file -par ./par/par_all27_prot_lipid.prm -exe $namdpath
		while {[gets $data_file  line] >= 1} {
      			set tempdata [string range $line 28 37]
      			if {[string range $line 0 5] >= 500} {
      				set sum [expr $sum+$tempdata]
      				set numfra [expr $numfra+1]
      			}
  		 	}
  		set average [expr 1.0*$sum/$numfra]
  		puts $output_flie "$i $average"
		}
		close $data_file 
}

close $output_flie
exit
