set filename "50-0.03testact-pol.lammpstrj"
set color "red"
set sel "type 1"
set nf 10
set cradius 0.08
set dradius 0.20

mol delete all
set molid [mol new $filename first 0 last $nf waitfor all]
proc draw_orientation {args} {
    global molid sel color cradius dradius
    global vmd_frame

    set filled yes
    graphics 0 delete all

    set now $vmd_frame($molid); set next [expr $now+1]
    set all [atomselect $molid all]
    #set box_length [expr [lindex [measure minmax $all] 1 0]-[lindex [measure minmax $all] 0 0] ]
    #set box_width  [expr [lindex [measure minmax $all] 1 1]-[lindex [measure minmax $all] 0 1] ]
    #set box_heigh  [expr [lindex [measure minmax $all] 1 2]-[lindex [measure minmax $all] 0 2] ]
    $all delete

    graphics $molid color red
    set atomlist [[atomselect $molid $sel] get index]
    foreach loop $atomlist {
        
        ############ we need deal with the PBC coordinate############# 
        set start [lindex [[atomselect $molid "index $loop" frame $now]  get {x y z}   ] 0]
        set end   [lindex [[atomselect $molid "index $loop" frame $now]  get {vx vy vz}] 0]
        set now_x  [[atomselect $molid "index $loop" frame $now]   get x ]
        set now_y  [[atomselect $molid "index $loop" frame $now]   get y ] 
        set now_z  [[atomselect $molid "index $loop" frame $now]   get z ] 
        set next_x [[atomselect $molid "index $loop" frame $now]   get vx]
        set next_y [[atomselect $molid "index $loop" frame $now]   get vy] 
        set next_z [[atomselect $molid "index $loop" frame $now]   get vz]
        set next_x [expr $now_x+$next_x];set next_y [expr $now_y+$next_y];set next_z [expr $now_z+$next_z]
        if {0} {
        set nx 0; set ny 0; set nz 0
        set dx [lindex [vecsub $end $start] 0]
        if {$dx > [expr $box_length / 2]} {incr nx -1}
        if {$dx < [expr -1 * $box_length / 2]} {incr nx}
        set dy [lindex [vecsub $end $start] 1]
        if {$dy > [expr $box_width / 2]} {incr ny -1}
        if {$dy < [expr -1 * $box_width / 2]} {incr ny}
        #if {$dz > [expr $box_heigh / 2]} {incr nz -1}
        #if {$dz < [expr -1 * $box_heigh / 2]} {incr nz}
        set next_x [expr $nx*$box_length+$next_x];set next_y [expr $next_y+$ny*$box_width];#set next_z [expr $next_z+$nz*$box_heigh]
        }
        
        set start {}
		set end   {}
        lappend start $now_x
        lappend start $now_y
        lappend start $now_z
        lappend end   $next_x
        lappend end   $next_y
        lappend end   $next_z
		set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
		graphics $molid cylinder $start $middle radius $cradius resolution 20 filled $filled
		graphics $molid cone $middle $end radius $dradius resolution 20
    }
    $atomlist delete
}

trace variable vmd_frame($molid) w draw_orientation
