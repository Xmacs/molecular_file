mol delete all
mol new all_f_07.lammpstrj first 0 last 10 waitfor all
proc draw_orientation {args} {
    global vmd_frame

    set res 6
    set gidlist {}
    set filled yes
    graphics 0 delete all

    set now $vmd_frame(0); set next [expr $now+1]
    set all [atomselect 0 all]
    set box_length [expr [lindex [measure minmax $all] 1 0]-[lindex [measure minmax $all] 0 0] ]
    set box_width  [expr [lindex [measure minmax $all] 1 1]-[lindex [measure minmax $all] 0 1] ]
    set box_heigh  [expr [lindex [measure minmax $all] 1 2]-[lindex [measure minmax $all] 0 2] ]
    $all delete

    lappend gidlist [graphics 0 color red]
    set atomlist [[atomselect 0 "type 1"] get index]
    foreach loop $atomlist {
        
        ############ we need deal with the PBC coordinate############# refer to http://simulation.haotui.com/viewthread.php?tid=39109&amp;amp;jdfwkey=olbei2
        set start [lindex [[atomselect 0 "index $loop" frame $now]  get {x y z}] 0]
        set end   [lindex [[atomselect 0 "index $loop" frame $next] get {x y z}] 0]
        set now_x  [[atomselect 0 "index $loop" frame $now]   get x]
        set now_y  [[atomselect 0 "index $loop" frame $now]   get y] 
        set now_z  [[atomselect 0 "index $loop" frame $now]   get z] 
        set next_x [[atomselect 0 "index $loop" frame $next]  get x]
        set next_y [[atomselect 0 "index $loop" frame $next]  get y] 
        set next_z [[atomselect 0 "index $loop" frame $next]  get z]
        set nx 0; set ny 0; set nz 0
        set dx [lindex [vecsub $end $start] 0]
        if {$dx > [expr $box_length / 2]} {incr nx -1}
        if {$dx < [expr -1 * $box_length / 2]} {incr nx}
        set dy [lindex [vecsub $end $start] 1]
        if {$dy > [expr $box_width / 2]} {incr ny -1}
        if {$dy < [expr -1 * $box_width / 2]} {incr ny}
        if {$dz > [expr $box_heigh / 2]} {incr nz -1}
        if {$dz < [expr -1 * $box_heigh / 2]} {incr nz}
        set next_x [expr $nx*$box_length+$next_x];set next_y [expr $next_y+$ny*$box_width];set next_z [expr $next_z+$nz*$box_heigh]
        
        set start {}
		set end   {}
        lappend start $now_x
        lappend start $now_y
        lappend start $now_z
        lappend end   $next_x
        lappend end   $next_y
        lappend end   $next_z
		set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
		graphics 0 cylinder $start $middle radius 0.2
		graphics 0 cone $middle $end radius 0.5
        #puts "fuck!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        
    }
}

trace variable vmd_frame(0) w draw_orientation
