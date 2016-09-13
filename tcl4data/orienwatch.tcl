set filename "50-0.03testact-pol.lammpstrj"
set color "red"
set sel "type 1"
set nf 10
set cradius 0.08
set dradius 0.20
set start   0
set end     100

mol delete all
set molid [mol new $filename first $start last $end step $nf waitfor all]
proc draw_orientation {args} {
    global molid sel color cradius dradius
    global vmd_frame

    set filled yes
    graphics $molid delete all

    set now $vmd_frame($molid); set next [expr $now+1]
    set all [atomselect $molid all]
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
        
        set start {}
	set end   {}
        lappend start $now_x
        lappend start $now_y
        lappend start $now_z
        lappend end   $next_x
        lappend end   $next_y
        lappend end   $next_z
		set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
		graphics $molid cylinder $start $middle radius $cradius resolution 30 filled $filled
		graphics $molid cone $middle $end radius $dradius resolution 30
    }
    $atomlist delete
}

trace variable vmd_frame($molid) w draw_orientation
