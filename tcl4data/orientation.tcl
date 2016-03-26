mol delete all
mol new all_f_07.lammpstrj first 0 last 10 waitfor all

#proc draw_orientation {args} {
    #global vmd_frame

    set res 6
    set gidlist {}
    set filled yes
    graphics 0 delete all

    set now $vmd_frame(0); set next [expr $now+1]

    lappend gidlist [graphics 0 color red]
    set atomlist [[atomselect 0 "type 1"] get index]
    foreach loop $atomlist {
        
        ############ we need deal with the PBC coordinate#############
        set now_x  [[atomselect 0 "index $loop" frame 1]  get x]
        set now_y  [[atomselect 0 "index $loop" frame 1]  get y] 
        set now_z  [[atomselect 0 "index $loop" frame 1]  get z] 
        set next_x [[atomselect 0 "index $loop" frame 2]  get x]
        set next_y [[atomselect 0 "index $loop" frame 2]  get y] 
        set next_z [[atomselect 0 "index $loop" frame 2]  get z]
        if {}
        
        set start {}
		set end   {}
        lappend start $now_x
        lappend start $now_y
        lappend start $now_z
        lappend end   $next_x
        lappend end   $next_y
        lappend end   $next_z
		set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
		graphics 0 cylinder $start $middle radius 0.05
		graphics 0 cone $middle $end radius 0.08
        puts "fuck!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        
    }
#}

#trace variable vmd_frame(0) w draw_orientation
