proc ::OrientationWatch::draw_orientation {mol sel {color red} {scale 1.0} {radius 0.2} {orientationidx 0}} {
    variable orientationvalue
    variable orientationcenter

    set res 6
    set gidlist {}
    set filled yes

    # perhaps this should use the center information
    if {[catch {measure center $sel weight mass} center]} {
        if {[catch {measure center $sel} center]} {
            puts stderr "problem computing orientation center: $center"
            return {}
        }
    }
    set orientationvalue($orientationidx) [format "%6.2f D" [veclength $vector]]
    set vechalf [vecscale [expr $scale * 0.5] $vector]
    set now [molinfo $mol get frame]
    set next [expr $now+1]

    lappend gidlist [graphics $mol color $color]
    set atomlist [[atomselect $mol $sel] get index]
    foreach loop $atomlist {
        set start [lindex [[atomselect $mol "index $loop" frame $now]  get {x y z}] 0]
		set end   [lindex [[atomselect $mol "index $loop" frame $next] get {x y z}] 0]
		set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
		lappend gidlist [graphics $mol cylinder $start $middle radius 0.05]
		lappend gidlist [graphics $mol cone $middle $end radius 0.08]
        
    }

    return $gidlist
}