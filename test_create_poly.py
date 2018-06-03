proc collision_tcl { i posx posy posz shield box_length } {
# Checks whether particle <i> at coordinates (<posx>, <posy>, <posz>) collides
# with any other particle j < i due to a minimum image distance smaller than <shield>.
    set min_dist [mindist $posx $posy $posz]
    if { [string is double $min_dist] }
    {
    if { $min_dist > $shield }
     { return 0 }
     else { return 1 }
    } else {
    set min_dist 30000
    for {set j 0} {$j<$i} {incr j} {
        set posj [part $j print pos]; set posjx [lindex $posj 0]; set posjy [lindex $posj 1]; set posjz [lindex $posj 2]
        set dx [expr ($posx-$posjx) - round(($posx-$posjx)/$box_length)*$box_length]
        set dy [expr ($posy-$posjy) - round(($posy-$posjy)/$box_length)*$box_length]
        set dz [expr ($posz-$posjz) - round(($posz-$posjz)/$box_length)*$box_length]
        set min_dist [min $min_dist [expr $dx*$dx + $dy*$dy + $dz*$dz]]
    }
    if { $min_dist > [sqr $shield] } { return 0 } else { return 1 }
    }
}



def collision_analyze(i, pos_, shield, box_):
    """Checks"""
    min_dist = min_dist [mindist $posx $posy $posz]

    if string is double min_dist:
        if min_dist > shield:
            return 0
        else:
            return 1
    else:
        min_dist = 30000
        for j in range(i):
            posj_ = position of particle j
            for _ in range(3):
                dx_[_] = (pos_[_]-posj_[i]) - round((pos_[_]-posj_[i])/box_[_])*box_[_]
            set min_dist = min(min_dist, np.einsum('i,i', dx_, dx_))

    return None
