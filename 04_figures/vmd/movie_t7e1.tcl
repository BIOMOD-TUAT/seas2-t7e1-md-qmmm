# movie_t7e1_v5.tcl  --  VMD script (run on your own PC, in the VMD Tk Console)
#
# v5: optional opening shot of the whole T7E1-DNA complex that zooms into the active site.
#     Needs rqB_complex_R.pdb (protein + DNA + Mn of the reactant, same coordinates) in the same folder.
#     Without that file the script behaves like v4.
# v4: no dashed lines (set DASHES 1 to bring them back). A broken bond simply disappears.
#
# Changes from v2 (kept from v3)
#   - closer zoom on the reaction site, smaller Mn spheres (they no longer hide the P atom)
#   - the DNA strand is coloured in two halves: BLUE = 3' side of the cut, ORANGE = 5' side of the cut
#   - the broken bond is drawn as thick red dashes (v2: thin dotted line)
#   - the distance read-out is placed at the upper left of the view, away from the atoms,
#     and says "bonded" / "BROKEN" / "forming"
#   - fewer surrounding atoms (within 5 A), far atoms fade into the background
#
# Usage (VMD Tk Console):
#   cd "path/to"
#   source movie_t7e1_v5.tcl
#   ... rotate / zoom with the mouse ...   (after rotating, move the frame slider once to refresh the text)
#   make_movie
# PowerShell, same folder:
#   ffmpeg -framerate 24 -i frames/f%04d.tga -c:v libx264 -pix_fmt yuv420p -crf 18 t7e1_cleavage.mp4

set PDB     "rqB_movie_traj.pdb"
set OUTDIR  "frames"
set BONDCUT 2.00    ;# P-O distance (A) below which the bond counts as "bonded"
set ZOOM    5.0     ;# initial zoom factor
set LABEL_UP    4.2 ;# position of the read-out, in A, measured on the screen from the reaction centre
set LABEL_RIGHT -7.0
set LABEL_STEP  0.9
set LABEL_SIZE  1.2
set OVPDB   "rqB_complex_R.pdb"  ;# whole complex (optional)
set ZOOM_OUT 14.0   ;# how much smaller the opening shot is than the close-up
set N_HOLD  48      ;# frames of the still opening shot (24 = 1 s)
set N_ZOOM  96      ;# frames of the zoom into the active site
set HIDE_OV_IN_CLOSEUP 1  ;# 1 = hide the whole-complex drawing during the reaction, 0 = keep it
set SHOW_MARKS 1
set DASHES 0        ;# 1 = draw dashed lines for a broken / not-yet-formed bond, 0 = draw nothing

mol delete all
mol new $PDB type pdb waitfor all
set MID [molinfo top]

display projection Orthographic
display resize 1280 720
axes location Off
color Display Background black
display depthcue on
display cuemode Linear
display ambientocclusion off
display shadows off

mol delrep 0 $MID

proc addrep {sel style colour {material Opaque}} {
    global MID
    mol representation {*}$style
    mol selection $sel
    mol color {*}$colour
    mol material $material
    mol addrep $MID
}

# ---- atoms of the breaking and the forming bond -------------------------------------------
set NF [molinfo $MID get numframes]
set sP  [atomselect $MID "beta > 1.5 and beta < 2.5 and name P"]
set sO3 [atomselect $MID "beta > 1.5 and beta < 2.5 and name \"O3'\""]
if {[$sP num] != 1 || [$sO3 num] != 1} {
    puts "WARNING: expected exactly one P and one O3' in the reaction core (found [$sP num] and [$sO3 num])"
}
set sW [atomselect $MID "beta > 1.5 and beta < 2.5 and name OH2" frame [expr {$NF - 1}]]
$sP frame [expr {$NF - 1}]
set pP [lindex [$sP get {x y z}] 0]
set best -1; set bestd 1e9
foreach i [$sW get index] xyz [$sW get {x y z}] {
    set d [veclength [vecsub $xyz $pP]]
    if {$d < $bestd} {set bestd $d; set best $i}
}
set sOw [atomselect $MID "index $best"]
puts [format "nucleophile O: index %d (P-Ow in the last frame %.2f A)" $best $bestd]

# ---- representations -----------------------------------------------------------------------
# surroundings (protein, water): thin dark-grey sticks, only residues within 5 A of the reaction core
addrep "beta < 0.5 and not hydrogen and not nucleic and same residue as (within 5 of (beta > 1.5 and beta < 2.5))" \
       {DynamicBonds 1.95 0.05 12} {ColorID 2}

# the DNA strand that is cut, in two colours
set segO3 [lindex [$sO3 get segname] 0]; set ridO3 [lindex [$sO3 get resid] 0]
set segP  [lindex [$sP  get segname] 0]; set ridP  [lindex [$sP  get resid] 0]
puts "cut between residue $ridO3 (segname $segO3, O3' side) and residue $ridP (segname $segP, P side)"
if {$segO3 eq $segP && $ridO3 != $ridP} {
    if {$ridO3 < $ridP} {set opA "<="; set opB ">="} else {set opA ">="; set opB "<="}
    addrep "segname $segO3 and resid $opA $ridO3 and not hydrogen and beta < 0.5" {DynamicBonds 1.95 0.20 16} {ColorID 0}
    addrep "segname $segP and resid $opB $ridP and not hydrogen and beta < 0.5"   {DynamicBonds 1.95 0.20 16} {ColorID 3}
} else {
    puts "WARNING: could not split the DNA strand by residue number; the strand is drawn in one colour"
    addrep "segname $segP and not hydrogen and beta < 0.5" {DynamicBonds 1.95 0.20 16} {ColorID 3}
}

# QM atoms: bonds are re-computed every frame
addrep "beta > 0.5 and beta < 2.5 and not hydrogen"  {DynamicBonds 1.95 0.15 20} {Name}
addrep "beta > 0.5 and beta < 2.5"                   {DynamicBonds 1.25 0.10 20} {Name}
addrep "beta > 0.5 and beta < 1.5 and not hydrogen"  {VDW 0.18 20}               {Name}
addrep "beta > 0.5 and beta < 2.5 and hydrogen"      {VDW 0.16 20}               {Name}
addrep "beta > 1.5 and beta < 2.5 and not hydrogen"  {VDW 0.28 24}               {Name}
addrep "beta > 2.5"                                  {VDW 0.32 24}               {ColorID 11}
color Name P orange

# ---- marks that follow every frame ---------------------------------------------------------
proc dashed {a b radius dash gap} {
    global MID
    set v [vecsub $b $a]; set L [veclength $v]
    if {$L < 1e-6} return
    set u [vecscale [expr {1.0 / $L}] $v]
    set t 0.0
    while {$t < $L} {
        set e [expr {min($t + $dash, $L)}]
        graphics $MID cylinder [vecadd $a [vecscale $t $u]] [vecadd $a [vecscale $e $u]] radius $radius resolution 16 filled yes
        set t [expr {$t + $dash + $gap}]
    }
}

set sCore0 [atomselect $MID "beta > 1.5 and beta < 2.5" frame 0]
set CORE_C [measure center $sCore0]

proc update_marks {args} {
    global SHOW_MARKS
    sync_view
    global MID sP sO3 sOw BONDCUT DASHES CORE_C LABEL_UP LABEL_RIGHT LABEL_STEP LABEL_SIZE
    set f [molinfo $MID get frame]
    foreach s [list $sP $sO3 $sOw] {$s frame $f}
    set p  [lindex [$sP  get {x y z}] 0]
    set o3 [lindex [$sO3 get {x y z}] 0]
    set ow [lindex [$sOw get {x y z}] 0]
    set d1 [veclength [vecsub $p $o3]]
    set d2 [veclength [vecsub $p $ow]]
    graphics $MID delete all
    if {!$SHOW_MARKS} return

    # breaking bond P-O3'
    if {$d1 < $BONDCUT} {
        graphics $MID color yellow
        graphics $MID cylinder $p $o3 radius 0.22 resolution 24 filled yes
        set c1 yellow; set t1 "bonded"
    } else {
        graphics $MID color red
        if {$DASHES} {dashed $p $o3 0.10 0.22 0.16}
        set c1 red; set t1 "BROKEN"
    }
    # forming bond P-Ow
    graphics $MID color green
    if {$d2 < $BONDCUT} {
        graphics $MID cylinder $p $ow radius 0.22 resolution 24 filled yes
        set t2 "bonded"
    } else {
        if {$DASHES} {dashed $p $ow 0.07 0.22 0.16}
        set t2 "forming"
    }

    # read-out at the upper left of the view (screen "right" and "up" expressed in molecule coordinates)
    set R [lindex [molinfo $MID get rotate_matrix] 0]
    set right [lrange [lindex $R 0] 0 2]
    set up    [lrange [lindex $R 1] 0 2]
    set base [vecadd $CORE_C [vecadd [vecscale $LABEL_UP $up] [vecscale $LABEL_RIGHT $right]]]
    graphics $MID color $c1
    graphics $MID text $base [format "P-O3'  %.2f A   %s" $d1 $t1] size $LABEL_SIZE thickness 2
    graphics $MID color green
    graphics $MID text [vecsub $base [vecscale $LABEL_STEP $up]] [format "P-Ow   %.2f A   %s" $d2 $t2] size $LABEL_SIZE thickness 2
}


# ---- whole complex (optional) ---------------------------------------------------------------
set HAVE_OV [file exists $OVPDB]
set OV -1
if {$HAVE_OV} {
    mol new $OVPDB type pdb waitfor all
    set OV [molinfo top]
    mol delrep 0 $OV
    mol representation NewCartoon 0.30 12 4.5
    mol selection "protein"
    mol color ColorID 6
    mol material Opaque
    mol addrep $OV
    mol representation NewCartoon 0.30 12 4.5
    mol selection "nucleic and not segname $segP"
    mol color ColorID 10
    mol addrep $OV
    mol representation NewCartoon 0.30 12 4.5
    mol selection "nucleic and segname $segP"
    mol color ColorID 3
    mol addrep $OV
    mol representation VDW 1.0 20
    mol selection "beta > 2.5 or name MN"
    mol color ColorID 11
    mol addrep $OV
    set sOVall [atomselect $OV "protein or nucleic"]
    set OV_C [measure center $sOVall]
    mol top $MID
    puts "whole complex loaded from $OVPDB ([$sOVall num] atoms)"
} else {
    puts "NOTE: $OVPDB not found - no opening shot of the whole complex"
}

proc sync_view {} {
    global MID OV HAVE_OV
    if {!$HAVE_OV} return
    foreach k {rotate_matrix center_matrix scale_matrix global_matrix} {
        molinfo $OV set $k [molinfo $MID get $k]
    }
}
proc set_cam {scale centre} {
    global MID
    molinfo $MID set scale_matrix  [list [list [list $scale 0 0 0] [list 0 $scale 0 0] [list 0 0 $scale 0] [list 0 0 0 1]]]
    molinfo $MID set center_matrix [list [transoffset [vecinvert $centre]]]
    sync_view
}
# show the opening shot / the close-up on screen (to check them before rendering)
proc view_far {} {
    global MID HAVE_OV OV OV_C ZOOM_OUT NEAR_SCALE SHOW_MARKS
    if {!$HAVE_OV} {puts "no whole-complex file"; return}
    mol on $OV; set SHOW_MARKS 0
    set_cam [expr {$NEAR_SCALE / $ZOOM_OUT}] $OV_C; update_marks
}
proc view_near {} {
    global MID HAVE_OV OV CORE_C NEAR_SCALE SHOW_MARKS HIDE_OV_IN_CLOSEUP
    set SHOW_MARKS 1
    if {$HAVE_OV && $HIDE_OV_IN_CLOSEUP} {mol off $OV}
    set_cam $NEAR_SCALE $CORE_C; update_marks
}

catch {trace remove variable vmd_frame($MID) write update_marks}
trace add variable vmd_frame($MID) write update_marks

# centre the view on the reaction core and zoom in
display resetview
molinfo $MID set center_matrix [list [transoffset [vecinvert $CORE_C]]]
scale by $ZOOM
animate goto 0
set NEAR_SCALE [lindex [molinfo $MID get scale_matrix] 0 0 0]
if {$HAVE_OV && $HIDE_OV_IN_CLOSEUP} {mol off $OV}
update_marks

proc save_view {{fn view.tcl}} {
    global MID
    set f [open $fn w]
    foreach k {rotate_matrix center_matrix scale_matrix global_matrix} {
        puts $f "molinfo $MID set $k [list [molinfo $MID get $k]]"
    }
    close $f
    puts "view written to $fn (load it later with: source $fn ; update_marks)"
}

proc make_movie {{renderer snapshot}} {
    global MID OV HAVE_OV OUTDIR ZOOM_OUT N_HOLD N_ZOOM CORE_C OV_C NEAR_SCALE SHOW_MARKS HIDE_OV_IN_CLOSEUP
    file mkdir $OUTDIR
    # the close-up is whatever zoom is on screen now (rotation is kept as it is)
    set NEAR_SCALE [lindex [molinfo $MID get scale_matrix] 0 0 0]
    set k 0
    if {$HAVE_OV} {
        mol on $OV
        set SHOW_MARKS 0
        animate goto 0
        set s0 [expr {$NEAR_SCALE / $ZOOM_OUT}]
        for {set i 0} {$i < $N_HOLD + $N_ZOOM} {incr i} {
            if {$i < $N_HOLD} {set t 0.0} else {set t [expr {double($i - $N_HOLD) / ($N_ZOOM - 1)}]}
            set t [expr {$t * $t * (3.0 - 2.0 * $t)}]
            set sc [expr {$s0 * pow($NEAR_SCALE / $s0, $t)}]
            set c [vecadd [vecscale [expr {1.0 - $t}] $OV_C] [vecscale $t $CORE_C]]
            set_cam $sc $c
            update_marks
            display update ui
            render $renderer [format "%s/f%04d.tga" $OUTDIR $k]
            incr k
        }
        if {$HIDE_OV_IN_CLOSEUP} {mol off $OV}
    }
    set SHOW_MARKS 1
    set_cam $NEAR_SCALE $CORE_C
    set n [molinfo $MID get numframes]
    for {set i 0} {$i < $n} {incr i} {
        animate goto $i
        update_marks
        display update ui
        render $renderer [format "%s/f%04d.tga" $OUTDIR $k]
        incr k
    }
    puts "$k frames written to $OUTDIR/  ([format %.1f [expr {$k / 24.0}]] s at 24 frames per second)"
}

puts "loaded $NF frames."
puts "yellow stick = P-O3' bond (red dashes once broken) ; green = new P-Ow bond ; blue / orange = the two sides of the cut DNA strand"
puts "Adjust the view, type  update_marks  to refresh the text, check with the frame slider, then run:  make_movie"
puts "With the whole complex: type  view_far  /  view_near  to check the opening shot and the close-up."
