# calculate SASA for every residue of a selection
# usage: getAllResSASA "selection" probe_radius <startframe <endframe>>

proc getAllResSASA { selection radius args } {
  puts "SASA value calculator for $selection residues, probe radius $radius"

  # get number of frames
  set numframes [molinfo top get numframes]
  puts "total number of frames in trajectory: $numframes"
  if {[llength $args] == 0} {
    set startframe 1
    set endframe $numframes
  }
  if {[llength $args] == 1} {
    set startframe [lindex $args 0]
    set endframe $numframes
    set stepframe 1
    if {$startframe <= 0 || $startframe > $numframes} {
      puts "illegal value of startframe, changing to first frame"
      set startframe 1
    }
  }
  if {[llength $args] >= 2} {
    set startframe [lindex $args 0]
    set endframe [lindex $args 1]
    if {$startframe <= 0 || $startframe > $numframes} {
      puts "illegal value of startframe, changing to first frame"
      set startframe 1
    }
    if {$endframe < $startframe || $endframe > $numframes} {
      puts "illegal value of endframe, changing to last frame"
      set endframe $numframes
    }
  }

  set totframes [expr ($endframe - $startframe + 1)]
  puts "analysis will be performed on $totframes frame(s) ($startframe to $endframe)"

  # get list of residues
  set allsel [atomselect top $selection]
  set residlist [lsort -unique -integer [ $allsel get residue ]]

  # resid map for every atom  
  set allResid [$allsel get residue]
  # resname map for every atom
  set allResname [$allsel get resname]
  # create resid->resname map
  foreach resID $allResid resNAME $allResname {
    set mapResidResname($resID) $resNAME
  }
  
  # set sasa for every residue to 0.0
  foreach r $residlist {
    lappend resSASA 0.0
  }
  
  # now subtract 1 from all frame indexes - numbering starts with 0
  set startframe [expr $startframe - 1]
  set endframe [expr $endframe - 1]

  puts "go and get coffee..."
  for {set i $startframe; set d 1} {$i <= $endframe} {incr i; incr d} {
    # show activity - one dot for every frame
    puts -nonewline "."
    if { [expr $d % 50] == 0 } {
      puts " "    
    }
    flush stdout
    
    $allsel frame $i
    $allsel update
  
    foreach r $residlist {
      set sel [atomselect top "residue $r" frame $i]
      set rsasa [measure sasa $radius $allsel -restrict $sel]
      # set user value for this frame
      $sel set user $rsasa
      $sel delete
      # puts "residue $r, sasa: $rsasa"
      # if sasa is above threshold, store it
      #if {$rsasa >= $sasa_limit} {
      lset resSASA $r [ expr {[lindex $resSASA $r] + $rsasa/$totframes} ]
      # }
    }
  }
  
  set fw [open "res_sasa.dat" w]
  foreach r $residlist {
    foreach {tmp resName} [split [array get mapResidResname $r] ] break
    puts $fw "$r $resName [lindex $resSASA $r]"
  }
  close $fw
  puts "done"

  # uncomment following code to change coloring method to "user based"
  #mol modcolor  0 [molinfo top] User
  #mol colupdate 0 [molinfo top] 1
  #mol scaleminmax [molinfo top] 0 auto
}

