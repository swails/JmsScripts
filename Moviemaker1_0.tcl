#A procedure to make a toplevel window
set xcom {frame .top.toprotate.frmx -relief groove -pady 5}
proc makeTop { } {
	global x y z t scalefac framea frameb rot zoom fade clip traj trans i fadearr fadei degx degy degz trajec scale time degtx degty newmat degtz
	set i 1
	set newmat 22
	set fadei 0
	#creates master gui
	toplevel .top
	#Put things in it
	#Make title, should be packed first
	label .top.lab -text "Advanced Movie Maker 1.0" -font "ansi 12 bold"
label .top.labe -text "This enables the user to create movies that utilize various actions 
that one can currently perform in VMD's OpenGL display(such as zooming and rotating). Once
'Start!' is pressed, these actions are defined in chronological order in this program. Once 
'done' is pressed, the movie will render the frames needed to make the movie. Once these 
frames are created, a frame merger, such as videomach, must be used to create the movie." -font "ansi 10 bold" -padx 0
#executed when render button is pressed
	proc push_button {} {
		global rot zoom trans fade clip traj i degx degy degz trajec scale degtx degty degtz newmat
		for {set j 1} {$j <= 25} {incr j 1} {
			set degx($j) 0
			set degy($j) 0
			set degz($j) 0
			set degtx($j) 0
			set degty($j) 0
			set degtz($j) 0
			set scale($j) 1
			set time($j) 0
			set framea($j) 0
			set frameb($j) 0
		}
		toplevel .top.topint
		label .top.topint.label -text "What do you want to do in this action?" -font "ansi 10 bold"
		button .top.topint.but -text "Continue" -command "push_button1"
		checkbutton .top.topint.rot -text "Rotate" -variable rot
		checkbutton .top.topint.zoom -text "Zoom" -variable zoom
		checkbutton .top.topint.trans -text "Translate" -variable trans
		checkbutton .top.topint.fade -text "Fade" -variable fade		
		checkbutton .top.topint.clip -text "Clip" -variable clip($i)
		checkbutton .top.topint.traj -text "Advance Trajectory" -variable traj($i)
		grid .top.topint.label -in .top.topint -row 1 -column 1 -columnspan 2
		grid .top.topint.rot -in .top.topint -row 2 -column 1
		grid .top.topint.zoom -in .top.topint -row 2 -column 2
		grid .top.topint.trans -in .top.topint -row 2 -column 3
		grid .top.topint.fade -in .top.topint -row 3 -column 1
		grid .top.topint.clip -in .top.topint -row 3 -column 2
		grid .top.topint.traj -in .top.topint -row 3 -column 3
		grid .top.topint.but -in .top.topint -row 4 -column 2
	}		
	proc push_button2 {} {

		destroy .top.toprotate
		global rot zoom trans fade clip traj i fadearr fadei fadedum
		incr fadei 1
		set fadedum 0
		set i [expr $i+1]
		toplevel .top.topint
		label .top.topint.label -text "What do you want to do in this action?" -font "ansi 10 bold"
		button .top.topint.but -text "Continue" -command "push_button1"
		checkbutton .top.topint.rot -text "Rotate" -variable rot
		checkbutton .top.topint.zoom -text "Zoom" -variable zoom
		checkbutton .top.topint.trans -text "Translate" -variable trans
		checkbutton .top.topint.fade -text "Fade" -variable fade		
		checkbutton .top.topint.clip -text "Clip" -variable clip($i)
		checkbutton .top.topint.traj -text "Advance Trajectory" -variable traj($i)
		grid .top.topint.label -in .top.topint -row 1 -column 1 -columnspan 2
		grid .top.topint.rot -in .top.topint -row 2 -column 1
		grid .top.topint.zoom -in .top.topint -row 2 -column 2
		grid .top.topint.trans -in .top.topint -row 2 -column 3
		grid .top.topint.fade -in .top.topint -row 3 -column 1
		grid .top.topint.clip -in .top.topint -row 3 -column 2
		grid .top.topint.traj -in .top.topint -row 3 -column 3
		grid .top.topint.but -in .top.topint -row 4 -column 2
	}	
	proc push_button1 {} {
		global x y z t scalefac cliprep start end framea frameb i rot zoom trans fade fadearr fademat fadedum fadei newmat clip traj trans time degx degy degz degtx degty degtz newmat fadearr
		destroy .top.topint
		set fadedum 1
		toplevel .top.toprotate
	#The following block containts all specifications for the action blocks.
	#The variables that will be utilized from this code are the frames created
		#Rotating Frames
		frame .top.toprotate.frmrot -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmx -relief groove -pady 5
		label .top.toprotate.labrotatex -text "Degrees around x-axis.(Horizontal axis)" -width 50
		entry .top.toprotate.entrotatex -textvariable degx($i)
		frame .top.toprotate.frmy -relief groove -pady 5
		label .top.toprotate.labrotatey -text "Degrees around y-axis.(Vertical axis)" -width 50
		entry .top.toprotate.entrotatey -textvariable degy($i)
		frame .top.toprotate.frmz -relief groove -pady 5
		label .top.toprotate.labrotatez -text "Degrees around z-axis.(Axis points out of screen)" -width 50
		entry .top.toprotate.entrotatez -textvariable degz($i)
		#Time Frame
		frame .top.toprotate.frmtime -relief groove -pady 5 -borderwidth 1
		#Translation Frames
		frame .top.toprotate.frmtrans -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmtx -relief groove -pady 5
		label .top.toprotate.labrotatetx -text "Units along x-axis.(Horizontal motion)" -width 50
		entry .top.toprotate.entrotatetx -textvariable degtx($i)
		frame .top.toprotate.frmty -relief groove -pady 5
		label .top.toprotate.labrotatety -text "Units along y-axis.(Vertical motion)" -width 50
		entry .top.toprotate.entrotatety -textvariable degty($i)
		frame .top.toprotate.frmtz -relief groove -pady 5
		label .top.toprotate.labrotatetz -text "Units along z-axis.(Forward and backward motion)" -width 50
		entry .top.toprotate.entrotatetz -textvariable degtz($i)
		label .top.toprotate.labtimerot -text "Duration in seconds" -width 50
		entry .top.toprotate.entrotatetime -textvariable time($i)
		frame .top.toprotate.frmzoom -relief groove -pady 5 -borderwidth 1
		label .top.toprotate.labzoom -text "Scaling factor for zoom" -width 50
		entry .top.toprotate.entzoom -textvariable scale($i)
		frame .top.toprotate.frmtraj -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.subtraja -relief groove
		frame .top.toprotate.subtrajb -relief groove
		label .top.toprotate.labtraj -text "Advance Trajectory" -width 25
		label .top.toprotate.labsuba -text "From frame number" -width 25
		label .top.toprotate.labsubb -text "To frame number" -width 25
		entry .top.toprotate.entsuba -textvariable framea($i)
		entry .top.toprotate.entsubb -textvariable frameb($i)
		#Fade Frames
		frame .top.toprotate.frmfade -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmfadea -relief groove -pady 5
		frame .top.toprotate.frmfadeb -relief groove -pady 5
		label .top.toprotate.labfadea -text "RepId to Fade(First rep is 0)" -width 50
		label .top.toprotate.labfadeb -text "Material of rep to fade" -width 50
		entry .top.toprotate.entfadea -textvariable fadeselect($i)
		entry .top.toprotate.entfadeb -textvariable fademat($i)
		#clip frames
		frame .top.toprotate.frmclip -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmclipa -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmclipb -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmclipba -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.frmclipbb -relief groove -pady 5 -borderwidth 1
		label .top.toprotate.labclipa -text "RepId to Clip(First rep is 0)"
		label .top.toprotate.labclipba -text "Start clip at" -width 30
		label .top.toprotate.labclipbb -text "End clip at" -width 30
		entry .top.toprotate.entclipa -textvariable cliprep($i) -width 20
		scale .top.toprotate.scalclipba -orient h -width 30 -digit 1 -from -100 -to 100 -variable start($i) -tickinterval 20
		scale .top.toprotate.scalclipbb -orient h -width 30 -digit 1 -from -100 -to 100 -variable end($i) -tickinterval 20
		#And now pack frames!
		pack .top.toprotate.labclipa -in .top.toprotate.frmclipa
		pack .top.toprotate.entclipa -in .top.toprotate.frmclipa
		pack .top.toprotate.labclipba -in .top.toprotate.frmclipba
		pack .top.toprotate.scalclipba -in .top.toprotate.frmclipba -fill x
		pack .top.toprotate.labclipbb -in .top.toprotate.frmclipbb
		pack .top.toprotate.scalclipbb -in .top.toprotate.frmclipbb -fill x
		pack .top.toprotate.frmclipba -in .top.toprotate.frmclipb
		pack .top.toprotate.frmclipbb -in .top.toprotate.frmclipb
		grid .top.toprotate.frmclipa -in .top.toprotate.frmclip -row 1 -column 1
		grid .top.toprotate.frmclipb -in .top.toprotate.frmclip -row 1 -column 2 -columnspan 2
		pack .top.toprotate.labsuba -in .top.toprotate.subtraja
		pack .top.toprotate.entsuba -in .top.toprotate.subtraja
		pack .top.toprotate.labsubb -in .top.toprotate.subtrajb
		pack .top.toprotate.entsubb -in .top.toprotate.subtrajb
		grid .top.toprotate.labtraj -in .top.toprotate.frmtraj -row 1 -column 1
		grid .top.toprotate.subtraja -in .top.toprotate.frmtraj -row 2 -column 1
		grid .top.toprotate.subtrajb -in .top.toprotate.frmtraj -row 2 -column 2
		pack .top.toprotate.labrotatex -in .top.toprotate.frmx
		pack .top.toprotate.entrotatex -in .top.toprotate.frmx
		pack .top.toprotate.labrotatey -in .top.toprotate.frmy
		pack .top.toprotate.entrotatey -in .top.toprotate.frmy
		pack .top.toprotate.labrotatez -in .top.toprotate.frmz
		pack .top.toprotate.entrotatez -in .top.toprotate.frmz
		pack .top.toprotate.labrotatetx -in .top.toprotate.frmtx
		pack .top.toprotate.entrotatetx -in .top.toprotate.frmtx
		pack .top.toprotate.labrotatety -in .top.toprotate.frmty
		pack .top.toprotate.entrotatety -in .top.toprotate.frmty
		pack .top.toprotate.labrotatetz -in .top.toprotate.frmtz
		pack .top.toprotate.entrotatetz -in .top.toprotate.frmtz
		pack .top.toprotate.labtimerot -in .top.toprotate.frmtime
		pack .top.toprotate.entrotatetime -in .top.toprotate.frmtime
		pack .top.toprotate.labzoom -in .top.toprotate.frmzoom
		pack .top.toprotate.entzoom -in .top.toprotate.frmzoom
		pack .top.toprotate.labfadea -in .top.toprotate.frmfadea
		pack .top.toprotate.entfadea -in .top.toprotate.frmfadea
		pack .top.toprotate.labfadeb -in .top.toprotate.frmfadeb
		pack .top.toprotate.entfadeb -in .top.toprotate.frmfadeb
		pack .top.toprotate.frmx -in .top.toprotate.frmrot
		pack .top.toprotate.frmy -in .top.toprotate.frmrot
		pack .top.toprotate.frmz -in .top.toprotate.frmrot
		pack .top.toprotate.frmtx -in .top.toprotate.frmtrans
		pack .top.toprotate.frmty -in .top.toprotate.frmtrans
		pack .top.toprotate.frmtz -in .top.toprotate.frmtrans
		pack .top.toprotate.frmfadea -in .top.toprotate.frmfade
		pack .top.toprotate.frmfadeb -in .top.toprotate.frmfade
		#The frames above are prepared and ready. 
		#Available frames are frmrot frmtx frmty frmtz frmtime frmzoom frmtraj frmfadea frmfadeb

		#The following are the buttons that will always appear
		frame .top.toprotate.allfrm -relief groove -pady 5 -borderwidth 1 -width 50
		frame .top.toprotate.prevfram -relief groove -pady 5 -borderwidth 1
		frame .top.toprotate.contfram -relief groove -pady 5 -borderwidth 1
		button .top.toprotate.button -text "Preview" -command "previes" -width 25
		button .top.toprotate.button2 -text "Undo Preview" -command "deprevies" -width 25
		button .top.toprotate.buttonnext -text "Next Action" -command "push_button2" -width 25
		button .top.toprotate.button3 -text "Done" -command "confirm" -width 25
		button .top.toprotate.button4 -text "Cancel" -command { destroy .top.toprotate } -width 40
		pack .top.toprotate.button -in .top.toprotate.prevfram
		pack .top.toprotate.button2 -in .top.toprotate.prevfram
		pack .top.toprotate.buttonnext -in .top.toprotate.contfram
		pack .top.toprotate.button3 -in .top.toprotate.contfram
		grid .top.toprotate.prevfram -in .top.toprotate.allfrm -row 1 -column 1
		grid .top.toprotate.contfram -in .top.toprotate.allfrm -row 1 -column 2
		proc confirm {} {
		global name pov snap
		toplevel .top.confirm
		frame .top.confirm.filename -relief groove -pady 5 -borderwidth 1 -width 50
		frame .top.confirm.filetype -relief groove -pady 5 -borderwidth 1 -width 50
		frame .top.confirm.exec -relief groove -pady 5 -borderwidth 1 -width 50
		label .top.confirm.filelabel -text "What filename base?" -width 50 -font "ansi 10 bold"
		entry .top.confirm.ent -textvariable name
		label .top.confirm.checklabel -text "Which renderer?" -width 50 -font "ansi 10 bold"
		checkbutton .top.confirm.pov -text "Pov-ray" -variable pov
		checkbutton .top.confirm.snap -text "Snapshot" -variable snap
		button .top.confirm.done -text "Create Movie" -command "setvariables2" -width 20
		button .top.confirm.cancel -text "Cancel" -command { destroy .top.confirm } -width 20
		pack .top.confirm.filelabel -in .top.confirm.filename
		pack .top.confirm.ent -in .top.confirm.filename
		grid .top.confirm.checklabel -in .top.confirm.filetype -row 1 -column 1 -columnspan 2
		grid .top.confirm.pov -in .top.confirm.filetype -row 2 -column 1
		grid .top.confirm.snap -in .top.confirm.filetype -row 2 -column 2
		grid .top.confirm.done -in .top.confirm.exec -row 1 -column 1
		grid .top.confirm.cancel -in .top.confirm.exec -row 1 -column 2
		pack .top.confirm.filename .top.confirm.filetype .top.confirm.exec
		}
		proc setvariables2 {} {
		global i frameb framea degx degy degz scale cliprep start end time traj degtx degty degtz clip name pov snap fadeselect fadedum fademat fade newmat fadearr
		set kk 1
		set fadedum 1
		for {set k 1} {$k <= $i} {incr k 1} {
			if { $traj($k) == 1 } {
				set frames [expr $frameb($k)-$framea($k)]
				set prevfram $framea($k)
			} else {
			set frames [expr $time($k)*30]
			}
			if { $fadearr([expr $fadedum-1]) == 1} {
				set hconst 0
				foreach {f} [material settings $fademat($fadedum)] {
 				set dumarray($hconst) $f
				incr hconst 1
				}
				set matfac [expr $dumarray(4)/$frames]
				set initopac $dumarray(4)
				material add copy $fademat($fadedum)
				#first new material is named Material20
				mol modmaterial $fadeselect($fadedum) top Material$newmat
			} else {
				set matfac 0
			}
			set normal_vector [vecnorm [vectrans [measure inverse [lindex [molinfo top get rotate_matrix] 0]] {0 0 -1}]]
			set clip_center [vecadd "0.0 0.0 0.0" [vecscale 0.0 $normal_vector]]
			set rotfacx [expr $degx($k)/$frames]
			set rotfacy [expr $degy($k)/$frames] 
			set rotfacz [expr $degz($k)/$frames]
			set transfacx [expr $degtx($k)/$frames]
			set transfacy [expr $degty($k)/$frames] 
			set transfacz [expr $degtz($k)/$frames]
			set scalfac [expr pow($scale($k),1.000/$frames)]
			if { $clip($fadedum) == 1 } {
				mol clipplane status 0 $cliprep($fadedum) top 2
				mol clipplane color 0 $cliprep($fadedum) top "0.25 0.75 0.75"
				mol clipplane normal 0 $cliprep($fadedum) top $normal_vector
				set dist $start($fadedum)
				set distfac [expr ($end($fadedum)-$start($fadedum))/$frames] 
				}
			for {set l 1} {$l <= $frames} {incr l 1} {
				if { $fadearr([expr $fadedum-1]) == 1 } {
				material change opacity Material$newmat $initopac
				set initopac [expr $initopac-$matfac]
				}
				scale by $scalfac
				rotate x by $rotfacx
				rotate y by $rotfacy
				rotate z by $rotfacz
				translate by $transfacx $transfacy $transfacz
				if { $traj($k) == 1 } {
					animate goto $prevfram
					set prevfram [expr $prevfram+1]
				}
				set filename $name.[format "%04d" $kk].bmp
				if { $pov == 1} {
					set filename $name.[format "%04d" $kk].pov
					render POV3 $filename
				} elseif { $snap == 1 } {
					set filename $name.[format "%04d" $kk].bmp
					render snapshot $filename
				} else {
					puts "Error something has gone horribly wrong!"
				}
				if { $clip($k) == 1 } {
					mol clipplane center 0 $cliprep($k) top "0 0 $dist"
					set dist [expr $dist+$distfac]
				}
				incr kk 1
				display update
			}
			if { $fadearr([expr $fadedum-1]) == 1} {
				 incr newmat 1
			}
			incr fadedum 1
			}
		}
		if { $rot == 1 } {
			pack .top.toprotate.frmrot
		}
		if { $trans == 1 } {
			pack .top.toprotate.frmtrans
		}
		if { $zoom == 1 } {
			pack .top.toprotate.frmzoom
		}
		set fadearr($fadei) $fade
                if { $fade == 1 } {
			pack .top.toprotate.frmfade
		}
		if { $clip($i) == 1 } {
			pack .top.toprotate.frmclip
		}
		if { $traj($i) == 1 } {
			pack .top.toprotate.frmtraj
		} else { pack .top.toprotate.frmtime 
		}
	pack .top.toprotate.allfrm .top.toprotate.button4
	}
	
	label .top.labrotate -text "Press button to begin"
	button .top.butrotate -text "start!" -command "push_button"
	#An option to close the window.
	button .top.butclose -text "Close" -command { destroy .top }
	#Pack everything
	pack .top.lab .top.labe .top.labrotate .top.butrotate .top.butclose
	proc previes {} {
		global i frameb framea degx cliprep start end clip degy degz scale time traj fadedum degtx degty degtz fadeselect fademat fade newmat fadearr
		set kk 1
		set fadedum 1
		for {set k 1} {$k <= $i} {incr k 1} {
			if { $traj($k) == 1 } {
				set frames [expr $frameb($k)-$framea($k)]
				set prevfram $framea($k)
			} else {
			set frames [expr $time($k)*30.0]
			}
			if { $fadearr([expr $fadedum-1]) == 1} {
				set hconst 0
				foreach {f} [material settings $fademat($fadedum)] {
 				set dumarray($hconst) $f
				incr hconst 1
				}
				set matfac [expr ($dumarray(4))/($frames)]
				set initopac $dumarray(4)
				material add copy $fademat($fadedum)
				#first new material is named Material20
				mol modmaterial $fadeselect($fadedum) top Material$newmat
			} else {
				set matfac 0
			}
			set normal_vector [vecnorm [vectrans [measure inverse [lindex [molinfo top get rotate_matrix] 0]] {0 0 -1}]]
			set clip_center [vecadd "0.0 0.0 0.0" [vecscale 0.0 $normal_vector]]
			set rotfacx [expr $degx($k)/$frames]
			set rotfacy [expr $degy($k)/$frames] 
			set rotfacz [expr $degz($k)/$frames]
			set transfacx [expr $degtx($k)/$frames]
			set transfacy [expr $degty($k)/$frames] 
			set transfacz [expr $degtz($k)/$frames]
			if { $clip($fadedum) == 1 } {
				mol clipplane status 0 $cliprep($fadedum) top 2
				mol clipplane color 0 $cliprep($fadedum) top "0.25 0.75 0.75"
				mol clipplane normal 0 $cliprep($fadedum) top $normal_vector
				set dist $start($fadedum)
				set distfac [expr ($end($fadedum)-$start($fadedum))/$frames] 
				}
			set scalfac [expr pow($scale($k),1.000/$frames)]
			for {set l 1} {$l <= $frames} {incr l 1} {
				if { $fadearr([expr $fadedum-1]) == 1 } {
				material change opacity Material$newmat $initopac
				set initopac [expr $initopac-$matfac]
				}
				scale by $scalfac
				rotate x by $rotfacx
				rotate y by $rotfacy
				rotate z by $rotfacz
				translate by $transfacx $transfacy $transfacz
				if { $traj($k) == 1 } {
					animate goto $prevfram
					set prevfram [expr $prevfram+1]
				}
				if { $clip($k) == 1 } {
					mol clipplane center 0 $cliprep($k) top "0 0 $dist"
					set dist [expr $dist+$distfac]
				}
				display update
			}
			if { $fadearr([expr $fadedum-1]) == 1} {
				 incr newmat 1
			}
			incr fadedum 1
			}
		}
	proc deprevies {} {
		global i frameb framea degx degy degz scale start end cliprep time traj degtx degty degtz fadeselect fademat clip fade newmat fadearr fadedum

		for {set k $i} {$k >= 1} {incr k -1} {
			if { $traj($k) == 1 } {
				set frames [expr $frameb($k)-$framea($k)]
				set prevfram $framea($k)
			} else {
			set frames [expr $time($k)*30.0]
			}
		incr fadedum -1
		incr newmat -1
		set rot_factorx2 [expr -1.0*($degx($k)/($frames))]
		set rot_factory2 [expr -1.0*($degy($k)/($frames))]
		set rot_factorz2 [expr -1.0*($degz($k)/($frames))]
		set transfacx2 [expr -1.0*($degtx($k)/$frames)]
		set transfacy2 [expr -1.0*($degty($k)/$frames)] 
		set transfacz2 [expr -1.0*($degtz($k)/$frames)]
		set scalfac2 [expr 1.0/(pow($scale($k),1.0/$frames))]
		if { $clip($k) == 1 } {
			mol clipplane status 0 $cliprep($k) top 0
			}
		for {set l 0} {$l < $frames} {incr l 1} {  
			if { $traj($k) == 1 } {
				animate goto $prevfram
				set prevfram [expr $prevfram-1]
			}
			translate by $transfacx2 $transfacy2 $transfacz2
			rotate z by $rot_factorz2 
			rotate y by $rot_factory2 
			rotate x by $rot_factorx2 
			scale by $scalfac2
		}
		if { $fadearr([expr $fadedum-1]) == 1 } {
			mol modmaterial $fadeselect($fadedum) top $fademat($fadedum)
			material delete Material$newmat
		}
		if { $fadearr([expr $fadedum-1]) == 1} {
		}
		}
	}
}

