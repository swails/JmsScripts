proc vmd_labelcb_xmgr { args } {

#source ~/devel/git/Projects/giambasu/VMDscripts/smooth_angles.tcl

global vmd_graph_label

set llength_vmd_graph_label [llength $vmd_graph_label]

if {$llength_vmd_graph_label==2} {
         
         puts "CORRELATION \n"
		set type_x [lindex $vmd_graph_label 0 0]
		set id_x   [lindex $vmd_graph_label 0 1]
		set type_y [lindex $vmd_graph_label 1 0]
		set id_y   [lindex $vmd_graph_label 1 1]
		set data_x [label graph $type_x $id_x]
		set info_x [label list $type_x]
		set data_y [label graph $type_y $id_y]
		set info_y [label list $type_y]
		set input "@type xy\n@  title \"Correlation Chart\" \n"

		switch [lindex $type_x] {
				"Bonds" {
					 set fin 2
					}
					  "Angles" {
					 set fin 3
						# set data_x [smooth_angles $data_x]
					}
				"Dihedrals" {
					 set fin 4
					 # set data_x [smooth_angles $data_x]
					}
					 }
	
		set label_x " "
	
		for {set deb 0} { $deb < $fin } {incr deb 1} {
				set atom [lindex [lindex [lindex $info_x $id_x] $deb] 1]
				set mol [lindex [lindex [lindex $info_x $id_x] $deb] 0]
				set sel [atomselect $mol "index $atom"]
				set atomname [$sel get name]
				set resid [$sel get resid]
				set resname [$sel get resname]
				append label_x "$resname $resid:$atomname - "
							 }

		set label_x [string trimright $label_x "- "]
														

		switch [lindex $type_y] {
				"Bonds" {
					set fin 2
					set units_y "\"Distance \(\\f\{Times-Roman\}\\cE\\C\)\""
					}
						 "Angles" {
					set fin 3
					set units_x "\"Angle  \(Deg.\) \""
				#	set data_y [smooth_angles $data_y]
					}
				"Dihedrals" {
					set fin 4
					set units_y "\"Angle  \(Deg.\) \""
				#	set data_y [smooth_angles $data_y]
					}
					}


																				
		set label_y " "

		for {set deb 0} { $deb < $fin } {incr deb 1} {
								set atom [lindex [lindex [lindex $info_y $id_y] $deb] 1]
								set mol [lindex [lindex [lindex $info_y $id_y] $deb] 0]
								set sel [atomselect $mol "index $atom"]
								set atomname [$sel get name]
								set resid [$sel get resid]
								set resname [$sel get resname]
								append label_y "$resname $resid:$atomname - "
									 }

		set label_y [string trimright $label_y "- "]


		append input "@  xaxis label \" $label_x \"  \n"
		append input "@  yaxis label \" $label_y \"  \n"
	
		foreach xdata $data_x ydata $data_y {
		append input "  $xdata  $ydata\n"
											}
	
		set xmgrout "correlation.agr"
	
			set rc [catch {exec xmgrace -maxpath 50000 -saveall $xmgrout -pipe << $input &} msg]
		if { $rc } {
				vmd_labelcb_save
			   }

} elseif { $llength_vmd_graph_label != 2} {
      puts "PLOT \n"
      #SET DEFAULTS
      set data { 1 2 3}
	  set data_hist {1 2 3}
	  set list_bins {1 2 3}
	  set data_mean 1
	  set data_stdv 1
	  set iq4r 1
	  set iq10r 1
	  set xmgrout "1.agr"
	  set xmin 0
	  set ymin 0
	  set xmax 100000
	  set ymax 360
	  set LABEL "NONE"
	  set title_xmgr "NONE"
	  #DONE DEFAULTS




	  #NOW GO OVER EACH LABEL
      foreach item $vmd_graph_label {
				
				set type [lindex $vmd_graph_label 0 0]
				set data [label graph $type [lindex $item 1]]
				set info [label list $type]
				set ind [lindex $item 1]
				set data_l [llength $data]
				
				switch [lindex $vmd_graph_label 0 0] {
						"Bonds" {
							set fin 2
							set YLABEL "Distance \(\\f\{Times-Roman\}\\cE\\C\) "
							set spacing 0.01
							set tick 1
							set xmgrout bond_
							}
						 "Angles" {
							set fin 3
							set YLABEL "Angle  \(Deg.\)"
							set spacing 1.0
							set tick 60
							set xmgrout angle_
						#	set data [smooth_angles $data]
							}
						"Dihedrals" {
							set fin 4
							set YLABEL "Torsion  \(Deg.\)"
							set spacing 1.0
							set tick 60
							set xmgrout dihed_
						#	set data [smooth_angles $data]
							}
							}
				set title_xmgr {}
				for {set deb 0} { $deb < $fin } {incr deb 1} {
	
							set atom [lindex [lindex [lindex $info $ind] $deb] 1]
							set mol [lindex [lindex [lindex $info $ind] $deb] 0]
							set sel [atomselect $mol "index $atom"]
							set atomname [$sel get name]
							set resid [$sel get resid]
							set resname [$sel get resname]
							append title_xmgr "$resname $resid:$atomname - "
							append xmgrout "$resname:$resid:$atomname:" 
											}
			
				set title_xmgr [string trimright $title_xmgr "- "]
				set xmgrout [string trimright $xmgrout ":"]
				append xmgrout ".agr"


				#####################
				#Do the calculations#
				#####################
				set data_mean [string range [mean $data] 0 5]
				set data_stdv [string range [stdev $data] 0 5]
				set data_med  [string range [median $data] 0 5]

                # compute histogram
				set hist_low  [expr round([expr [min $data] - 1])]
				set hist_high [expr round([expr [max $data] + 1])]
				set numb_bins [expr $hist_high - $hist_low]
				set numb_bins [expr $numb_bins/$spacing]
				
				set list_bins {}
				for {set i 0} {$i <= $numb_bins} {incr i} {
				set bin [expr $hist_low + $i*$spacing]
				lappend list_bins $bin
				}

				set data_hist_u [histogram $list_bins $data]
				set data_length [llength $data]
				set data_length [expr 1.0/$data_length]
				#normalize histogram
				set data_hist [vecscale $data_length $data_hist_u]
				set max_hist  [max $data_hist]
                # compute IQR's

                set iq4r   [expr [quantiles  $list_bins $data_hist_u 0.75] - [quantiles $list_bins $data_hist_u 0.25] ]
				set iq4r   [string range $iq4r 0 5]
                set iq10r  [expr [quantiles  $list_bins $data_hist_u 0.90] - [quantiles $list_bins $data_hist_u 0.10] ]
				set iq10r   [string range $iq10r 0 5]



                #PLOT THE DATA
				set    input "@g0 on \n"
				append input "@g0 type XY \n"
				append input "@with g0 \n"
				append input "@    world 0, $hist_low, $data_l, $hist_high \n"
				append input "@    view 0.150000, 0.531818, 1.143349, 0.850000 \n"
				append input "@    title \"$title_xmgr\" \n"
				append input "@    subtitle \"STATS: AV $data_mean, STD $data_stdv, MED $data_med, IQ4R $iq4r, IQ10R $iq10r \" \n"
				append input "@    xaxis  label \"Time(ps)\" \n"
				append input "@    xaxis  tick major 5000 \n"
				append input "@    yaxis label \" $YLABEL \" \n"
				append input "@    yaxis  tick major $tick \n"
				append input "@    s0 line type 1\n"
				append input "@    s0 line pattern 1\n"
				append input "@    s0 dropline off\n"
				append input "@g1 on \n"
				append input "@g1 type XY \n"
				append input "@with g1 \n"
				append input "@    world $hist_low, 0, $hist_high, $max_hist \n"
				append input "@    view 0.150000, 0.150000, 1.143349, 0.468182 \n"
				append input "@    xaxis  label \" $YLABEL \" \n"
				append input "@    xaxis  tick major $tick \n"
				append input "@    yaxis  label \"COUNT\" \n"
				append input "@    yaxis  tick major 0.01 \n"
				append input "@    s0 line type 2\n"
				append input "@    s0 line color 2\n"
				append input "@    s0 line pattern 1\n"
				append input "@    s0 dropline on\n"
				append input "@target G0.S0 \n"
				append input "@type xy \n"
				#DATA
				set i 0
				foreach elem $data {
					append input "  $i $elem\n"
					incr i
				}
				append input "& \n"
				append input "@target G1.S0 \n"
				append input "@type xy \n"
				#HIST
			
				set i 0
				foreach elem $data_hist {
				    set bin [expr $hist_low + $i*$spacing]
					append input " $bin $elem \n"
					set i [expr $i + 1.0]
				}
				append input "&"
				# make xmgrace pop up 
				set rc [catch {exec xmgrace -autoscale none -maxpath 50000 -saveall $xmgrout -pipe << $input &} msg]
				if { $rc } {
				vmd_labelcb_save
					   }
         
		 
		 
		 #foreach
		 }		 

#elseif
}

#end procedure
}
