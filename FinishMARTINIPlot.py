#!/usr/bin/python

import os,AutomatedPlot,HAConf, sys

gromacs = HAConf.configuration['gromacs_location']
vmd = HAConf.programs['vmd']
from HAConf import configuration
wd = os.getcwd()

try:
	initial_sequence = sys.argv[1]
except:
	initial_sequence = wd.split('/')[-2]

length = len(initial_sequence)

if length %2 == 0:
	middle_residues = str(length/2 - 1) + "-" + str(length/2 +2)
else:
	middle_residues = str(length/2 - 1) + "-" + str(length/2 +1)

make_ndx_command = 'echo "a b*\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n 16 & 17\n 16 & 18\n 16 & 19\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'

g_bundle_command = "echo '20\n21\n22\n'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "


os.system(make_ndx_command)
#print make_ndx_command
os.system(g_bundle_command)

AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",initial_sequence)
AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",initial_sequence)

AutomatedPlot.tilt_data_plot("bun_kink.xvg","kink_V_time.png",initial_sequence)
AutomatedPlot.hist_tilt_data_plot("bun_kink.xvg","hist_kink.png",initial_sequence)

first_trjconv_command = "echo '1\n0\n' | " + gromacs + "trjconv -f t_0.xtc -s t_0.tpr -center -pbc mol -o center.xtc"

trjconv_xy_command = "echo '1\n0\n' | " + gromacs + "trjconv_xy -f center.xtc -s t_0.tpr  -o xy_fit.xtc -fit rotxy+transxy"
os.system(first_trjconv_command)
os.system(trjconv_xy_command)
'''
vmd_command = vmd + " xy_fit.xtc t_0.gro -dispdev text -eofexit <" + configuration['vmd_state_files'] + "automated_analysis.tcl > vmd.log"
os.system(vmd_command)

AutomatedPlot.hist_com_difference_plot("bilayer_position.dat","bilayer_position.png",initial_sequence)
'''
vmd_command = vmd + " xy_fit.xtc t_0.gro -dispdev text -eofexit <" + configuration['vmd_state_files'] + "automated_analysis.tcl > vmd.log"
os.system(vmd_command)

AutomatedPlot.hist_com_difference_plot("bilayer_position.dat","bilayer_position.png",initial_sequence)

AutomatedPlot.bilayer_deformation_refined("up.dat","down.dat","both-bully.png",graph_title=initial_sequence)

AutomatedPlot.rotation_plot("rot.dat","rotation.png",initial_sequence)
