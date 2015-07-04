# module input_deck_m.py

import math
import string
import sys

#-----------------------------------------------------------------------
def read_deck(filename):

	print '\n   Reading input deck ...'
	file = open(filename,'r')

	EOF = False
	while not EOF:
		linestring = file.readline()
		linelist = string.split(linestring)
		print linelist
		if len(linestring) == 0:
			EOF = True
		else:
			if (linelist[0] == "dt"):
				dt = float(linelist[2])
			if (linelist[0] == "num_modes"):
				num_modes = int(linelist[2])
			if (linelist[0] == "num_timesteps"):
				num_timesteps = int(linelist[2])
			if (linelist[0] == "pod_mode_dir"):
				pod_mode_dir = linelist[2]
			if (linelist[0] == "rom_coefficient_dir"):
				rom_coefficient_dir = linelist[2]
	file.close()

	return dt, num_modes, num_timesteps, pod_mode_dir, rom_coefficient_dir

#--------------------------------------------------------------------
# EOF
#--------------------------------------------------------------------
