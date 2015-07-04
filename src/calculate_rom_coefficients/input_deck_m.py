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
			if (linelist[0] == "Re"):
				Re = float(linelist[2])
			if (linelist[0] == "correct_for_cell_volumes"):
				correct_for_cell_volumes = linelist[2]
			if (linelist[0] == "calibrate_coefficients"):
				calibrate_coefficients = linelist[2]
			if (linelist[0] == "interp_factor"):
				interp_factor = int(linelist[2])
			if (linelist[0] == "num_components"):
				num_components = int(linelist[2])
			if (linelist[0] == "num_dimensions"):
				num_dimensions = int(linelist[2])
			if (linelist[0] == "num_snapshots"):
				num_snapshots = int(linelist[2])
			if (linelist[0] == "num_modes"):
				num_modes = int(linelist[2])
			if (linelist[0] == "pod_mode_dir"):
				pod_mode_dir = linelist[2]

	if ( (calibrate_coefficients != 'none') and (calibrate_coefficients != 'constant') and (calibrate_coefficients != 'constant_and_linear') ):
	        print '\nFAILURE'
        	print "   calibrate_coefficients = 'none' or 'constant' or 'constant_and_linear'"
        	sys.exit()
		 
	if ( (correct_for_cell_volumes != 'true') and (correct_for_cell_volumes != 'false') ):
	        print '\nFAILURE'
        	print "   correct_for_cell_volumes = 'true' or 'false'"
        	sys.exit()
		 
	if ( (num_components > 3) or (num_components < 1) ):
	        print '\nFAILURE'
        	print "   1<= num_components <= 3"
        	sys.exit()

	if (num_snapshots < 1):
	        print '\nFAILURE'
        	print "   num_snapshots > 1"
        	sys.exit()

	if ( (num_components > 3) or (num_components < 1) ):
	        print '\nFAILURE'
        	print "   1 <= num_dimensions <= 3"
        	sys.exit()

	if (interp_factor < 1):
	        print '\nFAILURE'
        	print "   interp_factor >= 1"
        	sys.exit()

	file.close()

	return Re, correct_for_cell_volumes, calibrate_coefficients, interp_factor, num_snapshots, num_components, num_dimensions, num_modes, pod_mode_dir

#--------------------------------------------------------------------

