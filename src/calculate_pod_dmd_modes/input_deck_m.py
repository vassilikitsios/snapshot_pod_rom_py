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
			if (linelist[0] == "snapshot_dir"):
				snapshot_dir = linelist[2]
			if (linelist[0] == "snapshot_list_file_name"):
				snapshot_list_file_name = linelist[2]
			if (linelist[0] == "dt"):
				dt = float(linelist[2])
			if (linelist[0] == "cell_volume_file_name"):
				cell_volume_file_name = linelist[2]
			if (linelist[0] == "correct_for_cell_volumes"):
				correct_for_cell_volumes = linelist[2]
			if (linelist[0] == "subtract_mean"):
				subtract_mean = linelist[2]
			if (linelist[0] == "num_components"):
				num_components = int(linelist[2])
			if (linelist[0] == "x_min"):
				x_min = float(linelist[2])
			if (linelist[0] == "x_max"):
				x_max = float(linelist[2])
			if (linelist[0] == "y_min"):
				y_min = float(linelist[2])
			if (linelist[0] == "y_max"):
				y_max = float(linelist[2])
			if (linelist[0] == "z_min"):
				z_min = float(linelist[2])
			if (linelist[0] == "z_max"):
				z_max = float(linelist[2])
			if (linelist[0] == "num_modes_trunc"):
				num_modes_trunc = int(linelist[2])
			if (linelist[0] == "num_modes_to_write"):
				num_modes_to_write = int(linelist[2])
			if (linelist[0] == "write_matrices"):
				write_matrices = linelist[2]
			if (linelist[0] == "test_POD_orthogonality"):
				test_POD_orthogonality = linelist[2]
			if (linelist[0] == "test_DMD_biorthogonality"):
				test_DMD_biorthogonality = linelist[2]
			if (linelist[0] == "test_stochastic_mode_orthogonality"):
				test_stochastic_mode_orthogonality = linelist[2]
			if (linelist[0] == "tol_CN"):
				tol_CN = float(linelist[2])
			if (linelist[0] == "tol_ortho"):
				tol_ortho = float(linelist[2])
			if (linelist[0] == "restart_dir"):
				restart_dir = linelist[2]
			if (linelist[0] == "restart_flag"):
				restart_flag = linelist[2]
			if (linelist[0] == "var_name"):
				var_name = linelist[2]

	if ( (correct_for_cell_volumes != 'true') and (correct_for_cell_volumes != 'false') ):
	        print '\nFAILURE'
        	print "   correct_for_cell_volumes = 'true' or 'false'"
        	sys.exit()
		 
	if ( (subtract_mean != 'true') and (subtract_mean != 'false') ):
	        print '\nFAILURE'
        	print "   subtract_mean = 'true' or 'false'"
        	sys.exit()

	if ( (restart_flag != 'true') and (restart_flag != 'false') ):
	        print '\nFAILURE'
        	print "   restart_flag = 'true' or 'false'"
        	sys.exit()

	if ( (write_matrices != 'POD') and (write_matrices != 'DMD') and (write_matrices != 'both') and (write_matrices != 'none') ):
	        print '\nFAILURE'
        	print "   write_matrices = 'POD' or 'DMD' or 'both' or 'none'"
        	sys.exit()

	if ( (test_POD_orthogonality != 'true') and (test_POD_orthogonality != 'false') ):
	        print '\nFAILURE'
        	print "   test_POD_orthogonality = 'true' or 'false'"
        	sys.exit()

	if ( (test_DMD_biorthogonality != 'true') and (test_DMD_biorthogonality != 'false') ):
	        print '\nFAILURE'
        	print "   test_DMD_biorthogonality = 'true' or 'false'"
        	sys.exit()

	if ( (test_stochastic_mode_orthogonality != 'true') and (test_stochastic_mode_orthogonality != 'false') ):
	        print '\nFAILURE'
        	print "   test_stochastic_mode_orthogonality = 'true' or 'false'"
        	sys.exit()

	if ( (num_components < 1) or (num_components > 3) ):
	        print '\nFAILURE'
        	print "   num_components = 1, 2 or 3"
        	sys.exit()

	if (tol_CN < 0):
	        print '\nFAILURE'
        	print "   tol_CN > 0"
        	sys.exit()

	if (tol_ortho < 0):
	        print '\nFAILURE'
        	print "   tol_ortho > 0"
        	sys.exit()

	file.close()

	return snapshot_dir, snapshot_list_file_name, dt, cell_volume_file_name, correct_for_cell_volumes, subtract_mean, num_components, x_min, x_max, y_min, y_max, z_min, z_max, num_modes_trunc, num_modes_to_write, write_matrices, test_POD_orthogonality, test_DMD_biorthogonality, test_stochastic_mode_orthogonality, tol_CN, tol_ortho, var_name, restart_flag, restart_dir

#--------------------------------------------------------------------
# EOF
#--------------------------------------------------------------------

