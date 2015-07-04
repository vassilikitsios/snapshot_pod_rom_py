# module pod_m

#-----------------------------------------------------------------------
# import libraries
import numpy as np
import math
import string

#-----------------------------------------------------------------------
def calculate_correlation_matrix(num_snapshots,num_points,num_components,correct_for_cell_volumes,cell_volume,A,C):
	C[:,:] = 0.0
	if (correct_for_cell_volumes=="false"):
		C[:,:] = np.dot( A[:,0:num_snapshots].T, A[:,0:num_snapshots]) / num_snapshots
	elif (correct_for_cell_volumes=="true"):
		for i in range(0,num_snapshots):
			print "      Correlations with snapshot ", i+1, " of ", num_snapshots
			for j in range(0,num_snapshots):
				if (j>=i):
					for k in range(0,num_components):
						C[i,j] = C[i,j] + np.dot( A[k*num_points:(k+1)*num_points,i].T * cell_volume[:], A[k*num_points:(k+1)*num_points,j] ) / num_snapshots
				if (j<i):
					C[i,j] = C[j,i]

#-----------------------------------------------------------------------
def test_orthogonality(num_modes_trunc,num_components,num_points,spatial_modes_trunc,cell_volume,correct_for_cell_volumes,tol):
	print '\n   Testing orthogonality of POD spatial modes ...'
	ortho_POD = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
	if (correct_for_cell_volumes=="true"):
		for i in range(0,num_modes_trunc):
			for j in range(0,num_modes_trunc):
				if (j>=i):
					for k in range(0,num_components):
						ortho_POD[i,j] = ortho_POD[i,j] \
						+ np.dot( spatial_modes_trunc[k*num_points:(k+1)*num_points,i].T * cell_volume[:],\
								 spatial_modes_trunc[k*num_points:(k+1)*num_points,j] )
	else:
		ortho_POD[:,:] = np.dot(spatial_modes_trunc.T,spatial_modes_trunc)

	print ortho_POD
	num_ortho_modes = num_modes_trunc
	for i in range(0,num_modes_trunc):
		j = num_modes_trunc - 1 - i
		if (abs(ortho_POD[j,j]-1.0)>tol):
			num_ortho_modes = j

	return num_ortho_modes

#-----------------------------------------------------------------------
def write_correlation_matrix(num_snapshots,C):
	print '\n   Writing covariance matrix to'
	filename = './results/POD.covariance.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode i, mode j, covariance\n')
	file.write('#\n')
	for i in range(0,num_snapshots):
		for j in range(0,num_snapshots):
			file.write('%4.1d %4.1d %23.15e\n' % (i+1, j+1, C[i,j]) )
		file.write('\n')
	file.close()

#-----------------------------------------------------------------------
def read_correlation_matrix(pod_mode_dir,num_snapshots,C):
        print '\n   Reading correlation matrix ...'
        filename = pod_mode_dir+'/POD.covariance.dat'
        file = open(filename,'r')
        for i in range(0,3):
                linestring = file.readline()
	for i in range(0,num_snapshots):
        	for j in range(0,num_snapshots):
                        linestring = file.readline()
                        linelist = string.split(linestring)
                        C[i,j]= float(linelist[2])
                linestring = file.readline()
        file.close()

#-----------------------------------------------------------------------
def write_eigenvalues(num_valid_modes,num_snapshots,energy,filename):
	cumulative_energy = np.array(np.zeros((num_valid_modes), dtype=np.float64))
	cumulative_energy[0] = energy[0].real
	for i in range(1,num_valid_modes):
		cumulative_energy[i] = cumulative_energy[i-1] + energy[i].real
	total_energy = cumulative_energy[num_valid_modes-1]

	print '\n   Writing energies to'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, energy, cumulative, percenterage energy, percentage cumulative, condition number (absolute value if negative)\n')
	file.write('#		Note: cummulative energies are set to zero after first negative energy')
	file.write('#\n')
	for i in range(0,num_valid_modes):
		file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (i+1, energy[i].real, cumulative_energy[i], energy[i].real/total_energy*100.0, cumulative_energy[i]/total_energy*100.0, math.sqrt(energy[i].real / energy[0].real)) )
	for i in range(num_valid_modes,num_snapshots):
		file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (i+1, energy[i].real, 0.0, energy[i].real/total_energy*100.0, 0.0, math.sqrt(abs(energy[i].real / energy[0].real))) )
	file.close()

#-----------------------------------------------------------------------
def write_initial_conditions(num_valid_modes,initial_conditions):
	print '\n   Writing initial conditions to'
	filename = './results/POD.initial_conditions.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, amplitude\n')
	file.write('#\n')
	for i in range(0,num_valid_modes):
		file.write('%4.1d %18.10e\n' % (i+1, initial_conditions[i]) )
	file.close()

#-----------------------------------------------------------------------
def write_temporal_modes(num_valid_modes,num_snapshots,dt,temporal_modes):
	print '\n   Writing temporal modes to'
	for j in range(0,num_valid_modes):
		j_index = j+1
		output_file_name = './results/POD.temporal_mode_' + '%04d'%j_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'w')
		file.write('#\n')
		file.write('# time, amplitude\n')
		file.write('#\n')
		for i in range(0,num_snapshots):
			file.write('%18.10e %18.10e\n' % (i*dt, temporal_modes[i,j].real) )
		file.close()

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
