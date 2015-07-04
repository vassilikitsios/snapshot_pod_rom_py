# module rom_m.py

#-----------------------------------------------------------------------
# import libraries
import numpy as np
import string

#-----------------------------------------------------------------------
def read_interpolated_temporal_modes(rom_coefficient_dir,num_modes,num_timesteps,a_interp,da_dt_interp):
	print '\n   Reading temporal modes ...'
	for j in range(0,num_modes):
		j_index = j+1
		filename = rom_coefficient_dir+'/POD_interp.temporal_mode_' + '%04d'%j_index + '.dat'
		file = open(filename,'r')
		for i in range(0,3):
			linestring = file.readline()
		for i in range(0,num_timesteps):
			linestring = file.readline()
			linelist = string.split(linestring)
			a_interp[j,i]= float(linelist[1])
			da_dt_interp[j,i]= float(linelist[2])
		file.close()

#-----------------------------------------------------------------------
def read_initial_conditions(pod_mode_dir,num_modes,initial_conditions):
	print '\n   Reading initial conditions ...'
	filename = pod_mode_dir+'/POD.initial_conditions.dat'
	print '      ', filename
	file = open(filename,'r')
	for i in range(0,3):
		linestring = file.readline()
	for i in range(0,num_modes):
		linestring = file.readline()
		linelist = string.split(linestring)
		initial_conditions[i]= float(linelist[1])
	file.close()

#-----------------------------------------------------------------------
def read_constant_coefficients(rom_coefficient_dir,num_modes,constant):
	print '\n   Reading constant coefficients ...'
	filename = rom_coefficient_dir+'/ROM.constant_coefficients.dat'
	print '      ', filename
	file = open(filename,'r')
	for i in range(0,3):
		linestring = file.readline()
	for i in range(0,num_modes):
		linestring = file.readline()
		linelist = string.split(linestring)
		constant[i]= float(linelist[1])
	file.close()

#-----------------------------------------------------------------------
def read_linear_coefficients(rom_coefficient_dir,num_modes,linear):
	print '\n   Reading linear coefficients ...'
	filename = rom_coefficient_dir+'/ROM.linear_coefficients.dat'
	print '      ', filename
	file = open(filename,'r')
	for i in range(0,3):
		linestring = file.readline()
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			linestring = file.readline()
			linelist = string.split(linestring)
			linear[i,j]= float(linelist[2])
		linestring = file.readline()
	file.close()

#-----------------------------------------------------------------------
def read_quadratic_coefficients(rom_coefficient_dir,num_modes,quadratic):
	print '\n   Reading quadratic coefficients ...'
	filename = rom_coefficient_dir+'/ROM.quadratic_coefficients.dat'
	print '      ', filename
	file = open(filename,'r')
	for i in range(0,3):
		linestring = file.readline()
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				linestring = file.readline()
				linelist = string.split(linestring)
				quadratic[i,j,k]= float(linelist[3])
			linestring = file.readline()
	file.close()

#-----------------------------------------------------------------------
def calculate_time_derivative(num_modes,constant,linear,quadratic,a,da_dt):
	# Equation 5.29 of Kitsios (2012)
	for i in range(0,num_modes):
		quadratic_sum = 0.0
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				quadratic_sum = quadratic_sum + a[j]*a[k]*quadratic[i][j][k]
		linear_sum = 0.0
		for j in range(0,num_modes):
			linear_sum = linear_sum + a[j]*linear[i][j]
		da_dt[i] = linear_sum + quadratic_sum + constant[i]

#-----------------------------------------------------------------------
def write_fourier_coefficients(num_modes,num_timesteps,dt,a,da_dt):
        print '\n   Writing ROM fourier coefficients to'
        for j in range(0,num_modes):
                j_index = j+1
                output_file_name = './results/ROM.temporal_mode_' + '%04d'%j_index + '.dat'
                print '      ', output_file_name
                file = open(output_file_name,'w')
                file.write('#\n')
                file.write('# time, amplitude, da_dt\n')
                file.write('#\n')
                for i in range(0,num_timesteps-1):
                        file.write('%18.10e %18.10e %18.10e\n' % (i*dt, a[j,i], da_dt[j,i]) )
                file.close()

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
