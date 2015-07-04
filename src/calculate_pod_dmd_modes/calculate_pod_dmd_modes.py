#!/usr/bin/env python

#-----------------------------------------------------------------------
# Code name:	calculate_pod_dmd_modes.py
# Summary:	calculates POD and DMD/Koopman modes from a series of snapshots
# Author:	Vassili Kitsios

#-----------------------------------------------------------------------
# Summary:
#       The code calculates the proper orthogonal decomposition (POD) modes, (also known as principal component analysis, or empirical orthgonal functions, or singular value decomposition), and Koopman modes (also known as dynamic mode decomposiion, or principal oscillation patterns).
#	The POD modes are later to be used to build a reduced order model using these modes as basis.
# 
#       If you are using this code or elements from it please reference the following thesis:
#               Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne & l'Universite de Poitiers.
#               http://dtl.unimelb.edu.au/R/3UBLGMSICIDB6EAUJRYILBLG61QJPH79U8XSRU1FR5JQQL88XQ-00654?func=dbin-jump-full&object_id=265367&local_base=GEN01&pds_handle=GUEST
#
#       The module "input_deck_m" handles the reading of the input deck "calculate_pod_dmd_modes.in".
#       The module "pod_m" handles the calculation of the POD modes.
#	The module "dmd_m" handles the calculation of the DMD/Koopman modes.
#	The module "io_vtk_m" handles the input/output of the fields written in the visualisation tool kit (VTK) file format.

#-----------------------------------------------------------------------
# import libraries
import math
import sys
import os

import numpy as np
import scipy as Sci
import scipy.linalg as linalg
import vtk
import vtk.util.numpy_support as VN

import input_deck_m as input_deck
import pod_m as pod
import dmd_m as dmd
import io_vtk_m as io_vtk

#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Running', sys.argv[0]
print '   This code calculates DMD and POD modes from a series of snapshots.'

if not os.path.exists('./results'):
	os.mkdir('./results')

print 'Reading input deck ...'
snapshot_dir, snapshot_list_file_name, dt, cell_volume_file_name, correct_for_cell_volumes, subtract_mean, \
	num_components, x_min, x_max, y_min, y_max, z_min, z_max, num_modes_trunc, num_modes_to_write, write_matrices, \
	test_POD_orthogonality, test_DMD_biorthogonality, test_stochastic_mode_orthogonality, tol_CN, tol_ortho, var_name,\
	restart_flag, restart_dir = input_deck.read_deck("calculate_pod_dmd_modes.in")

print '\n--------------------------------------------------------------'
print 'Building snapshot matrix ...'
fin = open(snapshot_dir + snapshot_list_file_name,'r')
file_list = list()
file_list = fin.readlines()
fin.close()
num_snapshots = len(file_list) - 1
print '   Number of snapshots = ', num_snapshots

grid = vtk.vtkUnstructuredGrid()
if (correct_for_cell_volumes=="true"):
	io_vtk.get_grid_template(x_min,x_max,y_min,y_max,z_min,z_max,snapshot_dir+cell_volume_file_name,grid)
	num_points = grid.GetNumberOfPoints()
	cell_volume = np.array(np.zeros((num_points), dtype=np.float64))
	cell_volume = VN.vtk_to_numpy(grid.GetPointData().GetScalars("cell_volume"))
else:
	io_vtk.get_grid_template(x_min,x_max,y_min,y_max,z_min,z_max,snapshot_dir+file_list[0].strip(),grid)
	num_points = grid.GetNumberOfPoints()
	cell_volume = np.array(np.ones((num_points), dtype=np.float64))
print '   Number of points = ', num_points

A = np.array(np.zeros((num_points*num_components,num_snapshots+1), dtype=np.float64))
io_vtk.read_data_and_build_snapshot_matrix(x_min,x_max,y_min,y_max,z_min,z_max,\
					snapshot_dir,file_list,num_snapshots,num_points,num_components,var_name,A)

first_snapshot = np.array(np.zeros((num_points*num_components), dtype=np.float64))
first_snapshot = A[:,0]

print '\n   Calculating mean field ...'
mean_field = np.array(np.zeros((num_points*num_components), dtype=np.float64))
for j in range(0,num_snapshots):
	mean_field[:] = mean_field[:] + A[:,j].copy()
mean_field[:] = mean_field[:] / num_snapshots

if (subtract_mean=="true"):
	print '\n   Subtracting mean from instantaneous fields ...'
	for j in range(0,num_snapshots+1):
		A[:,j] = A[:,j] - mean_field[:]


#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Calculating POD modes ...'

C = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.float64))
if (restart_flag=='false'):
	print '\n   Calculating covariance matrix ...'
	pod.calculate_correlation_matrix(num_snapshots,num_points,num_components,correct_for_cell_volumes,cell_volume,A,C)
elif (restart_flag=='true'):
	pod.read_correlation_matrix(restart_dir,num_snapshots,C)

print '\n   Solving eigenvalue problem ...'
energy = np.array(np.zeros((num_snapshots), dtype=np.complex64))
temporal_modes = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.complex64))
energy,temporal_modes = linalg.eig(C)

num_valid_modes = 0
while ( (energy[num_valid_modes].real/energy[0].real>pow(tol_CN,2.0)) and (num_valid_modes<num_snapshots-1) \
		and(energy[num_valid_modes].real>0.0) ):
	num_valid_modes += 1
if ( (energy[num_valid_modes].real/energy[0].real > pow(tol_CN,2.0)) and (energy[num_valid_modes].real > 0.0) ):
	num_valid_modes += 1
print '      Number of valid POD modes with positive energies = ', num_valid_modes
if ( (num_modes_trunc < 0) or (num_modes_trunc > num_valid_modes) ):
	num_modes_trunc = num_valid_modes

print '\n   Scaling temporal modes ...'
for j in range(0,num_valid_modes):
	temporal_mode_mag = sum(temporal_modes[:,j].real * temporal_modes[:,j].real)/num_snapshots
	temporal_modes[:,j] = temporal_modes[:,j] * np.sqrt(energy[j].real/temporal_mode_mag)
initial_conditions = np.array(np.zeros((num_valid_modes), dtype=np.float64))
initial_conditions = temporal_modes[0,:].real

print '\n   Calculating truncated spatial modes ...'
energy_trunc_inv = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.float64))
energy_trunc_inv = np.diag(np.ones(num_modes_trunc)/energy[0:num_modes_trunc].real,0)
spatial_modes_trunc = np.array(np.zeros((num_components*num_points,num_modes_trunc), dtype=np.float64))
spatial_modes_trunc = np.dot( np.dot(A[:,0:num_snapshots], temporal_modes[:,0:num_modes_trunc].real), energy_trunc_inv) / num_snapshots

if (test_POD_orthogonality=="true"):
	num_ortho_modes = pod.test_orthogonality(num_modes_trunc,num_components,num_points,\
				spatial_modes_trunc,cell_volume,correct_for_cell_volumes,tol_ortho)
	print '      Number of POD modes that satisfy orthogonality = ', num_ortho_modes
	if (num_ortho_modes<num_modes_trunc):
		num_modes_trunc = num_ortho_modes
		del energy_trunc_inv
		energy_trunc_inv = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.float64))
		energy_trunc_inv = np.diag(np.ones(num_modes_trunc)/energy[0:num_modes_trunc].real,0)
		del spatial_modes_trunc
		spatial_modes_trunc = np.array(np.zeros((num_components*num_points,num_modes_trunc), dtype=np.float64))
		spatial_modes_trunc = np.dot( np.dot(A[:,0:num_snapshots], temporal_modes[:,0:num_modes_trunc].real), \
						energy_trunc_inv) / num_snapshots
print '      Number of POD modes used to truncate the DMD eigenvalue problem = ', num_modes_trunc


#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Writing POD output ...'
if ( (write_matrices=='POD') or (write_matrices=='both') ):
	pod.write_correlation_matrix(num_snapshots,C)
del C

io_vtk.write_mean_field(num_points,num_components,grid,mean_field)
del mean_field

pod.write_eigenvalues(num_valid_modes,num_snapshots,energy,'./results/POD.eigenvalues.dat')
pod.write_initial_conditions(num_valid_modes,initial_conditions)
pod.write_temporal_modes(num_valid_modes,num_snapshots,dt,temporal_modes)

num_POD_modes_to_write = num_modes_to_write
if ( (num_modes_to_write < 0) or (num_modes_to_write > num_modes_trunc) ):
	num_POD_modes_to_write = num_modes_trunc
print '\n   Number of spatial POD modes to write = ', num_POD_modes_to_write
io_vtk.write_spatial_POD_modes(num_POD_modes_to_write,num_points,num_components,grid,spatial_modes_trunc)

if (correct_for_cell_volumes=="true"):
	print '\n--------------------------------------------------------------'
	print 'Exiting program, only calculate DMD modes when there is no cell volume correction.'
	print '\nCode completed.'
	print '--------------------------------------------------------------\n'
	sys.exit()


#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Calculating DMD modes ...'

print '\n   Calculating linear operator ...'
DMD_trunc = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.float64))
DMD_trunc = np.dot( np.dot( np.dot(spatial_modes_trunc.T, A[:,1:num_snapshots+1]) , temporal_modes[:,0:num_modes_trunc].real ) , \
					energy_trunc_inv) / num_snapshots
for k in range(0,num_modes_trunc):
	DMD_trunc[k,k] = DMD_trunc[k,k] - 1.0
DMD_trunc = DMD_trunc / dt
del energy_trunc_inv
del A

print '\n   Solving eigenvalue problem ...'
w = np.array(np.zeros((num_modes_trunc), dtype=np.complex64))
vr = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
vl = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
w,vl,vr = linalg.eig(DMD_trunc, left=True, right=True, overwrite_a=False, overwrite_b=False)
w = w*1j

print '\n   Sorting eigenvalues ...'
mode_index = np.zeros((num_modes_trunc), dtype=np.integer)
num_valid_DMD_modes = dmd.sort_eigenvalues(num_modes_trunc,w,vl,vr,mode_index)
print '      Number of valid DMD modes ', num_valid_DMD_modes

print '\n   Scaling spatial modes ...'
for k in range(0,num_modes_trunc):
	mag = np.dot(vr[:,k].T.conj(),vr[:,k])		# scale consistent with POD spatial modes
	vr[:,k] = vr[:,k] / np.sqrt(mag)
	mag = np.dot(vl[:,k].T.conj(),vr[:,k])		# for bi-orthogonality
	vl[:,k] = vl[:,k] / mag.conj()
if (test_DMD_biorthogonality=="true"):
	dmd.test_biorthogonality(num_modes_trunc,vr,vl)

print '\n   Calculating energy in each mode ...'
energy_DMD = np.array(np.zeros((num_modes_trunc), dtype=np.float64))
energy_DMD = np.real(np.diag(np.dot( np.dot(vl.T.conj(),np.diag(energy[0:num_modes_trunc].real,0)), vr)))

print '\n   Calculating phase of each mode ...'
initial_conditions_DMD = np.array(np.zeros((num_modes_trunc), dtype=np.complex64))
initial_conditions_DMD = np.dot(np.dot(spatial_modes_trunc, vl).T.conj(),first_snapshot)
phase_DMD = np.array(np.zeros((num_modes_trunc), dtype=np.float64))
for k in range(0,num_modes_trunc):
	phase_DMD[k] = math.atan2(initial_conditions_DMD[k].imag,initial_conditions_DMD[k].real)

print '\n   Scaling temporal modes ...'
temporal_modes_DMD = np.array(np.zeros((num_snapshots,num_valid_DMD_modes), dtype=np.complex64))
dmd.calculate_temporal_modes(num_valid_DMD_modes,num_snapshots,dt,w,mode_index,energy_DMD,phase_DMD,temporal_modes_DMD)
#initial_conditions_DMD = temporal_modes_DMD[0,:]

print '\n   Calculating variance of stochastic noise for each DMD mode ...'
DMD_trunc_DMD_space = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
DMD_trunc_DMD_space = np.dot( np.dot(vl.T.conj(), DMD_trunc), vr)
noise_DMD = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
noise_DMD = -np.dot( DMD_trunc_DMD_space, np.dot(temporal_modes_DMD.T.conj(), temporal_modes_DMD) ) / num_snapshots
noise_DMD = noise_DMD + noise_DMD.T.conj() 


#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Calculating stochastic modes ...'

print '\n   Calculating stochastic noise projected onto the POD modes ...'
noise_POD = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.float64))
noise_POD = -np.dot( DMD_trunc, np.diag(energy[0:num_modes_trunc].real,0) )
noise_POD = noise_POD + noise_POD.T.conj() 

print '\n   Solving eigenvalue problem ...'
noise_energy = np.array(np.zeros((num_snapshots), dtype=np.complex64))
noise_evector = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.complex64))
noise_energy,noise_evector = linalg.eig(noise_POD)

num_valid_stochastic_modes = 0
while ( (noise_energy[num_valid_stochastic_modes].real/noise_energy[0].real > pow(tol_CN,2.0)) and (num_valid_stochastic_modes < num_modes_trunc-1) and (noise_energy[num_valid_stochastic_modes].real > 0.0) ):
	num_valid_stochastic_modes += 1
if ( (noise_energy[num_valid_stochastic_modes].real/noise_energy[0].real > pow(tol_CN,2.0)) and (noise_energy[num_valid_stochastic_modes].real > 0.0) ):
	num_valid_stochastic_modes += 1
print '      Number of valid Stochastic noise modes with positive energies = ', num_valid_stochastic_modes

print '\n   Scaling spatial modes ...'
for k in range(0,num_valid_stochastic_modes):
	mag = np.dot(noise_evector[:,k].T,noise_evector[:,k])		# scale consistent with POD spatial modes
	noise_evector[:,k] = noise_evector[:,k] / np.sqrt(mag)

if (test_stochastic_mode_orthogonality=="true"):
	dmd.test_stochastic_orthogonality(num_valid_stochastic_modes,noise_evector)


#-----------------------------------------------------------------------
print '\n--------------------------------------------------------------'
print 'Writing DMD output ...'
if ( (write_matrices=='DMD') or (write_matrices=='both') ):
	dmd.write_real_matrix(num_modes_trunc,DMD_trunc,'DMD.drain.pod_space.dat')
	dmd.write_complex_matrix(num_modes_trunc,DMD_trunc_DMD_space,'DMD.drain.dmd_space.dat')
	dmd.write_complex_matrix(num_modes_trunc,noise_DMD,'DMD.stochastic_noise.dmd_space.dat')
	dmd.write_real_matrix(num_modes_trunc,noise_POD,'DMD.stochastic_noise.pod_space.dat')

dmd.write_eigenvalues(num_valid_DMD_modes,w,mode_index,energy_DMD,phase_DMD,dt,noise_DMD,DMD_trunc_DMD_space)
dmd.write_initial_conditions(num_valid_DMD_modes,mode_index,initial_conditions_DMD)
dmd.write_temporal_modes(num_valid_DMD_modes,num_snapshots,dt,mode_index,temporal_modes_DMD)
pod.write_eigenvalues(num_valid_stochastic_modes,num_modes_trunc,noise_energy,'./results/DMD.eigenvalues_stochastic.dat')

num_stochastic_modes_to_write = num_modes_to_write
if ( (num_modes_to_write < 0) or (num_modes_to_write > num_valid_stochastic_modes) ):
	num_stochastic_modes_to_write = num_valid_stochastic_modes
print '\n   Number of spatial stochastic modes to write = ', num_stochastic_modes_to_write
io_vtk.write_spatial_stochastic_modes(num_stochastic_modes_to_write,num_points,num_components,grid,spatial_modes_trunc,noise_evector[:,:].real)

num_DMD_modes_to_write = num_modes_to_write
if ( (num_modes_to_write < 0) or (num_modes_to_write > num_valid_DMD_modes) ):
	num_DMD_modes_to_write = num_valid_DMD_modes
print '\n   Number of spatial DMD modes to write = ', num_DMD_modes_to_write
io_vtk.write_spatial_DMD_modes(num_DMD_modes_to_write,num_points,num_components,grid,spatial_modes_trunc,vr,vl,mode_index)

print '\nCode completed.'
print '--------------------------------------------------------------\n'

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------

