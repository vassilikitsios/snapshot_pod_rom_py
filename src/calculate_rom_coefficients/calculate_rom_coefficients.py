#!/usr/bin/env python

#-----------------------------------------------------------------------
# Code name:	calculate_rom_coefficients.py
# Author:	Vassili Kitsios
#
#-----------------------------------------------------------------------
# Summary:
#	The code calculates the constant, linear and quadratic coefficients required for a Proper Orthogonal Decomposition (POD) Reduced Order Model (ROM). Various calibration options are available.
# 
# 	If you are using this code or elements from it please reference the following thesis:
#		Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne & l'Universite de Poitiers.
#		http://dtl.unimelb.edu.au/R/3UBLGMSICIDB6EAUJRYILBLG61QJPH79U8XSRU1FR5JQQL88XQ-00654?func=dbin-jump-full&object_id=265367&local_base=GEN01&pds_handle=GUEST
#
#	The module "input_deck_m" handles the reading of the input deck "calculate_ROM_coefficients.in".
#	The module "rom_coefficients_m" handles the calculation of the POD ROM coefficients. The associated equation numbers from Kitsios (2010) are referenced in this module next to the appropriate sections in the code.
#
#-----------------------------------------------------------------------
# Notes:
#	The calculation of the constant and linear terms require the calculation of both a volume integral and a surface integral - see equations 5.13 and 5.14 of Kitsios (2010).
#	The calculation of the quadratic terms require only the calculation of a volume integral - see equation 5.15 of Kitsios (2010).
#	For the calculation of the surface integrals one must indentify of the cells on the boundary of the domain, determine the outward pointing normal vector, the velocity and velocity derivatives on these boundary cells.
#	This has been implemented in the present code for 3D data written in VTK format. For alternate formats this process of the boundary cells and their properties should be checked.
#	However, if one is to calibrate the linear and constant coefficients, then the issue concerning the surface integrals is mute.
#
#-----------------------------------------------------------------------
# import libraries
import sys
import os
import numpy as np
import input_deck_m as input_deck
import rom_coefficients_m as rom_coefficients

#-----------------------------------------------------------------------
# Main program
print '\n-----------------------------------------------------------------------'
print 'Running ', sys.argv[0]
print '-----------------------------------------------------------------------\n'
if not os.path.exists('./results'):
        os.mkdir('./results')

#-----------------------------------------------------------------------
print "Reading input deck ..."
Re, correct_for_cell_volumes, calibrate_coefficients, interp_factor, num_snapshots, num_components, num_dimensions, num_modes, pod_mode_dir = input_deck.read_deck("calculate_rom_coefficients.in")

#-----------------------------------------------------------------------
print "\nReading data ..."

constant = np.array(np.zeros((num_modes), dtype=np.float64))
linear = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
quadratic = np.array(np.zeros((num_modes,num_modes,num_modes), dtype=np.float64))

if (calibrate_coefficients!="temporal"):
	num_points, num_boundary_faces = rom_coefficients.get_number_of_points_and_boundary_faces(pod_mode_dir,calibrate_coefficients)
	print '   Number of points = ', num_points
	print '   Number of boundary faces = ', num_boundary_faces

	# cell values interpolated to points
	cell_volume = np.array(np.zeros((num_points), dtype=np.float64))
	u = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	v = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	w = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dudx = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dudy = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dvdx = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dvdy = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dudz = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dvdz = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dwdx = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dwdy = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))
	dwdz = np.array(np.zeros((num_points,num_modes+1), dtype=np.float64))

	# boundary face variables
	norm = np.array(np.zeros((num_boundary_faces,3), dtype=np.float64))
	uF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	vF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	wF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dudxF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dudyF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dvdxF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dvdyF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dudzF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dvdzF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dwdxF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dwdyF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))
	dwdzF = np.array(np.zeros((num_boundary_faces,num_modes+1), dtype=np.float64))

	rom_coefficients.read_data_and_differentiate(pod_mode_dir,num_modes,num_points,num_boundary_faces, \
						num_dimensions,num_components,correct_for_cell_volumes, \
						u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
						uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
						cell_volume,norm,calibrate_coefficients)

#-----------------------------------------------------------------------
print "\nCalculating ROM coefficients ..."
if (calibrate_coefficients=="none"):
	rom_coefficients.calculate_constant_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
							Re,cell_volume,norm,constant)
	rom_coefficients.calculate_linear_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
							Re,cell_volume,norm,linear)
	rom_coefficients.calculate_quadratic_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							cell_volume,quadratic)
else:
	a = np.array(np.zeros((num_modes,num_snapshots), dtype=np.float64))
	t = np.array(np.zeros((num_snapshots), dtype=np.float64))
	rom_coefficients.read_temporal_modes(pod_mode_dir,num_modes,num_snapshots,a,t)

	num_steps = (num_snapshots-1)*interp_factor+1
	t_interp = np.array(np.zeros((num_snapshots*interp_factor), dtype=np.float64))
	a_interp = np.array(np.zeros((num_modes,num_steps), dtype=np.float64))
	da_dt_interp = np.array(np.zeros((num_modes,num_steps), dtype=np.float64))
	if (interp_factor==1):
		t_interp = t
		a_interp = a
		da_dt_interp[:,0:num_snapshots-1] = (a[:,1:num_snapshots] - a[:,0:num_snapshots-1]) / (t[1]-t[0])
		da_dt_interp[:,num_snapshots-1] = da_dt_interp[:,num_snapshots-2]
	else:
		rom_coefficients.interpolate_temporal_modes(num_modes,num_snapshots,interp_factor,t,a,t_interp,a_interp,da_dt_interp)
	rom_coefficients.write_temporal_modes(num_modes,num_steps,t_interp[1]-t_interp[0],a_interp,da_dt_interp)

	if (calibrate_coefficients=="constant"):
		rom_coefficients.calculate_linear_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
							Re,cell_volume,norm,linear)
		rom_coefficients.calculate_quadratic_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							cell_volume,quadratic)
		rom_coefficients.calibrate_constant_coefficients(num_modes,num_steps,a_interp,da_dt_interp,quadratic,constant,linear)
	elif (calibrate_coefficients=="constant_and_linear"):
		rom_coefficients.calculate_quadratic_coefficients(num_modes, \
							u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
							cell_volume,quadratic)
		rom_coefficients.calibrate_constant_and_linear_coefficients(num_modes,num_steps,a_interp,da_dt_interp,quadratic,constant,linear)

#-----------------------------------------------------------------------
print "\nWriting coefficients to ..."
rom_coefficients.write_constant_coefficients(num_modes,constant)
rom_coefficients.write_linear_coefficients(num_modes,linear)
rom_coefficients.write_quadratic_coefficients(num_modes,quadratic)

print '\n-----------------------------------------------------------------------'
print 'Code completed.'
print '-----------------------------------------------------------------------\n'

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
