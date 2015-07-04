#!/usr/bin/env python

#-----------------------------------------------------------------------
# Code name:    integrate_rom.py
# Author:       Vassili Kitsios
#
#-----------------------------------------------------------------------
# Summary:
#       The code runs steps forward in time the Proper Orthogonal Decomposition (POD) Reduced Order Model (ROM) using the coefficients calculated from the code "calculate_rom_coefficients.py".
# 
#       If you are using this code or elements from it please reference the following thesis:
#               Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne & l'Universit\'e de Poitiers.
#               http://dtl.unimelb.edu.au/R/3UBLGMSICIDB6EAUJRYILBLG61QJPH79U8XSRU1FR5JQQL88XQ-00654?func=dbin-jump-full&object_id=265367&local_base=GEN01&pds_handle=GUEST
#
#       The module "input_deck_m" handles the reading of the input deck "calculate_ROM_coefficients.in".
#       The module "rom_m" handles the input/output operations and also the calculation of the time derivative required to setp forward the POD ROM. The associated equation numbers from Kitsios (2012) are referenced in this module next to the appropriate sections in the code.
#
#-----------------------------------------------------------------------
# Notes:
#       Currently the time stepping is done using a simple Euler method, but could be improved using a Runge Kutta scheme.
#
#-----------------------------------------------------------------------
# import libraries
import sys
import os
import numpy as np
import input_deck_m as input_deck
import rom_m as rom

#-----------------------------------------------------------------------
# Main program

print '\n-----------------------------------------------------------------------'
print 'Running ', sys.argv[0]
print '-----------------------------------------------------------------------\n'

if not os.path.exists('./results'):
        os.mkdir('./results')

#-----------------------------------------------------------------------
print "Reading input deck ...\n"
dt, num_modes, num_timesteps, pod_mode_dir, rom_coefficient_dir = input_deck.read_deck("integrate_rom.in")

a = np.array(np.zeros((num_modes,num_timesteps), dtype=np.float64))
da_dt = np.array(np.zeros((num_modes,num_timesteps), dtype=np.float64))
constant = np.array(np.zeros((num_modes), dtype=np.float64))
linear = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
quadratic = np.array(np.zeros((num_modes,num_modes,num_modes), dtype=np.float64))

#-----------------------------------------------------------------------
print "Reading initial conditions and coefficients ...\n"
rom.read_initial_conditions(pod_mode_dir,num_modes,a[:,0])
rom.read_constant_coefficients(rom_coefficient_dir,num_modes,constant)
rom.read_linear_coefficients(rom_coefficient_dir,num_modes,linear)
rom.read_quadratic_coefficients(rom_coefficient_dir,num_modes,quadratic)

#-----------------------------------------------------------------------
print "\nPerforming time integration ..."
for t in range(0,num_timesteps-1):
	if (t % 100 == 0):
		print "      timestep ", t, " of ", num_timesteps
	rom.calculate_time_derivative(num_modes,constant,linear,quadratic,a[:,t],da_dt[:,t])
	a[:,t+1] = a[:,t] + da_dt[:,t] * dt		# euler method, but could improve with Runge Kutta
rom.calculate_time_derivative(num_modes,constant,linear,quadratic,a[:,t+1],da_dt[:,t+1])

#-----------------------------------------------------------------------
print "\nWriting results ..."
rom.write_fourier_coefficients(num_modes,num_timesteps+1,dt,a,da_dt)

print '\n-----------------------------------------------------------------------'
print 'Code completed.'
print '-----------------------------------------------------------------------\n'
#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
