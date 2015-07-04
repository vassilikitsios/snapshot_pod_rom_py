# module dmd_m

#-----------------------------------------------------------------------
# import libraries
import numpy as np
import math

#-----------------------------------------------------------------------
def sort_eigenvalues(num_modes_trunc,w,vl,vr,mode_index):
	w_imag_sorted = np.array(np.zeros((num_modes_trunc), dtype=np.float64))
	num_valid_DMD_modes = 0
	for k in range(0,num_modes_trunc):
		mode_index[k] = k
		if ( (math.isnan(w[k].real)) or (math.isnan(w[k].imag)) ):
			w_imag_sorted[k] = -1.0e10
			vl[:,k] = 0.0
			vr[:,k] = 0.0
		else:
			w_imag_sorted[k] = np.imag(w[k])
			num_valid_DMD_modes += 1
	w_imag_sorted, mode_index[:] = zip(*sorted(zip(w_imag_sorted, mode_index),reverse=True))
	return num_valid_DMD_modes

#-----------------------------------------------------------------------
def test_biorthogonality(num_modes_trunc,vr,vl):
	DMD_product = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.complex64))
	print '\n   Amplitude of direct spatial modes ...'
	DMD_product = np.dot(vr.T.conj(),vr)
	print DMD_product
	print '\n   Amplitude of adjoint spatial modes ...'
	DMD_product = np.dot(vl.T.conj(),vl)
	print DMD_product
	print '\n   Testing biorthogonality of DMD eigenvectors in POD space ...'
	DMD_product = np.dot(vl.T.conj(),vr)
	#for k in range(0,num_modes_trunc):
	#	print k, DMD_product[k,k]
	print DMD_product

#-----------------------------------------------------------------------
def test_stochastic_orthogonality(N,v):
	product = np.array(np.zeros((N,N), dtype=np.float64))
	print '\n   Testing orthogonality of DMD stochastic eigenvectors in POD space ...'
	product = np.dot(v.T.conj(),v)
	#for k in range(0,N):
	#	print k, product[k,k]
	print product

#-----------------------------------------------------------------------
def write_real_matrix(num_modes_trunc,DMD_trunc,filename):
	print '\n   Writing real matrix to'
	filename_out = './results/' + filename
	print '      ', filename_out
	file = open(filename_out,'w')
	file.write('#\n')
	file.write('# mode i, mode j, real\n')
	file.write('#\n')
	for i in range(0,num_modes_trunc):
		for j in range(0,num_modes_trunc):
			file.write('%4.1d %4.1d %18.10e\n' % (i+1, j+1, np.real(DMD_trunc[i,j])) )
		file.write('\n')
	file.close()

#-----------------------------------------------------------------------
def write_complex_matrix(num_modes_trunc,noise_DMD,filename):
	print '\n   Writing complex matrix to'
	filename_out = './results/' + filename
	print '      ', filename_out
	file = open(filename_out,'w')
	file.write('#\n')
	file.write('# mode i, mode j, real, imag\n')
	file.write('#\n')
	for i in range(0,num_modes_trunc):
		for j in range(0,num_modes_trunc):
			file.write('%4.1d %4.1d %18.10e %18.10e\n' % (i+1, j+1, noise_DMD[i,j].real, noise_DMD[i,j].imag) )
		file.write('\n')
	file.close()

#-----------------------------------------------------------------------
def write_eigenvalues(num_valid_DMD_modes,w,mode_index,energy_DMD,phase_DMD,dt,noise_DMD,DMD_trunc):
	cumulative_energy_DMD = np.array(np.zeros((num_valid_DMD_modes), dtype=np.float64))
	cumulative_energy_DMD[0] = energy_DMD[mode_index[0]]
	for i in range(1,num_valid_DMD_modes):
		cumulative_energy_DMD[i] = cumulative_energy_DMD[i-1] + energy_DMD.real[mode_index[i]]
	total_energy = cumulative_energy_DMD[num_valid_DMD_modes-1]

	print '\n   Writing eigenvalues to'
	filename = './results/DMD.eigenvalues.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, ev-stability real, ev-stability imag, phase, energy, cumulative, %energy, %cumulative, ev-DMD real, ev-DMD imag, log ev-DMD real, log ev-DMD imag, noise variance, M real (dmd space), M imag (dmd space)\n')
	file.write('#\n')
	for k in range(0,num_valid_DMD_modes):
		ev = w[mode_index[k]]
		ev_DMD = 1.0 + dt*(-1j*ev)
		theta = math.atan2(np.imag(ev_DMD),np.real(ev_DMD))
		log_ev_DMD = ( math.log(np.abs(ev_DMD)) + 1j*theta ) / dt
		noise_var = np.real(noise_DMD[mode_index[k],mode_index[k]])
		M_real = -np.real(DMD_trunc[mode_index[k],mode_index[k]])
		M_imag = -np.imag(DMD_trunc[mode_index[k],mode_index[k]])
		file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (k+1, np.real(ev), np.imag(ev), phase_DMD[mode_index[k]], energy_DMD.real[mode_index[k]], cumulative_energy_DMD[k], energy_DMD.real[mode_index[k]] / total_energy * 100.0, cumulative_energy_DMD[k] / total_energy * 100.0, np.real(ev_DMD), np.imag(ev_DMD), np.real(log_ev_DMD), np.imag(log_ev_DMD), noise_var, M_real, M_imag))
	file.close()

#-----------------------------------------------------------------------
def write_initial_conditions(num_valid_DMD_modes,mode_index,initial_conditions_DMD):
	print '\n   Writing initial conditions to'
	filename = './results/DMD.initial_conditions.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, amplitude: real, imag\n')
	file.write('#\n')
	for k in range(0,num_valid_DMD_modes):
		file.write('%4.1d %18.10e %18.10e\n' % (k+1, np.real(initial_conditions_DMD[mode_index[k]]), np.imag(initial_conditions_DMD[mode_index[k]]) ) )
	file.close()

#-----------------------------------------------------------------------
def write_temporal_modes_old(num_valid_DMD_modes,num_snapshots,dt,w,mode_index,energy_DMD,phase_DMD):
	print '\n   Writing temporal modes to'
	for k in range(0,num_valid_DMD_modes):
		k_index = k+1
		output_file_name = './results/DMD.temporal_mode_' + '%04d'%k_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'w')
		file.write('#\n')
		file.write('# time, real, imag\n')
		file.write('#\n')
		temporal_mode_var = 0.0
		for i in range(0,num_snapshots):
			amp_real = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.cos(-i*dt*np.real(w[mode_index[k]]))
			temporal_mode_var = temporal_mode_var + amp_real * amp_real
		temporal_mode_var = temporal_mode_var / num_snapshots
		amp_energy = np.abs(energy_DMD[mode_index[k]].real)
		for i in range(0,num_snapshots):
			amp_real = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.cos(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]]) * np.sqrt(amp_energy/temporal_mode_var)
			amp_imag = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.sin(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]]) * np.sqrt(amp_energy/temporal_mode_var)
                        file.write('%18.10e %18.10e %18.10e\n' % (i*dt, amp_real, amp_imag))
                file.close()

#-----------------------------------------------------------------------
def calculate_temporal_modes(num_valid_DMD_modes,num_snapshots,dt,w,mode_index,energy_DMD,phase_DMD,temporal_modes_DMD):
	for k in range(0,num_valid_DMD_modes):
		temporal_mode_var = 0.0
		for i in range(0,num_snapshots):
			amp_real = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.cos(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]])
			amp_imag = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.sin(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]])
			temporal_mode_var = temporal_mode_var + amp_real*amp_real + amp_imag*amp_imag
		temporal_mode_var = temporal_mode_var / num_snapshots
		amp_energy = np.abs(energy_DMD[mode_index[k]].real)
		for i in range(0,num_snapshots):
			amp_real = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.cos(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]]) * np.sqrt(amp_energy/temporal_mode_var)
			amp_imag = np.exp(i*dt*np.imag(w[mode_index[k]])) * np.sin(-i*dt*np.real(w[mode_index[k]]) + phase_DMD[mode_index[k]]) * np.sqrt(amp_energy/temporal_mode_var)
                        temporal_modes_DMD[i,mode_index[k]] = amp_real + 1j*amp_imag

#-----------------------------------------------------------------------
def write_temporal_modes(num_valid_DMD_modes,num_snapshots,dt,mode_index,temporal_modes_DMD):
	print '\n   Writing temporal modes to'
	for k in range(0,num_valid_DMD_modes):
		k_index = k+1
		output_file_name = './results/DMD.temporal_mode_' + '%04d'%k_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'w')
		file.write('#\n')
		file.write('# time, real, imag\n')
		file.write('#\n')
		for i in range(0,num_snapshots):
                        file.write('%18.10e %18.10e %18.10e\n' % (i*dt, temporal_modes_DMD[i,mode_index[k]].real, temporal_modes_DMD[i,mode_index[k]].imag))
                file.close()

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
