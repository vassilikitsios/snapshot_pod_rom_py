# module rom_coefficients_m.py

#-----------------------------------------------------------------------
# import libraries
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import string
import math
import scipy.linalg as linalg

#-----------------------------------------------------------------------
def get_number_of_points_and_boundary_faces(pod_mode_dir,calibrate_coefficients):
	filename = pod_mode_dir+"/spatial_meanfield.vtk"
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.ReadAllTensorsOn()
	reader.SetFileName(filename)
	reader.Update()
	num_points = reader.GetOutput().GetNumberOfPoints()

	if( (calibrate_coefficients=="none") or (calibrate_coefficients=="constant") ):
		geom_filter = vtk.vtkDataSetSurfaceFilter()
		geom_filter.SetInput(reader.GetOutput())
		geom_filter.Update()
		boundary_faces = vtk.vtkPolyData()
		boundary_faces = geom_filter.GetOutput()

		point_to_cell_data = vtk.vtkPointDataToCellData()
		point_to_cell_data.SetInput(geom_filter.GetOutput())
		point_to_cell_data.Update()
		num_boundary_faces = point_to_cell_data.GetOutput().GetNumberOfCells()
	else:
		num_boundary_faces = 1

	return num_points, num_boundary_faces

#-----------------------------------------------------------------------
def read_data_and_differentiate(pod_mode_dir,num_modes,num_points,num_boundary_faces,num_dimensions,num_components,correct_for_cell_volumes, \
				u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,\
				uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF,cell_volume,norm,calibrate_coefficients):

	u_read = np.array(np.zeros((num_points,3), dtype=np.float64))
	dudx_read = np.array(np.zeros((num_points,3,3), dtype=np.float64))
	uF_read = np.array(np.zeros((num_boundary_faces,3), dtype=np.float64))
	dudxF_read = np.array(np.zeros((num_boundary_faces,3,3), dtype=np.float64))

	# ..................................................
	# Read meanfield and spatially differentiate
	filename = pod_mode_dir + "/spatial_meanfield.vtk"
	print '   Reading file ', filename.strip()
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.ReadAllTensorsOn()
	reader.SetFileName(filename)
	reader.Update()
	grid = vtk.vtkUnstructuredGrid()
	grid.DeepCopy(reader.GetOutput())

	# get cell volumes
	if (correct_for_cell_volumes=="true"):
	        cell_volume[:] = VN.vtk_to_numpy(grid.GetPointData().GetScalars("cell_volume"))
	else:
        	cell_volume[:] = 1.0

	# get mean velocity field within volume
	u_read = VN.vtk_to_numpy(grid.GetPointData().GetVectors("u"))
	u[:,0] = u_read[:,0].copy()
	v[:,0] = u_read[:,1].copy()
	w[:,0] = u_read[:,2].copy()
	grid.GetPointData().SetVectors(grid.GetPointData().GetVectors("u"))

	# get mean velocity derivative field within volume
	differentiator = vtk.vtkCellDerivatives()
	differentiator.SetTensorModeToComputeGradient()
	differentiator.SetInput(grid)
	differentiator.Update()
	cell_to_point_data = vtk.vtkCellDataToPointData()
	cell_to_point_data.SetInput(differentiator.GetOutput())
	cell_to_point_data.Update()
        grid.GetPointData().SetTensors(cell_to_point_data.GetOutput().GetPointData().GetTensors())
	dudx_read = VN.vtk_to_numpy(cell_to_point_data.GetOutput().GetPointData().GetTensors())
	dudx[:,0] = dudx_read[:,0].copy() 
	dvdx[:,0] = dudx_read[:,1].copy() 
	dwdx[:,0] = dudx_read[:,2].copy() 
	dudy[:,0] = dudx_read[:,3].copy() 
	dvdy[:,0] = dudx_read[:,4].copy() 
	dwdy[:,0] = dudx_read[:,5].copy() 
	dudz[:,0] = dudx_read[:,6].copy() 
	dvdz[:,0] = dudx_read[:,7].copy() 
	dwdz[:,0] = dudx_read[:,8].copy() 

	# extract boundary surface data
	if( (calibrate_coefficients=="none") or (calibrate_coefficients=="constant") ):
		geom_filter = vtk.vtkDataSetSurfaceFilter()
		geom_filter.SetInput(grid)
		geom_filter.Update()
		boundary_faces = vtk.vtkPolyData()
		boundary_faces = geom_filter.GetOutput()
		point_to_cell_data = vtk.vtkPointDataToCellData()
		point_to_cell_data.SetInput(geom_filter.GetOutput())
		point_to_cell_data.Update()

		# get mean velocity field on boundary surface
		uF_read = VN.vtk_to_numpy(point_to_cell_data.GetOutput().GetCellData().GetVectors("u"))
		uF[:,0] = uF_read[:,0]
		vF[:,0] = uF_read[:,1]
		wF[:,0] = uF_read[:,2]

		# get mean derivative velocity field on boundary surface
		dudxF_read = VN.vtk_to_numpy(point_to_cell_data.GetOutput().GetCellData().GetTensors())
		dudxF[:,0] = dudxF_read[:,0] 
		dvdxF[:,0] = dudxF_read[:,1] 
		dwdxF[:,0] = dudxF_read[:,2] 
		dudyF[:,0] = dudxF_read[:,3] 
		dvdyF[:,0] = dudxF_read[:,4] 
		dwdyF[:,0] = dudxF_read[:,5] 
		dudzF[:,0] = dudxF_read[:,6] 
		dvdzF[:,0] = dudxF_read[:,7] 
		dwdzF[:,0] = dudxF_read[:,8] 

		# get boundary face normals
		norm_filter = vtk.vtkPolyDataNormals()
		norm_filter.ComputeCellNormalsOn()
		norm_filter.SetInput(boundary_faces)
		norm_filter.Update()
		area_filter = vtk.vtkMeshQuality()
		area_filter.SetQuadQualityMeasureToArea()
		area_filter.SetTriangleQualityMeasureToArea()
		area_filter.SetInput(boundary_faces)
		area_filter.Update()
		for j in range(0,num_boundary_faces):
			area = area_filter.GetOutput().GetCellData().GetArray("Quality").GetComponent(j,0)
			norm[j,:] = norm_filter.GetOutput().GetCellData().GetNormals().GetTuple(j)
			norm[j,:] = norm[j,:] * area

	# ..................................................
	# Read modes and spatially differentiate
	for j in range(0,num_modes):
		j_index = j+1
		filename = pod_mode_dir + '/POD.spatial_mode_' + '%04d'%j_index + '.vtk'
		print '   Reading file ', filename.strip(), 'file number ', j_index, ' of ', num_modes
		reader.SetFileName(filename)
		reader.Update()

		# get mode velocity field within volume
		u_read = VN.vtk_to_numpy(reader.GetOutput().GetPointData().GetVectors("u"))
		u[:,j_index] = u_read[:,0]
		v[:,j_index] = u_read[:,1]
		w[:,j_index] = u_read[:,2]

		# get mode velocity derivative fields within volume
		grid.GetPointData().SetVectors(reader.GetOutput().GetPointData().GetVectors("u"))
		differentiator.SetInput(grid)
		differentiator.Update()
		cell_to_point_data.SetInput(differentiator.GetOutput())
		cell_to_point_data.Update()
        	grid.GetPointData().SetTensors(cell_to_point_data.GetOutput().GetPointData().GetTensors())
		dudx_read = VN.vtk_to_numpy(cell_to_point_data.GetOutput().GetPointData().GetTensors())
		dudx[:,j_index] = dudx_read[:,0] 
		dvdx[:,j_index] = dudx_read[:,1] 
		dwdx[:,j_index] = dudx_read[:,2] 
		dudy[:,j_index] = dudx_read[:,3] 
		dvdy[:,j_index] = dudx_read[:,4] 
		dwdy[:,j_index] = dudx_read[:,5] 
		dudz[:,j_index] = dudx_read[:,6] 
		dvdz[:,j_index] = dudx_read[:,7] 
		dwdz[:,j_index] = dudx_read[:,8]

		# extract boundary surface data
		if( (calibrate_coefficients=="none") or (calibrate_coefficients=="constant") ):
			geom_filter.SetInput(grid)
			geom_filter.Update()
			boundary_faces = geom_filter.GetOutput()
			point_to_cell_data.SetInput(geom_filter.GetOutput())
			point_to_cell_data.Update()

			# get mode velocity field on boundary surface
			uF_read = VN.vtk_to_numpy(point_to_cell_data.GetOutput().GetCellData().GetVectors("u"))
			uF[:,j_index] = uF_read[:,0]
			vF[:,j_index] = uF_read[:,1]
			wF[:,j_index] = uF_read[:,2]

			# get mean derivative velocity field on boundary surface
			dudxF_read = VN.vtk_to_numpy(point_to_cell_data.GetOutput().GetCellData().GetTensors())
			dudxF[:,j_index] = dudxF_read[:,0] 
			dvdxF[:,j_index] = dudxF_read[:,1] 
			dwdxF[:,j_index] = dudxF_read[:,2] 
			dudyF[:,j_index] = dudxF_read[:,3] 
			dvdyF[:,j_index] = dudxF_read[:,4] 
			dwdyF[:,j_index] = dudxF_read[:,5] 
			dudzF[:,j_index] = dudxF_read[:,6] 
			dvdzF[:,j_index] = dudxF_read[:,7] 
			dwdzF[:,j_index] = dudxF_read[:,8] 

	# ..................................................
	# zero appropriate coefficients
	if (num_dimensions<3):
		dudz[:,:] = 0.0
		dvdz[:,:] = 0.0
		dwdz[:,:] = 0.0
	if (num_dimensions<2):
		dudy[:,:] = 0.0
		dvdy[:,:] = 0.0
		dwdy[:,:] = 0.0
	if (num_components<3):
		w[:,:] = 0.0
		dwdx[:,:] = 0.0
		dwdy[:,:] = 0.0
		dwdz[:,:] = 0.0
	if (num_components<2):
		v[:,:] = 0.0
		dvdx[:,:] = 0.0
		dvdy[:,:] = 0.0
		dvdz[:,:] = 0.0
	if( (calibrate_coefficients=="none") or (calibrate_coefficients=="constant") ):
		if (num_dimensions<3):
			dudzF[:,:] = 0.0
			dvdzF[:,:] = 0.0
			dwdzF[:,:] = 0.0
		if (num_dimensions<2):
			dudyF[:,:] = 0.0
			dvdyF[:,:] = 0.0
			dwdyF[:,:] = 0.0
		if (num_components<3):
			wF[:,j_index] = 0.0
			dwdxF[:,:] = 0.0
			dwdyF[:,:] = 0.0
			dwdzF[:,:] = 0.0
		if (num_components<2):
			vF[:,:] = 0.0
			dvdxF[:,:] = 0.0
			dvdyF[:,:] = 0.0
			dvdzF[:,:] = 0.0

#-----------------------------------------------------------------------
def calculate_constant_coefficients(num_modes, \
				u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
				uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
				Re,cell_volume,norm,constant):
	print "   constant"
	# Equation 5.13 of Kitsios (2012)
	for i in range(0,num_modes):
		print "      mode: ", i+1
		constant[i] = - sum(cell_volume[:]*((u[:,0]*dudx[:,0] + v[:,0]*dudy[:,0] + w[:,0]*dudz[:,0])*u[:,i+1] \
						  + (u[:,0]*dvdx[:,0] + v[:,0]*dvdy[:,0] + w[:,0]*dvdz[:,0])*v[:,i+1] \
						  + (u[:,0]*dwdx[:,0] + v[:,0]*dwdy[:,0] + w[:,0]*dwdz[:,0])*w[:,i+1])) \
			      - sum(cell_volume[:]/Re*(dudx[:,i+1]*dudx[:,0] + dudy[:,i+1]*dudy[:,0] + dudz[:,i+1]*dudz[:,0] \
					 	      +dvdx[:,i+1]*dvdx[:,0] + dvdy[:,i+1]*dvdy[:,0] + dvdz[:,i+1]*dvdz[:,0] \
					 	      +dwdx[:,i+1]*dwdx[:,0] + dwdy[:,i+1]*dwdy[:,0] + dwdz[:,i+1]*dwdz[:,0])) \
		 	      - sum(1.0/Re*((dudxF[:,0]*norm[:,0] + dvdxF[:,0]*norm[:,1] + dwdxF[:,0]*norm[:,2])*uF[:,i+1] \
		                 	  + (dudyF[:,0]*norm[:,0] + dvdyF[:,0]*norm[:,1] + dwdyF[:,0]*norm[:,2])*vF[:,i+1] \
	               			  + (dudzF[:,0]*norm[:,0] + dvdzF[:,0]*norm[:,1] + dwdzF[:,0]*norm[:,2])*wF[:,i+1]))
	return

#-----------------------------------------------------------------------
def calculate_linear_coefficients(num_modes, \
				u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, \
				uF,vF,wF,dudxF,dudyF,dudzF,dvdxF,dvdyF,dvdzF,dwdxF,dwdyF,dwdzF, \
				Re,cell_volume,norm,linear):
	print "   linear"
	# Equation 5.14 of Kitsios (2012)
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			print "      mode: ", i+1, j+1
			linear[i,j] = - sum(cell_volume[:]*((u[:,0]*dudx[:,j+1] + v[:,0]*dudy[:,j+1] + w[:,0]*dudz[:,j+1] )*u[:,i+1] \
							 +  (u[:,0]*dvdx[:,j+1] + v[:,0]*dvdy[:,j+1] + w[:,0]*dvdz[:,j+1] )*v[:,i+1] \
							 +  (u[:,0]*dwdx[:,j+1] + v[:,0]*dwdy[:,j+1] + w[:,0]*dwdz[:,j+1] )*w[:,i+1])) \
				      - sum(cell_volume[:]*((u[:,j+1]*dudx[:,0] + v[:,j+1]*dudy[:,0] + w[:,j+1]*dudz[:,0] )*u[:,i+1] \
							  + (u[:,j+1]*dvdx[:,0] + v[:,j+1]*dvdy[:,0] + w[:,j+1]*dvdz[:,0] )*v[:,i+1] \
							  + (u[:,j+1]*dwdx[:,0] + v[:,j+1]*dwdy[:,0] + w[:,j+1]*dwdz[:,0] )*w[:,i+1])) \
				      - sum(cell_volume[:]/Re*( dudx[:,i+1]*dudx[:,j+1] + dudy[:,i+1]*dudy[:,j+1] + dudz[:,i+1]*dudz[:,j+1] \
							     +  dvdx[:,i+1]*dvdx[:,j+1] + dvdy[:,i+1]*dvdy[:,j+1] + dvdz[:,i+1]*dvdz[:,j+1] \
							     +  dwdx[:,i+1]*dwdx[:,j+1] + dwdy[:,i+1]*dwdy[:,j+1] + dwdz[:,i+1]*dwdz[:,j+1])) \
				      - sum(1.0/Re*((dudxF[:,j+1]*norm[:,0] + dvdxF[:,j+1]*norm[:,1] + dwdxF[:,j+1]*norm[:,2])*uF[:,i+1] \
	                	  		  + (dudyF[:,j+1]*norm[:,0] + dvdyF[:,j+1]*norm[:,1] + dwdyF[:,j+1]*norm[:,2])*vF[:,i+1] \
	     					  + (dudzF[:,j+1]*norm[:,0] + dvdzF[:,j+1]*norm[:,1] + dwdzF[:,j+1]*norm[:,2])*wF[:,i+1]))
	return

#-----------------------------------------------------------------------
def calculate_quadratic_coefficients(num_modes,u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,cell_volume,quadratic):
	print "   quadratic"
	# Equation 5.15 of Kitsios (2012)
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				print "      mode: ", i+1, j+1, k+1
				quadratic[i,j,k] = -sum(cell_volume[:]*((u[:,j+1]*dudx[:,k+1]+v[:,j+1]*dudy[:,k+1]+w[:,j+1]*dudz[:,k+1])*u[:,i+1] \
								      + (u[:,j+1]*dvdx[:,k+1]+v[:,j+1]*dvdy[:,k+1]+w[:,j+1]*dvdz[:,k+1])*v[:,i+1] \
								      + (u[:,j+1]*dwdx[:,k+1]+v[:,j+1]*dwdy[:,k+1]+w[:,j+1]*dwdz[:,k+1])*w[:,i+1]))
	return

#-----------------------------------------------------------------------
def write_constant_coefficients(num_modes,constant):
	filename = './results/ROM.constant_coefficients.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode i, coefficient\n')
	file.write('#\n')
	for i in range(0,num_modes):
		file.write('%4.1d %18.10e\n' % (i+1, constant[i]) )
	file.close()

#-----------------------------------------------------------------------
def write_linear_coefficients(num_modes,linear):
	filename = './results/ROM.linear_coefficients.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode i, mode j, coefficient\n')
	file.write('#\n')
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			file.write('%4.1d %4.1d %18.10e\n' % (i+1, j+1, linear[i,j]) )
		file.write('\n')
	file.close()

#-----------------------------------------------------------------------
def write_quadratic_coefficients(num_modes,quadratic):
	filename = './results/ROM.quadratic_coefficients.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode i, mode j, mode k, coefficient\n')
	file.write('#\n')
	for i in range(0,num_modes):
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				file.write('%4.1d %4.1d %4.1d %18.10e\n' % (i+1, j+1, k+1, quadratic[i,j,k]) )
			file.write('\n')
	file.close()

#-----------------------------------------------------------------------
def read_temporal_modes(pod_mode_dir,num_modes,num_snapshots,a,t):
	print '\n   Reading temporal modes ...'
	for j in range(0,num_modes):
		j_index = j+1
		filename = pod_mode_dir+'/POD.temporal_mode_' + '%04d'%j_index + '.dat'
		file = open(filename,'r')
		for i in range(0,3):
			linestring = file.readline()
		for i in range(0,num_snapshots):
			linestring = file.readline()
			linelist = string.split(linestring)
			a[j,i]= float(linelist[1])
			if (j==0):
				t[i]= float(linelist[0])
		file.close()

#-----------------------------------------------------------------------
def interpolate_temporal_modes(num_modes,num_snapshots,interp_factor,t,a,t_interp,a_interp,da_dt_interp):
	A = np.array(np.zeros((num_snapshots), dtype=np.float64))
	B = np.array(np.zeros((num_snapshots-1), dtype=np.float64))
	C = np.array(np.zeros((num_snapshots), dtype=np.float64))
	D = np.array(np.zeros((num_snapshots-1), dtype=np.float64))

	Q = np.array(np.zeros((num_snapshots), dtype=np.float64))
	diag1 = np.array(np.zeros((num_snapshots-1), dtype=np.float64))
	diag2 = np.array(np.zeros((num_snapshots), dtype=np.float64))
	diag3 = np.array(np.zeros((num_snapshots-1), dtype=np.float64))

	dt = t[1]-t[0]
	for i in range(0,(num_snapshots-1)*interp_factor+1):
		t_interp[i] = dt*i/interp_factor

	for j in range(0,num_modes):
		# assign A coefficients
		A[:] = a[j,:]

		# solve for C coefficients
		diag1[:] = dt
		diag1[num_snapshots-2] = 0.0
		diag2[:] = 4.0*dt
		diag2[0] = 1.0
		diag2[num_snapshots-1] = 1.0
		diag3[:] = dt
		diag3[0] = 0.0
		Q[0] = 0.0
		for i in range(1,num_snapshots-1):
			Q[i] = 3.0/dt*(A[i+1] - 2.0*A[i] + A[i-1])
		Q[num_snapshots-1] = 0.0

		for i in range(1,num_snapshots-1):
			diag2[i] = diag2[i] - diag3[i-1] * diag1[i-1] / diag2[i-1]
			Q[i] = Q[i] - Q[i-1] * diag1[i-1] / diag2[i-1]

	        C[num_snapshots-1] = 0.0
		for i in range(num_snapshots-2,0,-1):
			C[i] = (Q[i] - diag3[i] * C[i+1]) / diag2[i]

		# evaluate D coefficients
		for i in range(0,num_snapshots-1):
			D[i] = (C[i+1] - C[i]) / 3.0 / dt

		# evaluate B coefficients
		for i in range(0,num_snapshots-1):
			B[i] = (A[i+1] - A[i]) / dt - dt/3.0*(2.0*C[i] + C[i+1])
			#B[i] = (A[i+1] - A[i]) / dt - C[i]*dt - D[i]*dt*dt

		# check continuity of the function and the first and second derivatives
		tol = 1.0e-15
		#tol = 0.0
		for i in range(0,num_snapshots-2):
			errorA = A[i] + B[i]*dt + C[i]*dt*dt + D[i]*dt*dt*dt - A[i+1]
			if (errorA*errorA > tol):
				print "A", A[i] + B[i]*dt + C[i]*dt*dt + D[i]*dt*dt*dt, A[i+1], i
			errorB = B[i] + 2.0*C[i]*dt + 3.0*D[i]*dt*dt - B[i+1]
			if (errorB*errorB > tol):
				print "B", B[i] + 2.0*C[i]*dt + 3.0*D[i]*dt*dt, B[i+1], i
			errorC = C[i] + 3.0*D[i]*dt - C[i+1]
			if (errorC*errorC > tol):
				print "C", C[i] + 3.0*D[i]*dt, C[i+1], i

		# interpolate cubic spline
		for i in range(0,num_snapshots-1):
			for k in range(0,interp_factor):
				pos = i*interp_factor+k
				dt_i = t_interp[pos] - t[i]
				a_interp[j,pos] = A[i] + B[i]*dt_i + C[i]*dt_i*dt_i + D[i]*dt_i*dt_i*dt_i
				da_dt_interp[j,pos] = B[i] + 2.0*C[i]*dt_i + 3.0*D[i]*dt_i*dt_i
		a_interp[j,pos+1] = A[num_snapshots-1]
		da_dt_interp[j,pos+1] = B[num_snapshots-2] + 2.0*C[num_snapshots-2]*dt_i + 3.0*D[num_snapshots-2]*dt_i*dt_i

#-----------------------------------------------------------------------
def write_temporal_modes(num_valid_modes,num_snapshots,dt,a,da_dt):
	print '\n   Writing temporal modes to'
	for j in range(0,num_valid_modes):
		j_index = j+1
		output_file_name = './results/POD_interp.temporal_mode_' + '%04d'%j_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'w')
		file.write('#\n')
		file.write('# time, a, da_dt\n')
		file.write('#\n')
		for i in range(0,num_snapshots):
			file.write('%18.10e %18.10e %18.10e\n' % (i*dt, a[j,i], da_dt[j,i]) )
		file.close()

#-----------------------------------------------------------------------
def calibrate_constant_and_linear_coefficients(num_modes,num_timesteps,a,da_dt,quadratic,constant,linear):
	# Equation 5.21 of Kitsios (2012)
	H = np.array(np.zeros((num_modes), dtype=np.float64))
	for m in range(0,num_modes):
		for t in range(0,num_timesteps):
			H[m] = H[m] + a[m,t]

	# Equation 5.22 of Kitsios (2012)
	R = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for j in range(0,num_modes):
		for m in range(0,num_modes):
			for t in range(0,num_timesteps):
				R[j,m] = R[j,m] + a[j,t] * a[m,t]

	# Equation 5.27 of Kitsios (2012)
	P = np.array(np.zeros((num_modes,num_modes,num_modes), dtype=np.float64))
	for j in range(0,num_modes):
		for k in range(0,num_modes):
			for m in range(0,num_modes):
				for t in range(0,num_timesteps):
					P[j,k,m] = P[j,k,m] + a[j,t] * a[k,t] * a[m,t]

	# Equation 5.28 of Kitsios (2012)
	T = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for m in range(0,num_modes):
			for t in range(0,num_timesteps):
				T[l,m] = T[l,m] + da_dt[l,t] * a[m,t]

	# Equation 5.23 of Kitsios (2012)
	S = np.array(np.zeros((num_modes), dtype=np.float64))
	for m in range(0,num_modes):
		for t in range(0,num_timesteps):
			S[m] = S[m] + da_dt[m,t]

	# Component of right-hand-side of Equation 5.26 of Kitsios (2012)
	QZ = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for m in range(0,num_modes):
			for j in range(0,num_modes):
				for k in range(0,num_modes):
					QZ[l,m] = QZ[l,m] + quadratic[l,j,k] * (P[j,k,m] - R[j,k] * H[m] / num_timesteps)

	# Component of right-hand-side of Equation 5.26 of Kitsios (2012)
	SH = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for m in range(0,num_modes):
			SH[l,m] = S[l] * H[m] / num_timesteps
	
	# Component of right-hand-side of Equation 5.25 of Kitsios (2012)
	HH = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for m in range(0,num_modes):
			HH[l,m] = H[l] * H[m] / num_timesteps

	# Component of right-hand-side of Equation 5.20 of Kitsios (2012)
	QR = np.array(np.zeros((num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				QR[l] = QR[l] + quadratic[l,j,k] * R[j,k]

	# Solution of Equation 5.24 of Kitsios (2012)
	linear[:,:] = np.dot(T-SH-QZ, linalg.inv(R-HH))

	# Equation 5.20 of Kitsios (2012)
	constant[:] = (S - np.dot(linear,H) - QR) / num_timesteps

#-----------------------------------------------------------------------
def calibrate_constant_coefficients(num_modes,num_timesteps,a,da_dt,quadratic,constant,linear):
	# Equation 5.21 of Kitsios (2012)
	H = np.array(np.zeros((num_modes), dtype=np.float64))
	for m in range(0,num_modes):
		for t in range(0,num_timesteps):
			H[m] = H[m] + a[m,t]

	# Equation 5.22 of Kitsios (2012)
	R = np.array(np.zeros((num_modes,num_modes), dtype=np.float64))
	for j in range(0,num_modes):
		for m in range(0,num_modes):
			for t in range(0,num_timesteps):
				R[j,m] = R[j,m] + a[j,t] * a[m,t]

	# Equation 5.23 of Kitsios (2012)
	S = np.array(np.zeros((num_modes), dtype=np.float64))
	for m in range(0,num_modes):
		for t in range(0,num_timesteps):
			S[m] = S[m] + da_dt[m,t]

	# Component of right-hand-side of Equation 5.20 of Kitsios (2012)
	QR = np.array(np.zeros((num_modes), dtype=np.float64))
	for l in range(0,num_modes):
		for j in range(0,num_modes):
			for k in range(0,num_modes):
				QR[l] = QR[l] + quadratic[l,j,k] * R[j,k]

	# Calibration of only the constant coefficients and not the linear coefficients.
	constant[:] = (S - np.dot(linear,H) - QR) / num_timesteps

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
