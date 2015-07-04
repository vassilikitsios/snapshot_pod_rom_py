# module io_vtk_m.py 

#-----------------------------------------------------------------------
# import libraries
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import sys

#-----------------------------------------------------------------------
def get_number_of_points(x_min,x_max,y_min,y_max,z_min,z_max,filename):
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.SetFileName(filename)
	reader.Update()
	if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
		extract.SetInput(reader.GetOutput())
		extract.MergingOn()
		extract.Update()
		num_points = extract.GetOutput().GetNumberOfPoints()
	else:
		num_points = reader.GetOutput().GetNumberOfPoints()
	print '   Number of points = ', num_points

	return num_points

#-----------------------------------------------------------------------
def get_cell_volumes(x_min,x_max,y_min,y_max,z_min,z_max,filename,cell_volume):
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.SetFileName(filename)
	reader.Update()
	if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
		extract.SetInput(reader.GetOutput())
		extract.MergingOn()
		extract.Update()
		cell_volume = VN.vtk_to_numpy(extract.GetOutput().GetPointData().GetScalars("cell_volume"))
	else:
		cell_volume = VN.vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars("cell_volume"))

#-----------------------------------------------------------------------
def get_grid_template(x_min,x_max,y_min,y_max,z_min,z_max,filename,grid):
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.SetFileName(filename)
	reader.Update()
	if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
		extract.SetInput(reader.GetOutput())
		extract.MergingOn()
		extract.Update()
		grid.DeepCopy(extract.GetOutput())
	else:
		grid.DeepCopy(reader.GetOutput())

#-----------------------------------------------------------------------
def read_data_and_build_snapshot_matrix(x_min,x_max,y_min,y_max,z_min,z_max,snapshot_dir,file_list,\
						num_snapshots,num_points,num_components,var_name, A):
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	extract_region = "false"
	if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
		extract_region = "true"
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
		extract.MergingOn()
	u_temp_read = np.array(np.zeros((num_points,3), dtype=np.float64))
	print '\n   Reading data ...'
	for j in range(0,num_snapshots+1):
		print '      Reading file ', file_list[j].strip(), 'file number ', j, ' of ', num_snapshots
		reader.SetFileName(snapshot_dir + file_list[j].strip())
		reader.Update()
		if (extract_region=="true"):
			extract.SetInput(reader.GetOutput())
			extract.Update()
			u_temp_read = VN.vtk_to_numpy(extract.GetOutput().GetPointData().GetVectors(var_name))
		else:
			u_temp_read = VN.vtk_to_numpy(reader.GetOutput().GetPointData().GetVectors(var_name))
		for k in range(0,num_components):
			A[k*num_points:(k+1)*num_points,j] = u_temp_read[:,k]

#-----------------------------------------------------------------------
def write_field(num_points,num_components,grid,field,filename):
	u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
	for k in range(0,num_components):
		u_temp[:,k] = field[k*num_points:(k+1)*num_points]
	u_vtk = VN.numpy_to_vtk(u_temp)
	u_vtk.SetName('u')
	grid.GetPointData().AddArray(u_vtk)
	writer = vtk.vtkUnstructuredGridWriter()
	writer.SetFileTypeToBinary()
	writer.SetInput(grid)
	writer.SetFileName(filename)
	writer.Write()

#-----------------------------------------------------------------------
def write_mean_field(num_points,num_components,grid,mean_field):
	print '   Writing mean field to'
	filename = './results/spatial_meanfield.vtk'
	print '      ', filename
	write_field(num_points,num_components,grid,mean_field,filename)

#-----------------------------------------------------------------------
def write_spatial_POD_modes(num_modes_to_write,num_points,num_components,grid,spatial_modes_trunc):
	print '\n   Writing spatial POD modes to'
	for j in range(0,num_modes_to_write):
		j_index = j+1
		filename = './results/POD.spatial_mode_' + '%04d'%j_index + '.vtk'
		print '      ', filename
		write_field(num_points,num_components,grid,spatial_modes_trunc[:,j].real,filename)

#-----------------------------------------------------------------------
def write_spatial_stochastic_modes(num_modes_to_write,num_points,num_components,grid,spatial_modes_trunc,noise_evector):

	print '   Writing spatial stochastic modes to'
	spatial_modes = np.array(np.zeros((num_components*num_points), dtype=np.float64))
	u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
	writer = vtk.vtkUnstructuredGridWriter()
	writer.SetFileTypeToBinary()

	for k in range(0,num_modes_to_write):
		k_index = k+1
		filename = './results/DMD.spatial_stochastic_mode_' + '%04d'%k_index + '.vtk'
		print '      ', filename

		spatial_modes = np.dot(spatial_modes_trunc, noise_evector[:,k])
		for i in range(0,num_components):
			u_temp[:,i] = spatial_modes[i*num_points:(i+1)*num_points]

		u_vtk = VN.numpy_to_vtk(u_temp)
		u_vtk.SetName('u')
		grid.GetPointData().AddArray(u_vtk)

		writer.SetFileName(filename)
		writer.SetInput(grid)
		writer.Write()

#-----------------------------------------------------------------------
#def write_spatial_DMD_modes(num_modes_to_write,num_points,num_components,grid,spatial_modes_trunc,spatial_modes_trunc_w,vr,vl,mode_index):
def write_spatial_DMD_modes(num_modes_to_write,num_points,num_components,grid,spatial_modes_trunc,vr,vl,mode_index):

	print '   Writing spatial DMD modes to'
	spatial_modes_DMD = np.array(np.zeros((num_components*num_points), dtype=np.complex64))
	spatial_modes_DMD_adjoint = np.array(np.zeros((num_components*num_points), dtype=np.complex64))
	u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
	u_temp_imag = np.array(np.zeros((num_points,3), dtype=np.float64))
	u_temp_adjoint_real = np.array(np.zeros((num_points,3), dtype=np.float64))
	u_temp_adjoint_imag = np.array(np.zeros((num_points,3), dtype=np.float64))
	writer = vtk.vtkUnstructuredGridWriter()
	writer.SetFileTypeToBinary()

	for k in range(0,num_modes_to_write):
		k_index = k+1
		filename = './results/DMD.spatial_deterministic_mode_' + '%04d'%k_index + '.vtk'
		print '      ', filename

		spatial_modes_DMD = np.dot(spatial_modes_trunc, vr[:,mode_index[k]])
       		#spatial_modes_DMD_adjoint = np.dot(spatial_modes_trunc_w, vl[:,mode_index[k]])
       		spatial_modes_DMD_adjoint = np.dot(spatial_modes_trunc, vl[:,mode_index[k]])
		for i in range(0,num_components):
			u_temp[:,i] = spatial_modes_DMD[i*num_points:(i+1)*num_points].real
			u_temp_imag[:,i] = spatial_modes_DMD[i*num_points:(i+1)*num_points].imag
			u_temp_adjoint_real[:,i] = spatial_modes_DMD_adjoint[i*num_points:(i+1)*num_points].real
			u_temp_adjoint_imag[:,i] = spatial_modes_DMD_adjoint[i*num_points:(i+1)*num_points].imag

		u_vtk = VN.numpy_to_vtk(u_temp)
		u_vtk.SetName('u')
		grid.GetPointData().AddArray(u_vtk)

		u_imag_vtk = VN.numpy_to_vtk(u_temp_imag)
		u_imag_vtk.SetName('u_imag')
		grid.GetPointData().AddArray(u_imag_vtk)

		u_adjoint_real_vtk = VN.numpy_to_vtk(u_temp_adjoint_real)
		u_adjoint_real_vtk.SetName('u_adjoint_real')
		grid.GetPointData().AddArray(u_adjoint_real_vtk)

		u_adjoint_imag_vtk = VN.numpy_to_vtk(u_temp_adjoint_imag)
		u_adjoint_imag_vtk.SetName('u_adjoint_imag')
		grid.GetPointData().AddArray(u_adjoint_imag_vtk)

		writer.SetFileName(filename)
		writer.SetInput(grid)
		writer.Write()

#		filename = './results/DMD.spatial_deterministic_mode_' + '%04d'%k_index + '.dat'
#                file = open(filename,'w')
#                file.write('#\n')
#                file.write('# i, j, u_r, u_i, v_r, v_i, w_r, w_i, uA_r, uA_i, vA_r, vA_i, wA_r, wA_i\n')
#                file.write('#\n')
#                for i in range(0,num_points):
#			point_num = i
#                        file.write('%18.10e %18.10e %18.10e %18.10e %18.10e \
#					%18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e \n' \
#					% (	u_temp[point_num,0],              u_temp_imag[point_num,0], \
#						u_temp[point_num,1],              u_temp_imag[point_num,1],\
#						u_temp[point_num,2],              u_temp_imag[point_num,2],\
#						u_temp_adjoint_real[point_num,0], u_temp_adjoint_imag[point_num,0],\
#						u_temp_adjoint_real[point_num,1], u_temp_adjoint_imag[point_num,1],\
#						u_temp_adjoint_real[point_num,2], u_temp_adjoint_imag[point_num,2] ))
#                file.close()

#-----------------------------------------------------------------------
# EOF
#-----------------------------------------------------------------------
