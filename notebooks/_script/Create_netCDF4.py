import glob
import scipy.io
import numpy as np
import netCDF4 as nc
#from __future__ import print_function

## ORCA2
### ORCA2_SRC

#print('ORCA2_SRC')

## Ba
#Ba_name=glob.glob('../_data/Exchange/Ba_ORCA2_SRC.mat')
#Ba_Mat=scipy.io.loadmat(Ba_name[0])
#Ba_ORCA2_2deg=Ba_Mat['Ba_ORCA2_2deg']
#nav_lon=Ba_Mat['nav_lon']
#nav_lat=Ba_Mat['nav_lat']
## d18O
#d18O_name=glob.glob('../_data/Exchange/d18O_ORCA2_SRC.mat')
#d18O_Mat=scipy.io.loadmat(d18O_name[0])
#d18O_ORCA2_2deg=d18O_Mat['d18O_ORCA2_2deg']
#
time_counter = [  15., 46., 75., 106., 136., 167., 197., 228., 259., 289., 320., 350.]
#
#nc_obj = nc.Dataset('../_data/Exchange/TRC_ORCA2_SRC.nc', 'w') # format='NETCDF4'
#nc_obj.description = 'River sources of Barium concentration in Arctic'
#print(nc_obj.file_format)
#
#time=nc_obj.createDimension('time', None)
#lat=nc_obj.createDimension('lat', 149)
#lon=nc_obj.createDimension('lon', 182)
#
#time_counter_obj=nc_obj.createVariable('time_counter', np.float32, ('time',), zlib=True)
#BaVar_obj=nc_obj.createVariable('Ba_ORCA2_2deg', np.float64, ('time', 'lat', 'lon'), zlib=True)
#d18OVar_obj=nc_obj.createVariable('d18O_ORCA2_2deg', np.float64, ('time', 'lat', 'lon'), zlib=True)
#nav_lat_obj=nc_obj.createVariable('nav_lat', np.float64, ('lat', 'lon'), zlib=True)
#nav_lon_obj=nc_obj.createVariable('nav_lon', np.float64, ('lat', 'lon'), zlib=True)
#
#time_counter_obj[:]=time_counter
#BaVar_obj[:]=Ba_ORCA2_2deg
#d18OVar_obj[:]=d18O_ORCA2_2deg
#nav_lat_obj[:]=nav_lat
#nav_lon_obj[:]=nav_lon
#
#nc_obj.close()
#
## ORCA2_Nomask
#
#print('ORCA2_Nomask')
#
#MAT = scipy.io.loadmat('../_data/Exchange/NEMO_ORCA2_Ba.mat')
#Ba_orca = MAT['Ba_ini_orca'][:]
#MAT = scipy.io.loadmat('../_data/Exchange/NEMO_ORCA2_d18O.mat')
#d18O_orca = MAT['d18O_ini_orca'][:]
##Ba_orca[np.isnan(Ba_orca)] = -999
#d18O_name=glob.glob('../_data/Exchange/NEMO_ORCA2_d18O_Grid.mat')
#d18O_Mat=scipy.io.loadmat(d18O_name[0])
#d18O_grid=d18O_Mat['d18O_grid_orca']
#
#boundary_name=glob.glob('../_data/Exchange/TRC_BOUND.mat')
#Boundary_Mat=scipy.io.loadmat(boundary_name[0])
#Ba_boundary=Boundary_Mat['Ba']
#d18O_boundary=Boundary_Mat['d18O']
#
#ini_obj = nc.Dataset('../_data/Exchange/TRC_ORCA2_Nomask.nc', 'w') # format='NETCDF4'
#ini_obj.description = 'Barium Initial Field'
#print(ini_obj.file_format)
#
#time=ini_obj.createDimension('time', None)
#lev = ini_obj.createDimension('lev', 31)
#lat=ini_obj.createDimension('lat', 149)
#lon=ini_obj.createDimension('lon', 182)
#
#time_counter_obj=ini_obj.createVariable('time_counter', np.float32, ('time',), zlib=True)
#BaVar_obj=ini_obj.createVariable('Ba', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
#d18OVar_obj=ini_obj.createVariable('d18O', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
#d18OGVar_obj=ini_obj.createVariable('d18O_grid', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
#BaBVar_obj=ini_obj.createVariable('Ba_boundary', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
#d18OBVar_obj=ini_obj.createVariable('d18O_boundary', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
#nav_lat_obj=ini_obj.createVariable('nav_lat', np.float64, ('lat', 'lon'), zlib=True)
#nav_lon_obj=ini_obj.createVariable('nav_lon', np.float64, ('lat', 'lon'), zlib=True)
#
##
#Ba = np.empty([1, 31, 149, 182])
#Ba[0, :, :, :] = Ba_orca
#
#BaB = np.empty([1, 31, 149, 182])
#BaB[0, :, :, :] = Ba_boundary
##
#d18O = np.empty([1, 31, 149, 182])
#d18O[0, :, :, :] = d18O_orca
##
#d18OB = np.empty([1, 31, 149, 182])
#d18OB[0, :, :, :] = d18O_boundary
##
#d18OG = np.empty([1, 31, 149, 182])
#d18OG[0, :, :, :] = d18O_grid
#
#time_counter_obj[:]=3600.0
#BaVar_obj[:]=Ba
#BaBVar_obj[:]=BaB
#d18OVar_obj[:]=d18O
#d18OBVar_obj[:]=d18OB
#d18OGVar_obj[:]=d18OG
#nav_lat_obj[:]=nav_lat
#nav_lon_obj[:]=nav_lon
#
#BaVar_obj.units='1e-6 mol/L'
#
#ini_obj.close()
#
## ANHA4_SRC
#
print('ANHA4_SRC')

Ba_name=glob.glob('../_data/Exchange/Ba_ANHA4_SRC.mat')
Ba_Mat=scipy.io.loadmat(Ba_name[0])
Ba_ANHA4=Ba_Mat['Ba_ANHA4']
nav_lon=Ba_Mat['nav_lon']
nav_lat=Ba_Mat['nav_lat']

d18O_name=glob.glob('../_data/Exchange/d18O_ANHA4_SRC.mat')
d18O_Mat=scipy.io.loadmat(d18O_name[0])
d18O_ANHA4=d18O_Mat['d18O_ANHA4']

boundary_name=glob.glob('../_data/Exchange/Ba_boundary_ANHA4.mat')
Boundary_Mat=scipy.io.loadmat(boundary_name[0])
Ba_boundary=Boundary_Mat['Ba_boundary']

nc_obj = nc.Dataset('../_data/Exchange/TRC_ANHA4_SRC.nc', 'w') # format='NETCDF4'
nc_obj.description = 'River sources of Barium concentration in Arctic'
print(nc_obj.file_format)

time=nc_obj.createDimension('time', None)
lat=nc_obj.createDimension('lon', 544)
lon=nc_obj.createDimension('lat', 800)

time_counter_obj=nc_obj.createVariable('time_counter', np.float32, ('time',), zlib=True)
BaVar_obj=nc_obj.createVariable('Ba_ANHA4', np.float64, ('time', 'lat', 'lon'), zlib=True)
d18OVar_obj=nc_obj.createVariable('d18O_ANHA4', np.float64, ('time', 'lat', 'lon'), zlib=True)
#BaB_obj = nc_obj.createVariable('Ba_boundary', np.float64, ('lat', 'lon'), zlib=True)
nav_lat_obj=nc_obj.createVariable('nav_lat', np.float64, ('lat', 'lon'), zlib=True)
nav_lon_obj=nc_obj.createVariable('nav_lon', np.float64, ('lat', 'lon'), zlib=True)

time_counter_obj[:]=time_counter
BaVar_obj[:]=np.transpose(Ba_ANHA4, [0, 2, 1])
d18OVar_obj[:]=np.transpose(d18O_ANHA4, [0, 2, 1])
#BaB_obj[:]=Ba_boundary.T
nav_lat_obj[:]=nav_lat.T
nav_lon_obj[:]=nav_lon.T

nc_obj.close()

## ANHA4_Nomask

print('ANHA4_Nomask')

#MAT = scipy.io.loadmat('../_data/Exchange/NEMO_ANHA4_Ba.mat')
#Ba_ANHA4 = MAT['Ba_ini_ANHA4'][:]
#MAT = scipy.io.loadmat('../_data/Exchange/NEMO_ANHA4_d18O_Grid.mat')
#d18O_ANHA4 = MAT['d18O_grid_ANHA'][:]

MAT = scipy.io.loadmat('../_data/Exchange/ANHA4_ini_latest.mat')
Ba_ANHA4 = MAT['Ba_ini'][:]
d18O_ANHA4 = MAT['d18O_ini'][:]
print(d18O_ANHA4.shape) # should be 3D

ini_obj = nc.Dataset('../_data/Exchange/TRC_ANHA4_Nomask.nc', 'w') # format='NETCDF4'
ini_obj.description = 'Barium Initial Field'
print(ini_obj.file_format)

time=ini_obj.createDimension('time', None)
lev = ini_obj.createDimension('lev', 50)
lat=ini_obj.createDimension('lat', 800)
lon=ini_obj.createDimension('lon', 544)

time_counter_obj=ini_obj.createVariable('time_counter', np.float32, ('time',), zlib=True)
BaVar_obj=ini_obj.createVariable('Ba', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
d18OVar_obj=ini_obj.createVariable('d18O', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
BaBVar_obj=ini_obj.createVariable('Ba_boundary', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
d18OBVar_obj=ini_obj.createVariable('d18O_boundary', np.float64, ('time', 'lev', 'lat', 'lon'), zlib=True)
nav_lat_obj=ini_obj.createVariable('nav_lat', np.float64, ('lat', 'lon'), zlib=True)
nav_lon_obj=ini_obj.createVariable('nav_lon', np.float64, ('lat', 'lon'), zlib=True)

# Initial field
Ba = Ba_ANHA4
d18O = d18O_ANHA4
# domain-wide sensitivity test
Ba_4D = np.empty([1, 50, 800, 544])
d18O_4D = np.empty([1, 50, 800, 544])
Ba_4D[0, :, :, :] = Ba
d18O_4D[0, :, :, :] = d18O
print(d18O_4D.shape) # should be 4D
# Open boundary settings
## Ba
Ba_boundary=np.transpose(Ba_boundary, [0, 2, 1])
BaB = np.empty([1, 50, 800, 544])
for i in range(50):
    BaB[0, i, Ba_boundary[i, :, :]>0.5] =  Ba[i, Ba_boundary[i, :, :]>0.5]
## d18O
d18OB = np.empty([1, 50, 800, 544])
d18OB[0, :, :, :] = d18O

#---------------------------------------------------------------------

time_counter_obj[:]=3600.0
BaVar_obj[:]=Ba_4D
d18OVar_obj[:]=d18O_4D
nav_lat_obj[:]=nav_lat
nav_lon_obj[:]=nav_lon
BaBVar_obj[:]=BaB
d18OBVar_obj[:]=d18OB

ini_obj.close()

print('ANHA4_PTT_mask')

MAT = scipy.io.loadmat('../_data/Exchange/Ba_PTT_LP.mat')
temp_mask = MAT['mask_rf'][:]
mask_rf = np.zeros([12, 800, 544])
for i in range(12):
    mask_rf[i, :, :] = temp_mask.T
print(mask_rf.shape) # should be 3D

ini_obj = nc.Dataset('../_data/Exchange/PTT_mask_LP.nc', 'w') # format='NETCDF4'
ini_obj.description = 'Mackenzie, Coppermine =1, otherwise = 1'
print(ini_obj.file_format)

time=ini_obj.createDimension('time', None)
lat=ini_obj.createDimension('lat', 800)
lon=ini_obj.createDimension('lon', 544)

time_counter_obj=ini_obj.createVariable('time_counter', np.float32, ('time',), zlib=True)
mask_obj=ini_obj.createVariable('mask_NA', np.float64, ('time', 'lat', 'lon'), zlib=True)
nav_lat_obj=ini_obj.createVariable('nav_lat', np.float64, ('lat', 'lon'), zlib=True)
nav_lon_obj=ini_obj.createVariable('nav_lon', np.float64, ('lat', 'lon'), zlib=True)

time_counter_obj[:]=time_counter
mask_obj[:]=mask_rf
nav_lat_obj[:]=nav_lat.T
nav_lon_obj[:]=nav_lon.T

ini_obj.close()

