
import glob
import scipy.io
import numpy as np
import netCDF4 as nc

days_of_mon = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

ori_obj = nc.Dataset('../_data/precip_gdps_y2004.nc', 'r')

print('\tWork with daily files')

print('\t\tCreate new netCDF4')

nc_obj = nc.Dataset('daily.nc', 'w') # format='NETCDF4'

time = nc_obj.createDimension('time', None)
lat  = nc_obj.createDimension('lat', 368)
lon  = nc_obj.createDimension('lon', 801)

lat_obj  = nc_obj.createVariable('LAT', np.float64, ('lat'), zlib=True) # 1d array
lon_obj  = nc_obj.createVariable('LON', np.float64, ('lon'), zlib=True) # 1d array
precip_obj  = nc_obj.createVariable('precip', np.float64, ('time', 'lat', 'lon'), zlib=True)

print('\t\tWrite Lats/Lons')

lon_obj[:] = ori_obj.variables['LON'][:]
lat_obj[:] = ori_obj.variables['LAT'][:]

print('\t\tWrite precip')

for i in range(0, 365):
    precip_obj[i, :, :] = np.mean(ori_obj.variables['precip'][(i*24):(i+1)*24, :, :], 0)

nc_obj.close()

print('\tWork with monmean files')

daily_obj = nc.Dataset('daily.nc', 'r')

nc_obj = nc.Dataset('monmean.nc', 'w') # format='NETCDF4'

time = nc_obj.createDimension('time', None)
lat  = nc_obj.createDimension('lat', 368)
lon  = nc_obj.createDimension('lon', 801)

lat_obj  = nc_obj.createVariable('LAT', np.float64, ('lat'), zlib=True) # 1d array
lon_obj  = nc_obj.createVariable('LON', np.float64, ('lon'), zlib=True) # 1d array
precip_obj  = nc_obj.createVariable('precip', np.float64, ('time', 'lat', 'lon'), zlib=True)

print('\t\tWrite Lats/Lons')

lon_obj[:] = ori_obj.variables['LON'][:]
lat_obj[:] = ori_obj.variables['LAT'][:]

print('\t\tWrite precip')

for i in range(0, 12):
    temp = np.zeros([368, 801])
    for j in range(0, days_of_mon[i]):
        temp += daily_obj.variables['precip'][j, :, :]
        temp = temp/days_of_mon[i]
        precip_obj[i, :, :] = temp

nc_obj.close()



print('\tdone')

