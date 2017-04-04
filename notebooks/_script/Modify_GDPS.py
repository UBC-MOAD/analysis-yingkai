
import glob
import scipy.io
import numpy as np
import netCDF4 as nc

names = ['2004', '2008', '2012']

for nam in names:

    print(nam) 

    ori_obj = nc.Dataset('../_data/precip_gdps_y'+nam+'.nc', 'r')
    TIME = np.arange(1.0, 8785.0, 1.0)

    print('\tModify precip')

    precip_modified_1 = np.empty([1, 368, 801])
    precip_modified_1[:, :, :] = (ori_obj.variables['precip'][1439, :, :] + ori_obj.variables['precip'][1440, :, :])/2.0

    print('\tCreate new netCDF4')

    nc_obj = nc.Dataset('../_data/precip_gdps_y'+nam+'_new.nc', 'w') # format='NETCDF4'

    time = nc_obj.createDimension('time', None)
    lat  = nc_obj.createDimension('lat', 368)
    lon  = nc_obj.createDimension('lon', 801)

    time_obj = nc_obj.createVariable('TIME', np.float32, ('time',), zlib=True)
    lat_obj  = nc_obj.createVariable('LAT', np.float64, ('lat'), zlib=True) # 1d array
    lon_obj  = nc_obj.createVariable('LON', np.float64, ('lon'), zlib=True) # 1d array
    precip_obj  = nc_obj.createVariable('precip', np.float64, ('time', 'lat', 'lon'), zlib=True)

    time_obj[:]=TIME
    lon_obj[:] = ori_obj.variables['LON'][:]
    lat_obj[:] = ori_obj.variables['LAT'][:]

    print('\tWrite precip')

    precip_obj[0:1439, :, :] = ori_obj.variables['precip'][0:1439, :, :]
    for i in range(1439, 1464, 1):
    	precip_obj[i, :, :] = precip_modified_1
    precip_obj[1464:8784] = ori_obj.variables['precip'][1440:8760, :, :]

    nc_obj.close()

print('done')

