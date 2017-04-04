import glob
import scipy.io
import numpy as np
import netCDF4 as nc


def seasonal_decomp3d(data, method=0):
    '''
    =======================================================================
    Remove the seasonal cycle from 1D data
                            ----- created on 2015/06/15, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        data = seasonal_decomp(...)
    -----------------------------------------------------------------------
    Input:
            data: Time should be the first dim.
            method: removal done by anomaly (=0) or normalize (=1)
    ======================================================================= 
    '''
    for mon in range(12):
        temp_data = np.nanmean(data[mon:len(data):12, :, :], 0)
        if method == 0:
            data[mon:len(data):12, :, :] = data[mon:len(data):12, :, :]-temp_data
    return data

# ========================================================================================== #

names_u = sorted(glob.glob('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/GDPS/u10*monmean.nc'))
names_v = sorted(glob.glob('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/GDPS/v10*monmean.nc'))
names_slp = sorted(glob.glob('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/GDPS/slp*monmean.nc'))

L = len(names_u)

ua = np.empty([L*12, 368, 801])
va = np.empty([L*12, 368, 801])
slpa = np.empty([L*12, 368, 801])

print('Calculate the anomaly of GDPS atmos files')

for i in range(L):
 
    print('\tPatching file: {}'.format(i))
    
    u_obj = nc.Dataset(names_u[i])
    v_obj = nc.Dataset(names_v[i])
    slp_obj = nc.Dataset(names_slp[i])
    
    ua[i*12:(i+1)*12, :, :] = u_obj.variables['u_wind'][:]
    va[i*12:(i+1)*12, :, :] = v_obj.variables['v_wind'][:]
    slpa[i*12:(i+1)*12, :, :] = slp_obj.variables['atmpres'][:]
    
ua = seasonal_decomp3d(ua)
va = seasonal_decomp3d(va)
slpa = seasonal_decomp3d(slpa)

# ======================================================================================== #

print('Calculate the anomaly of ANHA4-EXH001 files')                             

obj_z = nc.Dataset('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/vozocrtx_monmean.nc')
obj_m = nc.Dataset('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/vomecrty_monmean.nc')
obj_h = nc.Dataset('/ocean/yingkai/GEOTRACES/FORCING/ANHA4/sossheig_monmean.nc')

print('\tvozocrtx')

zc = obj_z.variables['vozocrtx'][:]
zc[zc>100] = np.nan
zca = seasonal_decomp3d(zc)

print('\tvomecrty')

mc = obj_m.variables['vomecrty'][:]
mc[mc>100] = np.nan
mca = seasonal_decomp3d(mc)

print('\tsossheig')

ssh = obj_h.variables['sossheig'][:]
ssh[ssh>100] = np.nan
ssha = seasonal_decomp3d(ssh)

print('Save')
                             
save_var = {'u10':ua, 'v10':va, 'slp':slpa, 'zonalc':zca, 'meridc':mca, 'ssh':ssha}
scipy.io.savemat('../_data/Exchange/anomaly_for_AO.mat', mdict=save_var)
