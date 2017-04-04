import sys
sys.path.insert(0, '../_libs/')
import NEMO_tools as Nts

import glob
import scipy.io
import numpy as np
import netCDF4 as nc
import NEMO_tools as Nts

def interp_xy(x0, y0, data, x1, y1):
    data_int1 = Nts.reporj_NEMOgrid(x0, y0, data, x1, y1, method='linear')
    data_int2 = Nts.reporj_NEMOgrid(x0, y0, data, x1, y1, method='nearest')
    id1 = np.logical_and(np.isnan(data_int1), np.abs(x1)>90)
    id2 = np.logical_and(np.isnan(data_int1), (y1)>61)
    data_int1[id1]=data_int2[id1]
    data_int1[id2]=data_int2[id2]
    return data_int1

x = np.linspace(-180, 180, 180)
y = np.linspace(60, 90, 60)
lon, lat = np.meshgrid(x, y)

Ba_FName = '/ocean/yingkai/GEOTRACES/Simulations/BARIUM01_1m_20020101_20140103_ptrc_T.nc'
U_FName = '/ocean/yingkai/GEOTRACES/FORCING/ANHA4/vozocrtx_monmean.nc'
V_FName = '/ocean/yingkai/GEOTRACES/FORCING/ANHA4/vomecrty_monmean.nc'

ptrc_obj = nc.Dataset(Ba_FName)
nav_lat = ptrc_obj.variables['nav_lat'][:]
nav_lon = ptrc_obj.variables['nav_lon'][:]

u_obj = nc.Dataset(U_FName)
v_obj = nc.Dataset(V_FName)

Ba = np.zeros([144, 50, 60, 180])
u = np.zeros([144, 60, 180])
v = np.zeros([144, 60, 180])

for i in range(144):
    print('Time {}'.format(i))
    u[i, :, :]  = interp_xy(nav_lon, nav_lat, u_obj.variables['vozocrtx'][i, :, :], lon, lat)
    v[i, :, :]  = interp_xy(nav_lon, nav_lat, v_obj.variables['vomecrty'][i, :, :], lon, lat)
    for j in range(50):
        Ba[i, j, :, :] = interp_xy(nav_lon, nav_lat, ptrc_obj.variables['Ba'][i, j, :, :], lon, lat)

save_var = { 'lon': lon, 'lat': lat, 'Ba': Ba, 'u': u, 'v': v}
scipy.io.savemat('../_data/Exchange/REMAP_Ba_UV.mat', mdict=save_var)

