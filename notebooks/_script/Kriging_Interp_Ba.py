import sys
sys.path.insert(0, '../_libs/')
import glob
#import pyproj
import scipy.io
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from NEMO_tools import int_profile, bin_profile, find_inland, reporj_NEMOgrid
from pykrige.ok import OrdinaryKriging
import pykrige.kriging_tools as kt

def check_bounds(out, bounds_low, bounds_up):
    '''
    "What out of bounds becomes nan"
    '''
    out_adj = np.copy(out)
    for i in range(133):
        for j in range(17):
            if((out[i, j]<bounds_low[j]) | (out[i, j]>bounds_up[j])):
                out_adj[i, j]=np.nan
    return out_adj
def remove_nan(x, y, data):
    '''
    "Kriging does not like nan"
    '''
    data_temp = data[~np.isnan(data)]
    x_temp = x[~np.isnan(data)]
    y_temp = y[~np.isnan(data)]
    return x_temp, y_temp, data_temp

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/GIPY11_Ba.mat')
GIPY_lons = MAT['Ba_records'][:, 0]
GIPY_lats = MAT['Ba_records'][:, 1]
GIPY_deps = MAT['Ba_records'][:, 2]
GIPY_Ba   = MAT['Ba_records'][:, 3]
#
GIPY_lons[GIPY_lons>180]=GIPY_lons[GIPY_lons>180]-360

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/BGEP_Ba.mat')
x03 = MAT['Ba2003'][:, 1]; y03 = MAT['Ba2003'][:, 0]; z03 = MAT['Ba2003'][:, 2]; Ba03 = MAT['Ba2003'][:, 3]
x04 = MAT['Ba2004'][:, 1]; y04 = MAT['Ba2004'][:, 0]; z04 = MAT['Ba2004'][:, 2]; Ba04 = MAT['Ba2004'][:, 3]
x04[Ba04<0]=np.nan; y04[Ba04<0]=np.nan; z04[Ba04<0]=np.nan; Ba04[Ba04<0]=np.nan
x05 = MAT['Ba2005'][:, 1]; y05 = MAT['Ba2005'][:, 0]; z05 = MAT['Ba2005'][:, 2]; Ba05 = MAT['Ba2005'][:, 3]

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/HLY_Ba.mat')
xH = np.squeeze(MAT['lon'][:]); yH = np.squeeze(MAT['lat'][:]); zH = np.squeeze(MAT['dep'][:]); BaH = np.squeeze(MAT['Ba'][:])

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/CBL_Ba.mat')
xC = np.squeeze(MAT['lon'][:]); yC = np.squeeze(MAT['lat'][:]); zC = np.squeeze(MAT['dep'][:]); BaC = np.squeeze(MAT['Ba'][:])

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/ARK09_Ba.mat')
xA9 = np.squeeze(MAT['lon'][:]); yA9 = np.squeeze(MAT['lat'][:]); zA9 = np.squeeze(MAT['dep'][:]); BaA9 = np.squeeze(MAT['Ba'][:])

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/ARK14_Ba.mat')
xA14 = np.squeeze(MAT['lon'][:]); yA14 = np.squeeze(MAT['lat'][:]); zA14 = np.squeeze(MAT['dep'][:]); BaA14 = np.squeeze(MAT['Ba'][:])

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/NPEO_Ba.mat')
xNPEO = np.squeeze(MAT['lon'][:]); yNPEO = np.squeeze(MAT['lat'][:]); zNPEO = np.squeeze(MAT['dep'][:]); BaNPEO = np.squeeze(MAT['Ba'][:])

MAT = scipy.io.loadmat('../../../Arctic-obs/MATLAB files/ARC_Ba.mat')
xARC = np.squeeze(MAT['lon'][:]); yARC = np.squeeze(MAT['lat'][:]); zARC = np.squeeze(MAT['dep'][:]); BaARC = np.squeeze(MAT['Ba'][:])

#nc_name=glob.glob('/ocean/yingkai/GEOTRACES/NEMO-CODE/NEMOGCM/CONFIG/OFF_TEST/EXP00/EXP01_1m_00010101_00041001_ptrc_T.nc')
#nc_obj=nc.Dataset(nc_name[0])
#deptht = nc_obj.variables['deptht'][:]
ANHA4_MAT=scipy.io.loadmat('../_data/Exchange/coord_ANHA4.mat')
deptht = ANHA4_MAT['nav_lev'][:]

BGEP_lons = np.hstack((x03, x04, x05))
BGEP_lats = np.hstack((y03, y04, y05))
BGEP_deps = np.hstack((z03, z04, z05))
BGEP_Ba   = np.hstack((Ba03, Ba04, Ba05))

x_all = np.hstack((GIPY_lons, x03, x04, x05, xH, xC, xA14, xNPEO, xARC))
y_all = np.hstack((GIPY_lats, y03, y04, y05, yH, yC, yA14, yNPEO, yARC))
z_all = np.hstack((GIPY_deps, z03, z04, z05, zH, zC, zA14, zNPEO, zARC))
Ba_all = np.hstack((GIPY_Ba, Ba03, Ba04, Ba05, BaH, BaC, BaA14, BaNPEO, BaARC))

tar_dep = 5250
dep_surf = deptht[deptht<tar_dep]

x_surf = x_all[z_all<tar_dep]; y_surf = y_all[z_all<tar_dep]; z_surf = z_all[z_all<tar_dep]; Ba_surf = Ba_all[z_all<tar_dep]

locx, locy, out_surf = int_profile(x_surf, y_surf, z_surf, Ba_surf, dep_surf, thres=2000)

x  = {'GIPY':GIPY_lons, 'BGEP':BGEP_lons, 'HLY':xH , 'CBL':xC ,'ARK14':xA14 , 'NPEO':xNPEO, 'ARC':xARC}
y  = {'GIPY':GIPY_lats, 'BGEP':BGEP_lats, 'HLY':yH , 'CBL':yC , 'ARK14':yA14 , 'NPEO':yNPEO, 'ARC':yARC}
z  = {'GIPY':GIPY_deps, 'BGEP':BGEP_deps, 'HLY':zH , 'CBL':zC , 'ARK14':zA14 , 'NPEO':zNPEO, 'ARC':zARC}
Ba = {'GIPY':GIPY_Ba  , 'BGEP':BGEP_Ba  , 'HLY':BaH, 'CBL':BaC, 'ARK14':BaA14, 'NPEO':BaNPEO, 'ARC':BaARC}
x_cut   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
y_cut   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
z_cut   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
Ba_cut  = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
x_int   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
y_int   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
z_int   = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
Ba_int  = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
x_trans = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
y_trans = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
x_cord  = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}
y_cord  = {'GIPY':[], 'BGEP':[], 'HLY':[], 'CBL':[], 'ARK14':[], 'NPEO':[], 'ARC':[]}

keys = ['GIPY', 'BGEP', 'HLY', 'CBL', 'ARK14', 'NPEO', 'ARC']

for i in keys:
    x_cut[i] = x[i][z[i]<tar_dep]; y_cut[i] = y[i][z[i]<tar_dep]
    z_cut[i] = z[i][z[i]<tar_dep]; Ba_cut[i] = Ba[i][z[i]<tar_dep]
    x_int[i], y_int[i], Ba_int[i] = int_profile(x_cut[i], y_cut[i], z_cut[i], Ba_cut[i], dep_surf, thres=2000)
    
fig = plt.figure(figsize=(8, 10))
# Axis
ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)
ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)
AX = [ax1, ax2, ax3, ax4]
for i in range(4):
    AX[i].set_xlim(20, 140)
    AX[i].xaxis.set_tick_params(size=0)
    AX[i].yaxis.set_tick_params(size=0)
# Fig1
ax1.set_ylim(-10, 500); ax1.invert_yaxis()
for i in range(len(locx)):
    ax1.plot(Ba_all[x_all==locx[i]], z_all[x_all==locx[i]], 'o:')
ax1.set_title('(a.1) Original, surface part', fontsize=12)
# Fig2
ax2.set_ylim(-10, 500); ax2.invert_yaxis()
ax2.plot(out_surf.T, dep_surf, 'o:')
ax2.set_title('(b.1) Interpolated, surface part')
# Fig3
ax3.set_ylim(500, tar_dep); ax3.invert_yaxis()
for i in range(len(locx)):
    ax3.plot(Ba_all[x_all==locx[i]], z_all[x_all==locx[i]], 'o:')
ax3.set_title('(a.2) Original, deep ocean', fontsize=12)
# Fig4
ax4.set_ylim(500, tar_dep); ax4.invert_yaxis()
ax4.plot(out_surf.T, dep_surf, 'o:')
ax4.set_title('(b.2) Interpolated, deep ocean')

p1=Basemap(projection='npstere', resolution='l', boundinglat=0, lon_0=90)
# create frame
res=40
xylim = [1e7, 1.7e7]
listx_trans = np.linspace(xylim[0], xylim[1], res)
listy_trans = np.linspace(xylim[0], xylim[1], res)
gridx_trans, gridy_trans = np.meshgrid(listx_trans, listy_trans)
# convert frame back to lat/lon
listx, listy = p1(listx_trans, listy_trans, inverse=True)
gridx, gridy = p1(gridx_trans, gridy_trans, inverse=True)

# for CTD's locs
x_all_trans, y_all_trans = p1(locx, locy)
for i in keys:
    x_trans[i], y_trans[i] = p1(x_int[i], y_int[i])
    
num_layer = len(deptht)-4
Ba_ini = np.empty([res, res, num_layer])
for i in range(num_layer):
    x_temp, y_temp, Ba_temp = remove_nan(x_all_trans, y_all_trans, out_surf[:, i])
    # plot cruises saperately
    for j in keys:
        x_cord[j], y_cord[j], _ = remove_nan(x_trans[j], y_trans[j], Ba_int[j][:, i])
        
    print('data points participated: {}'.format(len(Ba_temp)))# Check the the number of data
    # Kriging
    OK = OrdinaryKriging(x_temp, y_temp, Ba_temp, variogram_model='linear', verbose=False, enable_plotting=False)
    Ba_ini[:, :, i], ss = OK.execute('grid', listx_trans, listy_trans) # (gridx_trans, gridy_trans, Ba_int)
    # Plot
    #fig=plt.figure(figsize=(10, 10))
    #ax=fig.gca()
    #ax.set_xlim(xylim[0], xylim[1])
    #ax.set_ylim(xylim[0], xylim[1])
    #CS = ax.pcolor(gridx_trans, gridy_trans, Ba_ini[:, :, i], vmin=30, vmax=100, cmap=plt.cm.gist_ncar_r)
    #ax.plot(x_cord['GIPY'] , y_cord['GIPY'] , 'kx', ms=6, mew=2.5, label='GIPY 2007')
    #ax.plot(x_cord['BGEP'] , y_cord['BGEP'] , 'k^', ms=6, mew=0.5, label='BGEP 2003-2005')
    #ax.plot(x_cord['NPEO'] , y_cord['NPEO'] , 'k*', ms=9, mew=1.0, label='NPEO 2001-2005')
    #ax.plot(x_cord['HLY']  , y_cord['HLY']  , 'k+', ms=8, mew=2.5, label='HLY0301 2003')
    #ax.plot(x_cord['CBL']  , y_cord['CBL']  , 'k<', ms=6, mew=0.5, label='CBL 2002')
    #ax.plot(x_cord['ARK14'], y_cord['ARK14'], 'k>', ms=6, mew=0.5, label='ARKTISXIV-II 1998')
    #ax.plot(x_cord['ARC'], y_cord['ARC'], 'kv', ms=6, mew=0.5, label='ARCSS107-1993')
    #CBar = plt.colorbar(CS, shrink=0.5)
    #CBar.set_label('Interpolated diss.Ba con. ( 1E-6 mol/L )', fontsize=12)
    #CBar.ax.tick_params(axis='y', length=0)
    #LG = ax.legend(numpoints=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.); LG.draw_frame(False)
    #ax.set_title('Interpolated | deptht('+str(i)+')', fontsize=15, y=1.025)    
    #fig.savefig('Barium_layer_'+str(i)+'modified.pdf', dpi=400, orientation='portrait', papertype='a4', format='pdf',
    #        transparent=False, bbox_inches='tight', pad_inches=0)

#coordinate_name=glob.glob('/ocean/yingkai/GEOTRACES/FORCING/ORCA2_LIM_nemo_v3.4/coordinates*.nc')
#coordinate_obj=nc.Dataset(coordinate_name[0])
#nav_lon=coordinate_obj.variables['nav_lon'][:]
#nav_lat=coordinate_obj.variables['nav_lat'][:]
#
lon_list = np.linspace(-180, 180, 360)
lat_list = np.linspace(60, 90, 30)
lonxy, latxy = np.meshgrid(lon_list, lat_list)

Ba_ini_xy = np.empty([np.size(lonxy, 0), np.size(lonxy, 1), num_layer])
#Ba_ini_orca = np.empty([np.size(nav_lon, 0), np.size(nav_lon, 1), num_layer])
#hit = find_inland(nav_lon, nav_lat)
#hit = np.zeros(latxy.shape)

for i in range(num_layer):
    int_temp = reporj_NEMOgrid(gridx, gridy, Ba_ini[:, :, i], lonxy, latxy, method='linear')
    int_temp_fill = reporj_NEMOgrid(gridx, gridy, Ba_ini[:, :, i], lonxy, latxy, method='nearest')
    #int_temp, lon_list = addcyclic(int_temp, lon_list)
    int_temp[np.isnan(int_temp)] = int_temp_fill[np.isnan(int_temp)]
    #int_temp[hit==1]=np.nan
    Ba_ini_xy[:, :, i] = int_temp
#for i in range(num_layer):
#    int_temp_orca = reporj_NEMOgrid(lonxy, latxy, Ba_ini_xy[:, :, i], nav_lon, nav_lat, method='linear')
#    Ba_ini_orca[:, :, i] = int_temp_orca
    
Ba_ini_masked = np.ma.masked_where(np.isnan(Ba_ini_xy), Ba_ini_xy)

#plt.pcolor(lonxy, latxy, Ba_ini_masked[:, :, 0])
#ax = plt.gca(); ax.set_xlim(-180, 180), ax.set_ylim(60, 90)

#ETOPO2_Arctic=scipy.io.loadmat('_libs/ETOPO2_Arctic.mat')
#lon_arctic=ETOPO2_Arctic['lon_arctic']
#lat_arctic=ETOPO2_Arctic['lat_arctic']
#topo_arctic=ETOPO2_Arctic['topo_arctic']

## Adjust resolution
#res_unit=5
#lon_arctic=lon_arctic[0:-1:res_unit, 0:-1:res_unit]
#lat_arctic=lat_arctic[0:-1:res_unit, 0:-1:res_unit]
#topo_arctic=topo_arctic[0:-1:res_unit, 0:-1:res_unit]*-1

#clevs=[1000, 2000, 3000]
#for num in range(num_layer):
#    # Locations of CTDs
#    for j in keys:
#        x_cord[j], y_cord[j], _ = remove_nan(x_int[j], y_int[j], Ba_int[j][:, num])
#    # Figures
#    fig=plt.figure(figsize=(10, 10)); ax=plt.gca()
#    proj=Basemap(projection='npstere', resolution='l', boundinglat=55, lon_0=90, round=True, ax=ax)
#    proj.drawmeridians(np.arange(0, 360, 60), labels=[1, 1, 1, 1], fontsize=10, latmax=90,linewidth=0)
#    proj.fillcontinents(color=[0.5, 0.5, 0.5], lake_color=None)
#    proj.drawcoastlines(linewidth=1.5, color='k')
#    gridx, gridy = proj(lonxy, latxy)
#    topox, topoy = proj(lon_arctic, lat_arctic)
#    for j in keys:
#        x_cord[j], y_cord[j] = proj(x_cord[j], y_cord[j])
#    CS = proj.pcolor(gridx, gridy, Ba_ini_masked[:, :, num], vmin=30, vmax=100, cmap=plt.cm.gist_ncar_r)
#    CS2 = proj.contour(topox, topoy, topo_arctic, clevs, linestyles='-', linewidths=1, colors=('k',), alpha=0.75)
#    CS2.collections[0].set_label('1, 2, 3 km bathymetry')
#    proj.plot(x_cord['GIPY'] , y_cord['GIPY'] , 'kx', ms=6, mew=2.5, label='GIPY11 2007')
#    proj.plot(x_cord['BGEP'] , y_cord['BGEP'] , 'k^', ms=6, mew=0.5, label='BGEP 2003-2005')
#    proj.plot(x_cord['NPEO'] , y_cord['NPEO'] , 'k*', ms=9, mew=1.0, label='NPEO 2001-2005')
#    proj.plot(x_cord['HLY']  , y_cord['HLY']  , 'k+', ms=8, mew=2.5, label='HLY0301 2003')
#    proj.plot(x_cord['CBL']  , y_cord['CBL']  , 'k<', ms=6, mew=0.5, label='CBL 2002')
#    proj.plot(x_cord['ARK14'], y_cord['ARK14'], 'k>', ms=6, mew=0.5, label='ARKTISXIV-II 1998')
#    proj.plot(x_cord['ARC']  , y_cord['ARC']  , 'kv', ms=6, mew=0.5, label='ARCSS107-1993')
#    #proj.plot(gridx.T, gridy.T, 'k-', lw=0.25)
#    #proj.plot(gridx, gridy, 'k-', lw=0.25)
#    cbaxes = fig.add_axes([0.975, 0.15, 0.02, 0.5]) 
#    CBar = plt.colorbar(CS, ax=ax, cax=cbaxes) 
#    CBar.set_label('Interpolated diss.Ba con. ( 1E-6 mol/L )', fontsize=12)
#    CBar.ax.tick_params(axis='y', length=0)
#    LG = ax.legend(numpoints=1, bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.); LG.draw_frame(False)
#    ax.set_title('Initial field in a regular lat/lon frame, depth='+str(deptht[num])+'m', fontsize=15, y=1.025)
#    fig.savefig('Initial_field'+str(num)+'.pdf', dpi=600, orientation='portrait', papertype='a4', format='pdf',
#                transparent=False, bbox_inches='tight', pad_inches=0)

# Save
save_var = { 'Ba_ini_xy': Ba_ini_xy, 'lon': lonxy, 'lat': latxy}
scipy.io.savemat('../_data/Exchange/Temp_Ba_int_ANHA4.mat', mdict=save_var)

