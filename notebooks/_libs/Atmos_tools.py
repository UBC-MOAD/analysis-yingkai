'''
============================================
Atmos_tools
--------------------------------------------
A collection of functions in Atmospheric Sciences
                ----- Author: Yingkai Sha
============================================
2014/12/21 File created
'''

import numpy as np

def central_diff(T):
    '''
    % =========================================================================== %
    % Central diffference for inner domain; forward/backward difference for edges 
    %   build-in function
    % =========================================================================== %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/02/25
    % =========================================================================== %
    Rewrite in Python 2.7, 2016/04/20
        yingkai@eos.ubc.ca
    '''
    # Get the length of x, y dimension
    M = np.size(T, axis=0)
    N = np.size(T, axis=1)
    # Allocation
    dx = np.zeros(T.shape)
    dy = np.zeros(T.shape)
    # dx
    dx[:, 0] = T[:, 1]-T[:, 0]
    dx[:,-1] = T[:,-1]-T[:,-2]
    for i in range(1, N-1):
        dx[:, i] = 0.5*(T[:, i+1]-T[:, i-1])
    # dy
    dy[0, :] = T[1, :]-T[0, :]
    dy[-1,:] = T[-1,:]-T[-2,:]
    for i in range(1, M-1):
        dy[i, :] = 0.5*(T[i+1, :]-T[i-1, :])
    return dx, dy

def dx_atmos(lon, lat):
    '''
    % ======================================================================= %
    % Calculate the grid distance on the earth,
    %   dy=dx*h1, h1=R*sin(phi)=R*cos(lat)
    %       (h1 is the scale factor in spherical coordinate)
    %   build-in function
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/02/25
    % ======================================================================= %
    Rewrite in Python 2.7, 2016/04/20
        yingkai@eos.ubc.ca
    '''
    R = 6.3781e6 # earth radius (m)
    dx, _ = central_diff(lon)
    dx = dx * (np.pi/180.0)*R*np.cos(lat*np.pi/180.0)
    return dx


def dy_atmos(lat):
    '''
    % ======================================================================= %
    % Calculate grid distance on Earth, 
    %   dy=dy*h2, h2=R
    %       (h2 is scale factor in spherical coordinate)
    %   build-in function
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/02/25
    % ======================================================================= %
    Rewrite in Python 2.7, 2016/04/20
        yingkai@eos.ubc.ca
    '''
    R = 6.3781e6 # earth radius (m)
    _, dy = central_diff(lat)
    dy = dy*(np.pi/180.0)*R
    return dy

def curlz_atmos(lon, lat, u, v):
    '''
    % curlz=curlz_atmos(longitude, latitude, u, v)
    %   Calculate vorticity in vertical direction
    % ======================================================================= %
    % Input
    %   longitude: Longitude, deg
    %   latitude: Latitude, deg
    %   u: Zonal Wind, m/s
    %   v: Meditorial Wind, m/s
    % Output
    %   curlz: Vorticity in vertical direction, s^-1
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/2/27
    % ======================================================================= %
    Rewrite in Python 2.7, 2016/04/25
        yingkai@eos.ubc.ca
    '''
    R = 6.3781e6; # earth's radius (m)
    dx = dx_atmos(lon, lat)
    dy = dy_atmos(lat)
    _, du = central_diff(u)
    dv, _ = central_diff(v)
    curlz = dv/dx-du/dy+u*np.tan(lat*np.pi/180.0)/R
    return curlz

def divh_atmos(lon, lat, u, v):
    '''
    % divh=divh_atmos(longitude, latitude, u, v)
    %   Calculate horizontal divergence
    % ======================================================================= %
    % Input
    %   longitude: Longitude, deg
    %   latitude: Latitude, deg
    %   u: Zonal Wind, m/s
    %   v: Meditorial Wind, m/s
    % Output
    %   divh: Horizontal Divergence, s^-1
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/2/27
    % ======================================================================= %
        Rewrite in Python 2.7, 2016/05/21
            yingkai@eos.ubc.ca
    '''
    R=6.3781e6; # earth radius (m)
    dx = dx_atmos(lon, lat);
    dy = dy_atmos(lat);
    du, _ = central_diff(u);
    _, dv = central_diff(v);
    divh = du/dx+dv/dy-v*np.tan(lat*np.pi/180.0)/R;
    #divh[np.abs(lat)==90]=np.nan;

    return divh


def advh_atmos(lon, lat, u, v, T):
    '''
    % ======================================================================= %
    % Calculate horizontal advection
    % Input
    %   longitude: Longitude, deg
    %   latitude: Latitude, deg
    %   u: Zonal Wind, m/s
    %   v: Meditorial Wind, m/s
    %   T: Scalar field
    % Output
    %   advh: Horizontal Advection of T, unit(T)*s^-1
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/02/27
    % ======================================================================= %
    Rewrite in Python 2.7, 2016/04/20
        yingkai@eos.ubc.ca
    '''
    R=6.3781e6 # earth radius (m)
    dx = dx_atmos(lon, lat)
    dy = dy_atmos(lat)
    dTx, dTy = central_diff(T)
    advh=-(v*dTy/dy+u*dTx/dx)
    return advh

from Atmos_tools import *
def grad_atmos(lon, lat, H):
    '''
    % [grad_x grad_y]=grad_atmos(longitude, latitude, H)
    %   Calculate gradient of a 2-D variable.
    % ======================================================================= %
    % Input
    %   longitude: Longitude, deg
    %   latitude: Latitude, deg
    %   H: a scalar field
    % Output
    %   grad_x: \frac{\partial H}{\partial x}, unit(H)/m
    %   grad_y: \frac{\partial H}{\partial y}, unit(H)/m
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/2/27
    % ======================================================================= %
    Rewrite in Python 2.7, 2016/07/07
        yingkai@eos.ubc.ca
    '''
    dx = dx_atmos(lon, lat);
    dy = dy_atmos(lat);
    #for i in range(np_layer+1, np.size(lat, 0)):
    #    dy[i, :] = dy[np_layer, :]
    grad_x, grad_y = central_diff(H);
    grad_x = grad_x/dx;
    grad_y = grad_y/dy;
    return grad_x, grad_y

def geo_wind(lon, lat, H):
    '''
    % [ug vg]=geostrophic_wind(longitude, latitude, H)
    %   Calculate geostrophic wind (in quasi-geostrophic balance).
    % ======================================================================= %
    % Input
    %   longitude: Longitude, deg
    %   latitude: Latitude, deg
    %   H: Geopotential Height, m
    % Output
    %   [ug vg]: Geostrophic Wind components, [zonal meditorial], m/s
    % Note
    %   More information about the arrangement of input data, or if you get ill
    %   results (e.g. NaN or Inf value), see "check_input".
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/2/27
    % ======================================================================= %
    '''
    omega = 7.292*1e-5; # rotational velocity
    g = 9.8; # m/s^2 % gravity
    f = 2*omega*np.sin(lat*np.pi/180.0); # geostrophic parameter
    grad_x, grad_y = grad_atmos(lon, lat, H);
    ug = -(g/f)*grad_y;
    vg = (g/f)*grad_x;
    #ug[np.abs(lat)>88.5] = np.nan;
    #vg[np.abs(lat)>88.5] = np.nan;
    return ug, vg


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def equi(centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    return np.array(X), np.array(Y)


    
def EOF(H, nmode=10, ndim=3, reverse=1):
    '''
    Converted from MATLAB to Python 2.7 code @ 2015/06/15 - YKS
     + ndim: [LAT, LON, TIME] data (=3) or [MAP, TIME] data (=2)
     + reverse: normalized spatial pattern (=0), normalized PC (=1)
    % ======================================================================= %
    % Input
    %   H: Variable required for EOF comutation, H(LAT, LON, Time) 
    %       or H(Space, Time) is accepted.
    %   nmode: Number of modes output
    % Output
    %   EOFs: EOF Spatial Pattern
    %   PC: Timeseries cooresponding to each EOFs
    %   expvar: Explained variance
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/3/18
    % ======================================================================= %
    '''
    ##import scipy.linalg.eig as eig
    # Get the size of array
    if ndim == 3:
        LAT, LON, T = H.shape
    elif ndim == 2:
        LON, T = H.shape
        LAT = 1
    # Covarience
    H = np.reshape(H, [LAT*LON, T]).T
    R=np.dot(H, H.T); N = np.size(R, 0)
    # Allocation
    PC     = np.zeros([nmode, T]);
    expvar = np.zeros([nmode]);
    eof    = np.zeros([N, LAT*LON]);
    EOFs   = np.zeros([LAT, LON, nmode]);
    # Eigvector analysis
    L, E = np.linalg.eig(R)
    # Get modes
    E    = np.dot(H.T, E)
    #sq   = (np.sqrt(np.diag(L))).T
    #sq   = sq[0, :]
    sq = np.sqrt(L)
    E    = E/sq
    Z    = np.dot(E.T, H.T)
    for i in range(nmode):
        eof[i, :] = np.squeeze(E[:, i]).T
        PC[i, :]  = np.squeeze(Z[i, :])
    # Get expvar
    L = np.abs(L)
    dsum = np.sum(np.abs(L))
    # Output
    for i in range(nmode):
        expvar[i] = L[i]/dsum
        EOFs[:, :, i] = np.reshape(eof[i, :], [LAT, LON])
    if reverse==1:
        EOFs, PC = reverse_std(EOFs, PC, nmode)
    return EOFs, PC, expvar

def reverse_std(EOFs, PC, nmode):
    for i in range(nmode):
        STD = np.nanstd(PC[i, :])
        PC[i, :] = PC[i, :]/STD
        EOFs[:, :, i] = EOFs[:, :, i]*STD
    return EOFs, PC

def seasonal_cycle(data):
    '''
    Compute the average of 12 months in all years.
    Input monmean data, time must be the 1st dim.
    '''

    out = np.empty([12, np.size(data, 1), np.size(data, 2)])
    for mon in range(12):
        out[mon, :, :] = np.nanmean(data[mon:len(data):12, :, :], 0)
    return out


def seasonal_decomp1d(data, method=0):
    '''
    =======================================================================
    Remove the seasonal cycle from 1D data
                            ----- created on 2015/06/15, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        data = seasonal_decomp(...)
    -----------------------------------------------------------------------
    Input:
            data
            method: removal done by anomaly (=0) or normalize (=1)
    ======================================================================= 
    '''
    data2 = np.empty(data.shape)
    N = len(data)
    for mon in range(12):
        temp_data = np.nanmean(data[mon:N:12], 0)
        if method == 0:
            data2[mon:N:12] = data[mon:N:12]-temp_data
    return data2

def seasonal_decomp3d(data, method=0):
    '''
    =======================================================================
    Remove the seasonal cycle from 3D data
                            ----- created on 2015/06/15, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        data = seasonal_decomp(...)
    -----------------------------------------------------------------------
    Input:
            data: Time should be the first dim.
            method: removal done by anomaly (=0) or normalize (=1)
    ======================================================================= 
    '''
    data2 = np.empty(data.shape)
    N = len(data)
    for mon in range(12):
        temp_data = np.nanmean(data[mon:N:12, :, :], 0)
        if method == 0:
            data2[mon:N:12, :, :] = data[mon:N:12, :, :]-temp_data
    return data2

   
