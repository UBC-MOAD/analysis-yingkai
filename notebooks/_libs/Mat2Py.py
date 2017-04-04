'''
========================================================
Mat2Py
--------------------------------------------------------
A collection of useful Python functions for MATLAB users
                ----- CCAR Modelling Team
                -----   Yingkai (Kyle) Sha
========================================================
2014/12/24 File created
2015/02/22 Convert listed Python datetime to MATLAB datenum
'''
import datetime
import numpy as np

def datenum2datetime(matlab_datenum):
    '''
    =======================================================================
    Convert MATLAB datenum value to Python datetime object
                            ----- created on 2014/12/24, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        python_datetime = datenum2datetime(matlab_datenum)
    -----------------------------------------------------------------------
    Input:
            matlab_datenum: 1D array
    Output:
            datetime objects
    ======================================================================= 
    '''
    matlab_datenum=matlab_datenum.astype(int)
    python_datetime=[]
    for i in range(len(matlab_datenum)):
        temp = datetime.datetime.fromordinal(matlab_datenum[i]) + \
                            datetime.timedelta(days=matlab_datenum[i]%1) - datetime.timedelta(days = 366)
        python_datetime.append(temp)
    return python_datetime

def datetime2datenum_single(dt):
    '''
    =======================================================================
    Convert a single Python datetime object to MATLAB datenum
                            ----- created on 2014/12/24, Yingkai (Kyle) Sha
    ======================================================================= 
    '''
    #toord = dt.toordinal()
    mdn = dt + datetime.timedelta(days = 366)
    frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac
    
def datetime2datenum(dt):
    '''
    =======================================================================
    Convert a list of Python datetime object to MATLAB datenum
                            ----- created on 2015/2/22, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        matlab_datenum = datenum2datetime(python_datetime)
    -----------------------------------------------------------------------
    See datenum2datetime
    ======================================================================= 
    '''
    #dt = np.squeeze(dt)
    temp = np.zeros(len(dt))
    for i in range(len(dt)):
        temp[i] = datetime2datenum_single(dt[i])
    return temp

