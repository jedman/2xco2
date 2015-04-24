import numpy as np
import matplotlib.pyplot as plt
from vstats import *
import netCDF4 as nc 

class TimeSeries:
    def __init__(self, time, z = []):
        '''minimum thing to define a time series. can also include a z dimension'''
        self.time = time
        self.z = z

    def read1d(self,vars,filepath,nfiles ):
        '''vars: list of variables to read in
        filepath: beginning of file string
        nfiles: list of file numbers to read in '''
        for index, var in enumerate(vars):
            for i in range(nfiles):
                filename = filepath + str(i) + ".nc"
                tmp = nc.Dataset(filename)
                tmp_var = getvar(var, tmp)
                if i == 0:
                    prevar = tmp_var
                else:
                    prevar = np.append(prevar, tmp_var) # collect the variable
            if (index == 0):
                self.vars_1d = {var:prevar} # create the dictionary the variables
            else:
                self.vars_1d.update({var:prevar}) # update dictionary with each new variable

    def read2d(self,vars,filepath, nfiles):
        '''vars: list of variables to read in
          filepath: beginning of file string
          nfiles: list of file numbers to read in '''
        for index, var in enumerate(vars):
            tmp_data = np.zeros((len(self.time), len(self.z)))
            for i in range(nfiles):
                filename = filepath + str(i) + ".nc"
                tmp = nc.Dataset(filename)
                tmp_time = getvar('time',tmp)
                tmp_var = getvar(var, tmp)
                if i == 0:
                    chunk = range(0,len(tmp_time))
                    tmp_data[chunk, :] = tmp_var
                else:
                    chunk = range(max(chunk)+1, max(chunk)+len(tmp_time)+1)
                    tmp_data[chunk,:] = tmp_var
            if (index == 0):
                self.vars_2d = {var:tmp_data} # create the dictionary the variables
            else:
                self.vars_2d.update({var:tmp_data})


def readtime(filepath,nfiles, ndays = 10., offset = 0. ):
    '''Special case of read1d.
    filepath: beginning of file string
    nfiles: list of file numbers to read in
    ndays: number of days in each file
    offset: day at start of time series'''

    time = np.zeros(0)
    for i in range(nfiles):
        filename = filepath + str(i) + ".nc"
        tmp = nc.Dataset(filename)
        time_tmp = getvar('time',tmp) + ndays*i + offset
        time = np.append(time,time_tmp)
    return time
