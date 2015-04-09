# a bunch of routines for handling vertical stats files from DAM
# and calculating derivatives
import numpy as np 
def getvar(var,data): 
    '''grab var from a netcdf data set (data); read into memory'''
    variable = data.variables[var][:]
    return variable
def profile(var,tb=0, te=-1):
    '''make a profile averaged in time from tb to te'''
    prof_var = np.mean(var[tb:te,:],axis=0)  
    return prof_var 

def prime(var1):
    '''calculate the transient anomaly of var1, eg, var1-profile(var1)'''
    tmp = profile(var1)
    transient = var1-tmp 
    return transient 

def ddzp(prof1, z, sdo=False):
    ''' calculate ddz of some scalar profile on the scalar levels,\
       interpolating to the surface if sdo =True'''
    dz = np.zeros(len(z))  
    dz[0] = 0.5*(z[0]+z[1])
    for i in range(1,len(z)-1):
       dz[i] = 0.5*(z[i+1]-z[i-1])
    dz[-1]= dz[-2] # fudge for the top level -- don't know dzi[-1] 
    vflux = np.zeros(len(z)+1)
    for k in xrange(1,len(z)):
        vflux[k] = 0.5*(prof1[k]+prof1[k-1]) # value of prof1 at the interface k
        if(sdo): # for EvRTdv and Eldl, the surface flux is nonzero, so use the value at the interface
            vflux[0] = prof1[0]-(prof1[1]-prof1[0])/(z[1]-z[0])*z[0]  
    ddz_something = (vflux[1:] - vflux[0:-1])/dz[0:] 
    return ddz_something

def ddzi(prof1, z, sdo =False):
    ''' calculate ddz of some interface profile on the scalar levels'''
    dz = np.zeros(len(z))  
    dz[0] = 0.5*(z[0]+z[1])
    for i in range(1,len(z)-1):
       dz[i] = 0.5*(z[i+1]-z[i-1])
    dz[-1]= dz[-2] # fudge for the top level -- don't know dzi[-1] 
    vflux = np.zeros(len(z)+1)
    vflux[0:-1] = prof1
    ddz_something = (vflux[1:] - vflux[0:-1])/dz[0:] 
    return ddz_something

def scalep(prof,z):
    '''interpolate interface var to scalar levels defined by z '''
    dz = np.zeros(len(z))  
    dz[0] = 0.5*(z[0]+z[1])
    for i in range(1,len(z)-1):
       dz[i] = 0.5*(z[i+1]-z[i-1])
    dz[-1]= dz[-2]
    prof_scalar = np.zeros(len(z))
    for k in xrange(1,len(z)-1):
        prof_scalar[k] = prof[k]+((prof[k+1]-prof[k])/dz[k])*0.5*(z[k]-z[k-1])
            
    prof_scalar[0] = prof[0]+((prof[1]-prof[0])/dz[0])*(z[0]) 
    k = len(z)-1
    prof_scalar[len(z)-1] = prof[k]+((0.-prof[k])/dz[k])*(0.5*(z[k]-z[k-1]))
                                        
    return prof_scalar

def scale2int(var,z): 
    '''interpolate a scalar profile to the interfaces'''
    dz = makedz(z)
    varint = np.zeros(len(z))
    for k in range(1,len(z)):
        varint[k] = 0.5*(var[k]+var[k-1])
    varint[0] = var[0] - (var[1]-var[0])/(z[1]-z[0])*z[0]    
    return varint

def makedz(z):
    '''from some z vector of scalar levels, make the dz vector'''
    dzvec = np.zeros(len(z))  
    dzvec[0] = 0.5*(z[0]+z[1])
    for i in range(1,len(z)-1):
       dzvec[i] = 0.5*(z[i+1]-z[i-1])
    dzvec[-1]= dzvec[-2]
    return dzvec
