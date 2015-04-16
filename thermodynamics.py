##  THERMODYNAMICS  #####################################################
##                                                                     ##
##  These functions and constants are meant to be consistent with      ##
##  the definitions in Das Atmospharische Modell (DAM) (specifically,  ## 
##  thermo_mod.f90 and params_mod.f90). See also Romps (2008).         ##
##                                                                     ##
##  Python implementation J. Seeley, 2014                              ##
##                                                                     ## 
##  Revision history:                                                  ## 
##  August 2014 -- first version                                       ## 
##                                                                     ## 
##                                                                     ## 
#########################################################################

## Setup  ###############################################################
                                                                       ##
import numpy as np                                                     ##
import scipy.optimize as opt                                           ##
                                                                       ##
#########################################################################

##  Constants  ##########################################################
                                                                       ##
##  Dry air                                                            ##
c_vd             = 719.             # J/kg/K                           ##
R_d              = 287.04           # J/kg/K                           ##
c_pd             = c_vd + R_d       # J/kg/K                           ##
                                                                       ##
##  Water vapor                                                        ##
c_vv             = 1418.            # J/kg/K                           ##
R_v              = 461.4            # J/kg/K                           ##
c_pv             = c_vv + R_v       # J/kg/K                           ##
                                                                       ##
##  Liquid water                                                       ##
c_vl             = 4216.            # J/kg/K                           ##
                                                                       ##
##  Solid water                                                        ##
c_vs             = 2106.            # J/kg/K                           ##    
                                                                       ##
##  Reference temperatures and pressures                               ##
T_0              = 273.16           # K                                ## 
p_0              = 1.e5             # Pa                               ##
e_0              = 611.65           # Pa                               ##
                                                                       ##
##  Energies, enthalpies, entropies                                    ##
L_v              = 2260000.         # J/kg                             ##
E_0v             = 2374000.         # J/kg                             ##
E_0s             = 333700.          # J/kg                             ##
s_0v             = E_0v/T_0 + R_v   # J/kg/K                           ##
s_0s             = E_0s/T_0         # J/kg/K                           ##
                                                                       ##
##  other                                                              ##
gg               = 9.81             # m/s^2                            ##
eps              = R_v/R_d - 1.     # Unitless                         ##
                                                                       ##
#########################################################################
 
##  Saturation vapor pressure with respect to liquid  ###################
def pvstar_l(T):
    return (e_0*(T/T_0)**((c_pv - c_vl)/R_v))* \
           np.exp(((E_0v - T_0*(c_vv - c_vl))/R_v)* \
           (1./T_0 - 1./T))

##  Saturation vapor pressure with respect to solid  ####################
def pvstar_s(T):
    return (e_0*(T/T_0)**((c_pv - c_vs)/R_v))* \
           np.exp(((E_0v + E_0s - T_0*(c_vv - c_vs))/R_v)* \
           (1./T_0 - 1./T))

##  Saturation vapor pressure with transition at triple point temp  #####
def pvstar(T):
    if np.isscalar(T):
        if T > T_0:
            return pvstar_l(T)
        else:
            return pvstar_s(T)
    elif isinstance(T,np.ndarray):
        shape = T.shape
        pvstar_array = np.nan*np.empty(shape)
        pvstar_array[T>T_0]  = pvstar_l(T[T>T_0])
        pvstar_array[T<=T_0] = pvstar_s(T[T<=T_0]) 
    	return pvstar_array
    else:
        print 'Error: supplied temperature data not array or scalar'
        return
         	
##  Saturation specific humidity  #######################################
def qvstar(p,T,ql=0.,qi=0.,qt=0.,liquid=False,solid=False):
    if liquid:
        dcv = c_vv - c_vl
        dcp = c_pv - c_vl
        E0  = E_0v
    elif solid:
        dcv = c_vv - c_vs
        dcp = c_pv - c_vs
        E0  = E_0v + E_0s
    else:
        warm = 0.5 + 0.5*np.sign(T-T_0)
        dcv  = (c_vv - c_vs) + warm*((c_vv - c_vl) - (c_vv - c_vs))
        dcp  = (c_pv - c_vs) + warm*((c_pv - c_vl) - (c_pv - c_vs))
        E0   = (E_0v + E_0s) - warm*E_0s
    if qt != 0.:
        return (1. - qt)*(R_d/R_v)/(p/e_0 * (T_0/T)**(dcp/R_v) * \
             np.exp(-(E0 - dcv*T_0)/R_v*(1./T_0 - 1./T)) - 1.)
    else:
        ql2 = ql
        qi2 = qi
        return (1. - ql2 - qi2)/(R_v*p/(R_d*e_0)*(T_0/T)**(dcp/R_v)* \
               np.exp(-(E0 - dcv*T_0)/R_v*(1./T_0 - 1./T)) - eps)

##  Equivalent potential temperature  ###################################
def theta_e(p,T,qv,ql=0.,qi=0.,p_0=p_0):
    qa = 1. - qv - ql - qi
    rv = qv/qa
    rl = ql/qa
    ri = qi/qa
    pa = (R_d/(R_d + rv*R_v))*p
    pv = rv*(R_v/(R_d + rv*R_v))*p
    theta_e_val = T*(p_0/pa)**(R_d/c_pd)*(T/T_0)** \
                  ((rv*c_pv + rl*c_vl + ri*c_vs)/c_pd)* \
                  (e_0/pv)**(rv*R_v/c_pd)*np.exp((rv*s_0v - ri*s_0s)/c_pd)
    return theta_e_val

##  Begin adiabatic state  ##############################################
def dry_tabs_from_p(p,qt,thetae0):                                     
    rv = qt/(1.-qt)                                                    
    pa = (1.-qt)*(R_d/((1.-qt)*R_d + qt*R_v))*p                        
    pv = qt*(R_v/((1.-qt)*R_d + qt*R_v))*p                             
    return (thetae0*T_0**(rv*c_pv/c_pd)*(pa/p_0)** \
           (R_d/c_pd)*(pv/e_0)**(rv*R_v/c_pd)* \
           np.exp(-rv*s_0v/c_pd))**(1./(1.+rv*c_pv/c_pd))             
                                                                       
def dry_qv_surplus(p,qt,thetae0):                                      
    return qt - qvstar(p,dry_tabs_from_p(p,qt,thetae0))                
                                                                       
# Theta_e surpluses at the triple point								   
def tp_thetae_surplus_no_solid(p,qt,thetae0):						   
    qv = (1.-qt)*(R_d/R_v)*e_0/(p-e_0)								   
    ql = qt - qv													   
    qi = 0.															   
    return theta_e(p,T_0,qv,ql,qi)-thetae0							   	
    																   
def tp_thetae_surplus_no_liquid(p,qt,thetae0):						   
    qv = (1.-qt)*(R_d/R_v)*e_0/(p-e_0)								   
    ql = 0.															   
    qi = qt - qv													   
    return theta_e(p,T_0,qv,ql,qi)-thetae0							   
																	   
# Objective functions for root solver								   
def thetae_func_vl(T,p,qt,thetae0):									   
    # Find pv with no condensate									   
    pv = qt*(R_v/((1.-qt)*R_d + qt*R_v))*p							   
    if pv > pvstar(T):												   
        qv = qvstar(p,T,qt=qt)										   
    else:															   
        qv = qt														   
    ql = qt - qv													   
    qi = 0.															   
    return theta_e(p,T,qv,ql,qi)-thetae0							   
																	   
def thetae_func_vls(ql,p,qt,thetae0,plow,phigh):					   
    if (not (p >= plow and p <= phigh)):							   
        print 'Error in thetae_func_vls: p not in isothermal range'    
        return NA													   
    qv = (1. - qt)*(R_d/R_v)*(e_0/(p - e_0))						   
    qi = qt - qv - ql												   
    return theta_e(p,T_0,qv,ql,qi)-thetae0							   
																	   
def thetae_func_vs(T,p,qt,thetae0):									   
    # Find pv with no condensate									   
    pv = qt*(R_v/((1.-qt)*R_d + qt*R_v))*p							   
    if pv > pvstar(T):												   
        qv = qvstar(p,T,qt=qt)										   
    else:															   
        qv = qt														   
    ql = 0.															   
    qi = qt - qv													   
    return theta_e(p,T,qv,ql,qi)-thetae0							   
																	   
def find_saturated_regime(plow,phigh,p,qt,thetae0):					   
    na_val = -9999.													   
    if (plow == na_val and phigh == na_val):						   
        def obj_func(x):											   
            return thetae_func_vl(x,p,qt,thetae0)					   
        T = opt.fsolve(obj_func,10.)[0]								   
        if (T > T_0):												   
            regime = 'vl'											   
        else:														   
            regime = 'vs'											   
    elif plow == na_val:											   
        if p > phigh:												   
            regime = 'vl'											  
        else:										
            regime = 'vls'
    elif phigh == na_val:
        if p > plow:
            regime = 'vls'
        else:
            regime = 'vs'
    else:
        if p > phigh:
            regime = 'vl'
        elif p > plow:
            regime = 'vls'
        else:
            regime = 'vs'
    return regime

##  Actual function  ##
def adiabatic_state(thetae0,qt,p,verbose=False):
    pmin = 1.e3
    pmax = 2.e5
    na_val = -9999.
    if (p <= pmin or p >= pmax):
        print 'Warning in adiabatic_state: p='+str(p)+' is out of bounds'
    # Find psat for where the parcel becomes saturated
    always_unsat = False
    always_sat   = False
    if (dry_qv_surplus(pmin,qt,thetae0) < 0.):
        always_unsat = True
        psat  = na_val
        plow  = na_val
        phigh = na_val
    elif (dry_qv_surplus(pmax,qt,thetae0) > 0.):
        always_sat = True
        psat = na_val
    else:
        def obj_func(x):
            return dry_qv_surplus(x,qt,thetae0)
        psat = opt.fsolve(obj_func,p_0)[0]
    if (not always_unsat):
        if (not psat == na_val):
            ptop = psat
        else:
            ptop = pmax
        if (np.sign(tp_thetae_surplus_no_liquid(ptop,qt,thetae0)) != 
            np.sign(tp_thetae_surplus_no_liquid(pmin,qt,thetae0))):
            def obj_func(x):
                return tp_thetae_surplus_no_liquid(x,qt,thetae0)
            plow = opt.fsolve(obj_func,pmin)[0]
        else:
            plow = na_val
        if (np.sign(tp_thetae_surplus_no_solid(ptop,qt,thetae0)) != 
            np.sign(tp_thetae_surplus_no_solid(pmin,qt,thetae0))):
            def obj_func(x):
                return tp_thetae_surplus_no_solid(x,qt,thetae0)
            phigh = opt.fsolve(obj_func,p_0)[0]
        else:
            phigh = na_val
    else:
        phigh = na_val
        
    # Find the regime for the given p
    if always_unsat:
        regime = 'v'
    elif always_sat:
        regime = find_saturated_regime(plow,phigh,p,qt,thetae0)
    elif p > psat:
        regime = 'v'
    else:
        regime = find_saturated_regime(plow,phigh,p,qt,thetae0)
    
    # Solve for the state
    if (regime == 'v' or regime == 'vl'):
        def obj_func(x):
            return thetae_func_vl(x,p,qt,thetae0)
        T = opt.fsolve(obj_func,10.)[0]
        # Find pv with no condensate
        pv = qt*R_v/((1.-qt)*R_d + qt*R_v)*p
        if pv > pvstar(T):
            qv = qvstar(p,T,qt=qt)
        else:
            qv = qt
        ql = qt - qv
        qi = 0.
        rho = p/(((1.-qt)*R_d + qv*R_v)*T)
    elif regime == 'vls':
        T = T_0
        def obj_func(x):
            return thetae_func_vls(x,p,qt,thetae0,plow,phigh)
        ql = opt.fsolve(obj_func,0.)[0]
        qv = (1. - qt)*(R_d/R_v)*e_0/(p - e_0)
        qi = qt - qv - ql
        rho = p/(((1.-qt)*R_d + qv*R_v)*T)
    else:
        def obj_func(x):
            return thetae_func_vs(x,p,qt,thetae0)
        T = opt.fsolve(obj_func,10.)[0]
        # Find pv with no condensate
        pv = qt*R_v/((1.-qt)*R_d + qt*R_v)*p
        if pv > pvstar(T):
            qv = qvstar(p,T,qt=qt)
        else:
            qv = qt
        ql = 0.
        qi = qt - qv
        rho = p/(((1.-qt)*R_d + qv*R_v)*T)
     
    if verbose:
        print ''
        print 'T  = '+str(T) 
        print 'qv = '+str(qv)
        print 'ql = '+str(ql)
        print 'qi = '+str(qi)
        print 'qt = '+str(qv + ql + qi)
        print 'rho= '+str(rho)
        print ''
    # Check the solution
    if (np.abs(theta_e(p,T,qv,ql,qi)/thetae0 - 1.) > 1.e-6):
        print 'Error in adiabatic_state: reconstructed state gives '\
               +'theta_e = '+str(theta_e(p,T,qv,ql,qi))+' instead of '\
               +str(thetae0)
        print p
        print T
        print qv
        print ql
        print qi
    if (np.abs((qv + ql + qi)/qt - 1.) > 1.e-6):
        print 'Error in adiabatic_state: water not conserved.'
        print 'initial thetae: '+str(thetae0)
        print 'initial qt: '+str(qt)
        print 'initial p: '+str(p)
    return (T,qv,ql,qi,rho)
##  End adiabatic state  ################################################


def verbose_fallout(dqv,dql,dqi):
    if dqv < 0:
        source = ' vapor '
    else:
        source = ''
    if dql < 0:
        source = source + 'liquid'
    else:
        source = source + ''
    if dqi < 0:
        source = source + ' solid '
    else:
        source = source + ''
    if dqv > 0:
        sink = ' vapor '
    else:
        sink = ''
    if dql > 0:
        sink = sink + 'liquid'
    else:
        sink = sink + ''
    if dqi > 0:
        sink = sink + ' solid '
    else:
        sink = sink + ''
    return '['+source+'] converted into ['+sink+']'
        
    
def fallout(gamma,p,T,qv_series,ql_series,qi_series,verbose=False):
    thetae_init = theta_e(p,T,qv_series[1],ql_series[1],qi_series[1])
    qt_init = qv_series[1]+ql_series[1]+qi_series[1]
    qt_min = qvstar(p,T)
    dqv = qv_series[1]-qv_series[0]
    dql = ql_series[1]-ql_series[0]
    dqi = qi_series[1]-qi_series[0]
    dqcon = dql + dqi
    
    if verbose:
        print verbose_fallout(dqv,dql,dqi)
        print 'theta_e before fallout: '+str(thetae_init)
        print 'qt before fallout: '+str(qt_init)
        print 'Minimum qt: '+str(qt_min)
        print 'dqv: '+str(dqv)
        print 'dql: '+str(dql)
        print 'dqi: '+str(dqi)
        print 'new condensates: '+str(dqcon)
    # If no new condensates, no fallout
    if (dqcon <= 0.):
        if verbose:
            print 'No new condensates.'
        return (ql_series[1],qi_series[1],thetae_init,qt_init)
    else:
        # more liquid and ice
        if (dql >= 0. and dqi >= 0.):
            dql_new = (1. - gamma)*dql
            dqi_new = (1. - gamma)*dqi
        # less liquid, more ice
        elif (dql < 0. and dqi > 0.):
            dql_new = dql
            dqi_new = dqi - gamma*dqcon
        # more liquid, less ice
        elif (dql > 0. and dqi < 0.):   
            dql_new = dql - gamma*dqcon
            dqi_new = dqi_new
        else:
            print 'warning: other case'
        qv_new = qv_series[1]
        ql_new = ql_series[0] + dql_new
        qi_new = qi_series[0] + dqi_new
        qt_new = qv_series[1] + ql_new + qi_new
        dqcon_new = ql_new - ql_series[0] + qi_new - qi_series[0]
        if (qt_new < qt_min):
            if verbose:
                print 'warning: parcel now subsaturated'
            qt_new = qt_min
            qv_new = qt_min
            ql_new = 0.
            qi_new = 0.
        thetae_new = theta_e(p,T,qv_new,ql_new,qi_new)
        if verbose:
            print 'New ql: '+str(ql_new)
            print 'New qi: '+str(qi_new)
            print 'New thetae: '+str(thetae_new)
            print 'New qt: '+str(qt_new)
            print 'Fraction fallout: '+str(1. - dqcon_new/dqcon)
        return (ql_new,qi_new,thetae_new,qt_new) 

## Compute adiabat  #####################################################
def adiabat(T0,qv0,p_prof,ql0=0.,qi0=0.,fallout_factor=0.):
    thetae0 = theta_e(p_prof[0],T0,qv0,ql0,qi0)
    n = len(p_prof)
    T_prof       = np.zeros(n)
    qv_prof      = np.zeros(n)
    ql_prof      = np.zeros(n)
    qi_prof      = np.zeros(n)
    qt_prof      = np.zeros(n)
    rho_prof     = np.zeros(n)
    thetae_prof  = np.zeros(n)
    thetae_prof[0] = thetae0
    qt_prof[0]     = qv0 + ql0 + qi0
    T_prof[0],qv_prof[0],ql_prof[0],qi_prof[0],rho_prof[0] = \
    adiabatic_state(thetae0,qv0+ql0+qi0,p_prof[0])
    for i in range(1,n):
        T_prof[i],qv_prof[i],ql_prof[i],qi_prof[i],rho_prof[i] = \
        adiabatic_state(thetae_prof[i-1],qv_prof[i-1]+ql_prof[i-1]+\
                        qi_prof[i-1],p_prof[i])
        if (fallout_factor > 0.):
            ql_prof[i],qi_prof[i],thetae_prof[i],qt_prof[i] = \
            fallout(fallout_factor,p_prof[i],T_prof[i],\
            qv_prof[(i-1):(i+1)],ql_prof[(i-1):(i+1)],\
            qi_prof[(i-1):(i+1)])

        else:
            thetae_prof[i] = theta_e(p_prof[i],T_prof[i],qv_prof[i],\
                             ql_prof[i],qi_prof[i])
            qt_prof[i] = qv_prof[i] + ql_prof[i] + qi_prof[i]

        rho_prof[i] = adiabatic_state(thetae_prof[i],qt_prof[i],p_prof[i])[4]

    return (T_prof,qv_prof,ql_prof,qi_prof,qt_prof,rho_prof,thetae_prof)

## Compute CAPE  ########################################################
def calc_cape(b,z,dz=1.,LCL=-99.,LNB=-99.):
    
    # check that LCL and LNB, if supplied, are in range
    if (LCL != -99. and LCL < z[0]):
        print 'Error in calc_cape: LCL out of z range'
        return
    if (LNB != -99. and LCL > z[-1]):
        print 'Error in calc_cape: LNB out of z range'
        return
    
    #interpolate data to high res
    z_hr = np.arange(start=z[0],stop=z[-1],step=dz)
    b_hr = np.interp(z_hr,z,b)
    
    # find LCL and LNB
    if LNB != -99.:
        lnb_k = np.argmin(np.abs(z_hr-LNB))
        if b_hr[lnb_k] < 0.:
            print 'Error in calc_cape: negative buoyancy at supplied LNB'
            print b_hr[(lnb_k-2):(lnb_k+2)]
            return
    else:
        lnb_k = max(loc for loc, val in enumerate(np.sign(b_hr)) if val > 0.)
    if LCL != -99.:
        lcl_k = np.argmin(np.abs(z_hr-LCL))
        if b_hr[lcl_k] < 0.:
            print 'Error in calc_cape: negative buoyancy at supplied LCL'
            return
    else:
        lcl_k = min(loc for loc, val in enumerate(np.sign(b_hr[0:(lnb_k+1)])) if val > 0.)
        
    # compute CAPE with trapezoidal approximation
    cape = np.trapz(b_hr[lcl_k:(lnb_k+1)],z_hr[lcl_k:(lnb_k+1)])
        
    return (cape,b_hr[lcl_k:(lnb_k+1)],z_hr[lcl_k:(lnb_k+1)])

##  Moist static energy  ################################################
def moist_stat_en(z,T,qv):
	return gg*z + c_pd*T + L_v*qv
