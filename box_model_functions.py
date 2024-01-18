#!/bin/env python 

def co2_emissions(yr, escheme):

    from scipy.interpolate import interp1d
    import numpy as np

    ## historical emissions
    time = np.arange(1764, 2006, step=1)
    spinup_time = np.arange(1, 5001, step=1)
    ideal_emit_hist = np.array([0] * len(spinup_time))

    emit_hist = [0,0.003,0.003,0.003,0.003,0.003,0.003,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,
                 0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.006,0.006,0.006,0.006,0.006,0.006,0.007,
                 0.007,0.007,0.008,0.008,0.010,0.009,0.009,0.009,0.010,0.010,0.010,0.010,0.010,0.011,0.011,0.011,0.011,
                 0.012,0.013,0.014,0.014,0.014,0.014,0.014,0.015,0.016,0.016,0.017,0.017,0.018,0.018,0.018,0.024,0.023,
                 0.023,0.024,0.024,0.025,0.029,0.029,0.030,0.031,0.033,0.034,0.036,0.037,0.039,0.043,0.043,0.046,0.047,
                 0.050,0.054,0.054,0.057,0.059,0.069,0.071,0.076,0.077,0.078,0.083,0.091,0.095,0.097,0.104,0.112,0.119,
                 0.122,0.130,0.135,0.142,0.147,0.156,0.173,0.184,0.174,0.188,0.191,0.195,0.196,0.210,0.236,0.243,0.256,
                 0.272,0.275,0.277,0.281,0.295,0.327,0.327,0.356,0.372,0.374,0.370,0.383,0.406,0.419,0.440,0.465,0.507,
                 0.534,0.552,0.566,0.617,0.624,0.663,0.707,0.784,0.750,0.785,0.819,0.836,0.879,0.943,0.850,0.838,0.901,
                 0.955,0.936,0.806,0.932,0.803,0.845,0.970,0.963,0.975,0.983,1.062,1.065,1.145,1.053,0.940,0.847,0.893,
                 0.973,1.027,1.130,1.209,1.142,1.192,1.299,1.334,1.342,1.391,1.383,1.160,1.238,1.392,1.469,1.419,1.630,
                 1.768,1.796,1.841,1.865,2.043,2.178,2.270,2.330,2.462,2.577,2.594,2.700,2.848,3.008,3.145,3.305,3.411,
                 3.588,3.800,4.076,4.231,4.399,4.635,4.644,4.615,4.883,5.029,5.105,5.387,5.332,5.168,5.127,5.110,5.290,
                 5.444,5.610,5.753,5.964,6.089,6.144,6.235,6.118,6.124,6.242,6.372,6.510,6.619,6.588,6.569,6.735,6.896,
                 6.949,7.286,7.672,7.971]
    
    emit_hist_gcb = [0.72, 0.74, 0.75, 0.77, 0.79, 0.79, 0.81, 0.81, 0.82, 0.84,
                     0.81, 0.78, 0.78, 0.79, 0.79, 0.78, 0.79, 0.79, 0.79, 0.80,
                     0.85, 0.89, 0.92, 0.95, 0.95, 0.97, 0.99, 1.00, 1.02, 1.04,
                     1.08, 1.09, 1.12, 1.15, 1.16, 1.16, 1.17, 1.20, 1.23, 1.24,
                     1.30, 1.34, 1.36, 1.37, 1.39, 1.43, 1.45, 1.48, 1.51, 1.57,
                     1.63, 1.66, 1.69, 1.74, 1.77, 1.82, 1.85, 1.96, 1.94, 1.97,
                     2.01, 2.02, 2.04, 2.11, 2.02, 2.00, 2.05, 2.09, 2.09, 1.95,
                     2.12, 2.03, 2.09, 2.21, 2.22, 2.25, 2.23, 2.32, 2.32, 2.40,
                     2.32, 2.22, 2.14, 2.20, 2.29, 2.33, 2.44, 2.52, 2.46, 2.53,
                     2.69, 2.77, 2.77, 2.80, 2.82, 2.57, 2.68, 2.80, 2.86, 2.80,
                     3.12, 3.27, 3.33, 3.41, 3.49, 3.72, 3.91, 4.03, 4.15, 4.36,
                     4.35, 4.24, 4.27, 4.35, 4.45, 4.47, 4.56, 4.69, 4.88, 5.03,
                     5.32, 5.48, 5.68, 5.88, 5.82, 5.85, 6.09, 6.25, 6.36, 6.46,
                     6.44, 6.39, 6.36, 6.50, 6.82, 6.93, 7.03, 7.20, 7.41, 7.45,
                     7.54, 7.64, 7.50, 7.56, 7.71, 7.85, 8.04, 8.57, 8.11, 8.25,
                     8.34, 8.30, 8.59, 9.09, 9.28, 9.37, 9.73, 9.81, 10.02, 9.99,
                     10.43, 10.77, 10.87, 10.90, 11.05, 11.17, 10.93, 11.03, 11.19,
                     11.36, 10.73, 11.21]
    
    
    #         https://c4mip.net/zecmip-protocol
    bell_curve_1000PgC = [0.20873014, 0.25276203, 0.30488921, 0.3663328, 0.43844296,
                          0.52270172, 0.62072365, 0.73425378, 0.86516239, 1.01543611,
                          1.18716509, 1.38252556, 1.6037577, 1.8531385, 2.13294934,
                          2.44543847, 2.79277839, 3.17701853, 3.60003364, 4.06346858,
                          4.56868053, 5.11667948, 5.70806844, 6.34298476, 7.0210441,
                          7.74128883, 8.50214249, 9.30137222, 10.1360608, 11.0025899,
                          11.8966362, 12.8131814, 13.746537, 14.6903849, 15.6378333,
                          16.5814888, 17.5135425, 18.4258706, 19.3101466, 20.1579639,
                          20.9609659, 21.7109814, 22.400162, 23.0211173, 23.5670474,
                          24.0318658, 24.4103126, 24.6980536, 24.8917628, 24.9891865,
                          24.9891865, 24.8917628, 24.6980536, 24.4103126, 24.0318658,
                          23.5670474, 23.0211173, 22.400162, 21.7109814, 20.9609659,
                          20.1579639, 19.3101466, 18.4258706, 17.5135425, 16.5814888,
                          15.6378333, 14.6903849, 13.746537, 12.8131814, 11.8966362,
                          11.0025899, 10.1360608, 9.30137222, 8.50214249, 7.74128883,
                          7.0210441, 6.34298476, 5.70806844, 5.11667948, 4.56868053,
                          4.06346858, 3.60003364, 3.17701853, 2.79277839, 2.44543847,
                          2.13294934, 1.8531385, 1.6037577, 1.38252556, 1.18716509,
                          1.01543611, 0.86516239, 0.73425378, 0.62072365, 0.52270172,
                          0.43844296, 0.3663328, 0.30488921, 0.25276203, 0.20873014]
    

    if escheme == "rcp85":
        time2 = np.arange(2006, 2101, step=1)
        time = np.concatenate([time, time2])
        emit_future = [8.162,8.352,8.543,8.735,8.926,9.187,9.448,9.709,9.970,10.232,10.493,10.754,11.015,
        11.276,11.538,11.768,11.998,12.228,12.458,12.688,12.918,13.149,13.379,13.609,13.839,
        14.134,14.429,14.723,15.018,15.313,15.608,15.902,16.197,16.492,16.787,17.128,17.470,
        17.812,18.154,18.496,18.837,19.179,19.521,19.863,20.205,20.544,20.883,21.222,21.561,
        21.900,22.240,22.579,22.918,23.257,23.596,23.833,24.069,24.306,24.543,24.779,25.016,
        25.252,25.489,25.726,25.962,26.107,26.251,26.395,26.540,26.684,26.829,26.973,27.117,
        27.262,27.406,27.499,27.592,27.685,27.778,27.871,27.964,28.058,28.151,28.244,28.337,
        28.377,28.417,28.458,28.498,28.538,28.579,28.619,28.659,28.700,28.740]

        emit = np.concatenate([emit_hist, emit_future])
    
    elif escheme == "hist_gcb":
        time_future = np.arange(0, 1000, step=1)
        time = np.concatenate([spinup_time, time_future+spinup_time[-1]])

        emit_future = np.zeros(len(time_future))
        emit_future[0:172] = emit_hist_gcb
        emit = np.concatenate([ideal_emit_hist, emit_future])

    elif escheme == "pulse": # done
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+spinup_time[-1]])

        emit_future = np.zeros(len(time_future))
        emit_future[0] = 1000
        emit = np.concatenate([ideal_emit_hist, emit_future])
        
    elif escheme == "pulse500": # done
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+spinup_time[-1]])

        emit_future = np.zeros(len(time_future))
        emit_future[0] = 500
        emit = np.concatenate([ideal_emit_hist, emit_future])
        
    elif escheme == "flat10_zec": # done
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+spinup_time[-1]])

        emit_future = np.zeros(len(time_future))
        emit_future[0:100] = 10
        emit = np.concatenate([ideal_emit_hist, emit_future])
#         emit = emit_future
#         time = time_future
        
    elif escheme == "CDR-pi-pulse":
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+spinup_time[-1]])
        
        emit_future = np.zeros(len(time_future))
        emit_future[0] = 1000
        emit = np.concatenate([ideal_emit_hist, emit_future])

    elif escheme == "pulse-reversal": #done
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])

        emit_future = np.zeros(len(time_future))
        emit_future[0] = 100
        emit_future[100] = -100
        emit = np.concatenate([ideal_emit_hist, emit_future])
    
    elif escheme == "esm-bell-1000PgC":
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])

        emit_future = np.zeros(len(time_future))                
        emit_future[0:100] = bell_curve_1000PgC
        emit = np.concatenate([ideal_emit_hist, emit_future])

    elif escheme == "esm-restoration": #HAVE NOT EDITED THIS ONE
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])
        
        emit_future = np.zeros(len(time_future)) 
        emit_future[0:99] = np.diff(bell_curve_1000PgC)
        emit = np.concatenate([ideal_emit_hist, emit_future])
        
    elif escheme == "CDR-pi-pulse-500":
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])
        
        emit_future = np.zeros(len(time_future))
        emit_future[0] = 500
        emit = np.concatenate([ideal_emit_hist, emit_future])
        
    elif escheme == "pulse-reversal-500-100":
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])
        
        emit_future = np.zeros(len(time_future))
        emit_future[0] = 500
        emit_future[100] = -100
        emit = np.concatenate([ideal_emit_hist, emit_future])

    elif escheme == "pulse-reversal-500-500":
        time_future = np.arange(1, 1001, step=1)
        time = np.concatenate([spinup_time, time_future+len(spinup_time)])
        
        emit_future = np.zeros(len(time_future))
        emit_future[0] = 500
        emit_future[100] = -500
        emit = np.concatenate([ideal_emit_hist, emit_future])

    else:
        print('unspecified emissions scenario')

# GEMS commented
#     time = [-1e6, time, 1e6]
#     time = np.array([-1e6, time, 1e6])
    time = np.insert(time, 0, -1e6) #GEMS
    time = np.append(time, 1e6) #GEMS

# #     emit = [0, emit, emit[-1]]
    emit = np.insert(emit, 0, 0) #GEMS
    emit = np.append(emit, emit[-1]) #GEMS
# GEMS commented


    # FF=interp1(time,emit,yr);         
    #FF = interp1d(time, emit, yr)
#     FF = np.interp(yr, time, emit) #GEMS
    FF_fctn = interp1d(time, emit) #GEMS commented
    FF = FF_fctn(yr) #GEMS commented

#     return(FF)
    return(FF)

def calc_pco2(t, s, ta, c, phg):
    '''
    this function calculates the partial pressure of co2
    '''

    import numpy as np

    pt = 0e-3
    sit = 40.0e-3
    tk = 273.15 + t
    tk100 = tk / 100.0
    tk1002 = tk100**2
    invtk = 1.0 / tk
    dlogtk = np.log(tk)
    ### note this variable has to change names since "is" is inbuilt in python 
    # is    = 19.924*s./(1000.-1.005*s);
    iss = 19.924 * s / (1000. - 1.005 * s)
    # is2   =is.*is;
    iss2 = iss**2
    sqrtis = np.sqrt(iss)
    s2 = s**2
    sqrts = np.sqrt(s)
    s15 = s ** 1.5
    scl = s / 1.80655

    fflocal = (np.exp(-162.8301 + 218.2968 / tk100  + 
               90.9241 * np.log(tk100) - 1.47696 * tk1002 + 
               s * (.025695 - .025225 * tk100 + 
               0.0049867 * tk1002)))

    k0local = (np.exp(93.4517 / tk100 - 60.2409 + 
               23.3585 * np.log(tk100) + 
               s * (0.023517 - 0.023656 * tk100 + 
               0.0047036 * tk1002)))

    k1local = 10**((-1 * (3670.7 * invtk - 
              62.008 + 9.7944 * dlogtk - 
              0.0118 * s + 0.000116 * s2)))

    k2local = 10**(-1 * (1394.7 * invtk + 4.777 - 
              0.0184 * s + 0.000118 * s2))

    kblocal = np.exp((-8966.90 - 2890.53 * sqrts - 77.942 * s + 
                     1.728 * s15 - 0.0996 * s2) * invtk + 
                     (148.0248 + 137.1942 * sqrts + 1.62142 * s) + 
                     (-24.4344 - 25.085 * sqrts - 0.2474 * s) * 
                     dlogtk + 0.053105 *sqrts * tk)

    k1plocal = np.exp(-4576.752 * invtk + 115.525 - 
                      18.453 * dlogtk + 
                     (-106.736 * invtk + 0.69171) * sqrts + 
                     (-0.65643 * invtk - 0.01844) * s)

    k2plocal = np.exp(-8814.715 * invtk + 172.0883 - 
                      27.927 * dlogtk + 
                     (-160.340 * invtk + 1.3566) * sqrts +  
                     (0.37335 * invtk - 0.05778) * s)

    k3plocal = np.exp(-3070.75 * invtk - 18.141 + 
                      (17.27039 * invtk + 2.81197) *
                      sqrts + (-44.99486 * invtk - 0.09984) * s)

    ksilocal = np.exp(-8904.2 * invtk + 117.385 - 
                      19.334 * dlogtk + 
                      (-458.79 * invtk + 3.5913) * sqrtis + 
                      (188.74 * invtk - 1.5998) * iss + 
                      (-12.1652 * invtk + 0.07871) * iss2 +  
                      np.log(1.0 - 0.001005 * s))

    kwlocal = np.exp(-13847.26 * invtk + 148.9652 -
                     23.6521 * dlogtk +
                    (118.67 * invtk - 5.977 + 1.0495 * dlogtk) * 
                    sqrts - 0.01615 * s)

    kslocal = np.exp(-4276.1 * invtk + 141.328 - 
                     23.093 * dlogtk + 
                     (-13856 * invtk + 324.57 - 47.986 * dlogtk) *sqrtis + 
                     (35474 * invtk - 771.54 + 114.723 * dlogtk) *iss - 
                     2698 * invtk * iss**1.5 + 1776 * invtk * iss2 + 
                     np.log(1.0 - 0.001005 * s))

    kflocal = np.exp(1590.2 * invtk - 12.641 + 1.525 * sqrtis + 
                     np.log(1.0 - 0.001005 * s) +
                     np.log(1.0 + (0.1400 / 96.062) * (scl) / kslocal))

    btlocal = 0.000232 * scl/10.811
    stlocal = 0.14 * scl/96.062
    ftlocal = 0.000067 * scl/18.998

    pHlocal = phg
    permil =1.0 / 1024.5
    pt = pt * permil
    sit = sit * permil
    ta = ta * permil
    c = c * permil

    ####################
    ## start iteration ##
    ####################

    phguess = pHlocal
    hguess = 10.0**(-phguess)
    bohg = btlocal*kblocal / (hguess + kblocal)
    stuff = (hguess * hguess * hguess 
             + (k1plocal * hguess * hguess)
             + (k1plocal * k2plocal * hguess) 
             + (k1plocal * k2plocal * k3plocal))
    h3po4g = (pt * hguess * hguess * hguess) / stuff
    h2po4g = (pt * k1plocal * hguess * hguess) / stuff
    hpo4g = (pt * k1plocal * k2plocal * hguess) / stuff
    po4g = (pt * k1plocal * k2plocal * k3plocal) / stuff

    siooh3g = sit * ksilocal / (ksilocal + hguess);

    cag = (ta - bohg - (kwlocal / hguess) + hguess 
           - hpo4g - 2.0*po4g + h3po4g - siooh3g)

    gamm  = c / cag
    hnew = (0.5 * (-k1local * (1 - gamm) + np.sqrt((k1local**2) * (1 - gamm)**2 
            +4 * k1local * k2local * (2 * gamm - 1) ) ))

    pHlocal_new = -np.log10(hnew)
    pHlocal = pHlocal_new

    pco2local = (c / fflocal / (1.0 + (k1local / hnew) +  
                 (k1local * k2local / (hnew**2))))
    fflocal = fflocal / permil

    return(pco2local, pHlocal, fflocal)

def get_matrix_index(arr_row_num, arr_col_num, row_ind, col_ind):
    import numpy as np
    pool_indices = []
    element_nums = np.arange(0, 9*5).reshape(arr_col_num, arr_row_num).transpose()
    for ind in range(0, len(row_ind)):
        # print(element_nums[row_ind[ind], col_ind[0][ind]])
        pool_indices.append(element_nums[row_ind[ind], col_ind[0][ind]])
    return(pool_indices)

def carbon_climate_derivs(t, y, PE, PS, PL, PO):
    '''
    this is the main function for the box model
    '''

    import numpy as np
    from scipy.interpolate import interp1d
    #import seawater as sw
    # added the necessary seawater functions to their own .py module
    from seawater_functions import dens0, dens, seck, T68conv

    Tloc = y[PE['Jtmp']].transpose()
    Nloc = y[PE['Jnut']].transpose()
    Dloc = y[PE['Jcoc']].transpose()
    Cloc = y[PE['Jcla']]
    patm = y[PE['Jatm']]


    ## special cases for ocean carbon pumps

    # homogenize T,S if no solubility pump (for pCO2 only)
    ############################### NOTE: Need to add T from whatever dict it's coming from ####################################
    if PS['DoOcn'] == 1:
        Tsol = PO['T']
        Ssol = PO['S']
        if PS['DoOcnSol'] == 0: 
            Ttmp=Tsol.flatten()
            Stmp=Ssol.flatten()
            Tsol[0,PO['Isfc']] = np.sum(Ttmp[PO['Isfc']] * PO['A'][PO['Isfc']]) / np.sum(PO['A'][PO['Isfc']])
            Ssol[0,PO['Isfc']] = np.sum(Stmp[PO['Isfc']] * PO['A'][PO['Isfc']]) / np.sum(PO['A'][PO['Isfc']])

        # homogenize alkalinity if no bio pump
        TAsol = PO['TA']
        if PS['DoOcnBio'] == 0:
            TAsol[PO['Isfc']] = np.sum(PO['TA'][PO['Isfc']] * PO['A'][PO['Isfc']]) / np.sum(PO['A'][PO['Isfc']])

    ## update basic quantities

    # time
    ymod = t / PE['spery'] # year in model time (starting from 0)
    ycal = ymod - PS['yspin'] + PS['ypert'] # calendar year (negative means "BCE")
    if ycal < PS['ypert']:
        doAtm = 0 # hold atmospheric co2 constant to equilibrate
    else:
        doAtm = 1 # allow atmospheric co2 to evolve

    # interp1d example
    # matlab: interp1(x, y, xn, 'linear')
    # python: yn_f2 = interp1d(x[::-1], y[::-1])
    # python: yn_py2 = yn_f2(xn)

    # atmosphere + climate
    FF = co2_emissions(ycal, PS['escheme']) # fossil fuel co2 emissions (Pg/yr)	
    # [ycal FF]
    FF = FF * 1e15 / 12 / PE['spery'] # convert to molC/s
    RFco2 = 5.35 * np.log(patm / PE['patm0']) * PS['DoRadCO2'] # radiative forcing from CO2
    RFsto=np.interp(round(ycal),PS['Yint'].transpose(), PS['RFint'].transpose())
    RF = (RFco2 + np.nansum(RFsto)) * doAtm
    dTbar = np.sum(Tloc[PO['Isfc']] * PO['A'][PO['Isfc']]) / np.sum(PO['A'][PO['Isfc']])

    #------ terrestrial
    NPPfac = 1 + np.interp(ycal,PS['Yint'].transpose(), PS['NPPint'].transpose())
    
    NPP = PL['NPP_o'] * NPPfac * (1 + PS['CCC_LC'] * PL['beta_fert'] * np.log(patm / PE['patm0'])) # perturbation NPP
    #krate = np.diag(PL['kbase']) * PL['Q10_resp']**(PS['CCC_LT'] * dTbar / 10)  # scaled turnover rate
    krate = PL['kbase'] * PL['Q10_resp']**(PS['CCC_LT'] * dTbar / 10)  # scaled turnover rate (vector)
    ## create a matrix version of krate with values on the diagonal 
    krate_diag = np.zeros((krate.shape[0], krate.shape[0]))
    krate_diag_row, krate_diag_col = np.diag_indices(krate_diag.shape[0])
    krate_diag[krate_diag_row, krate_diag_col] = np.squeeze(krate) # matrix version
    Rh = np.sum(-np.sum(PL['acoef'],0) * np.transpose(krate) * Cloc) # Heterotrophic respiration
    
    # To get back to PgC for land pools we take Cloc*(land area)*12e-15. This means that Cloc is in mol/km2
    NEE = (NPP - Rh) * PL['Ala'] # total carbon pool tendency (mol/s)
    
    # set fluxes to 0 in ocean-only case
    if PS['DoTer'] == 0:
        NEE = 0
        krate = 0
        NPP = 0
        Rh = 0

    #------ ocean
    if PS['DoOcn'] == 1:
        Qbio = PO['Qup'] + PO['Qrem']
        pco2loc, pHloc, Ksol = calc_pco2(Tsol + PS['CCC_OT'] * Tloc, Ssol, TAsol, Dloc, PO['pH0']) # CO2 chemistry
        pco2Cor = patm * PS['CCC_OC'] + PE['patm0'] * (1 - PS['CCC_OC']) # switch for ocean carbon-carbon coupling
        Fgasx = PO['kwi'] * PO['A'] * Ksol * (pco2loc - pco2Cor) # gas exchange rate

        # circulation change
        #rho = sw.dens(PO['S'], PO['T'] + Tloc, PO['T'] * 0).flatten() # density
        rho = dens(PO['S'], PO['T'] + Tloc, PO['T'] * 0).flatten() # density
        bbar = PO['rho_o'][6] - PO['rho_o'][2]
        db = (rho[6] - rho[2]) - bbar
        Psi = PO['Psi_o'] * (1 - PS['CCC_OT'] * PO['dPsidb'] * db / bbar)
        
        #------ Compute Tendencies - should have units mol/s
        dNdt = np.matmul(Psi + Qbio, Nloc.transpose()) ######!!!! There is likely a problem with dNdt - need to check with matlab
        dDdt = np.matmul(Psi, Dloc.transpose()) + PO['Rcp'] * np.matmul(Qbio, Nloc.transpose()) - Fgasx / PO['V'].transpose()
        
        
    # set fluxes to 0 in land-only case
    if PS['DoOcn'] == 0:
        Fgasx = 0
        Psi = PO['Psi_o'] # this probably gets set somewhere else when the ocn is turned on.. check

    # [ycal/yend]

    #------ Compute Tendencies - should have units mol/s
    dTdt = np.matmul(Psi,Tloc.transpose()) -((PO['lammbda'] / PO['V']) * Tloc).transpose() + RF / PO['cm'].transpose() ###!!! problem here too?
    
    dAdt = (1 / PE['ma']) * (np.sum(Fgasx) - NEE + FF) 

    # land tendencies
    dCdt = np.matmul(np.matmul(PL['acoef'],krate_diag),  Cloc.reshape(9, 1)) + NPP * PL['bcoef']
    
    ## matrix of derivatives

    dydtmat = np.copy(PE['m0']) #initialize with a matrix of zeros. Making a copy here to avoid overwriting the values in PE
    if PS['DoOcn'] == 1: 
        dydtmat[0:PE['nb'],1] = dNdt.flatten()
        dydtmat[0:PE['nb'],2] = dDdt.flatten()
        
    dydtmat[0:PE['nb'],0] = dTdt.flatten()
    dydtmat[0, 4] = dAdt * doAtm;
    
    if PS['DoTer'] == 1:
        dydtmat[0:PE['np'],3] = dCdt.flatten()

    temporary = np.transpose(dydtmat).flatten()
    dydt=temporary[PE['Ires']]


    return(dydt)

