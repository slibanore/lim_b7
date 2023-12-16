from EBL_forecast import *

rho_crit = (2.77536627e11*(u.Msun.to(u.kg)*u.kg*u.Mpc**-3) * cu.c**2 ).to(u.erg/u.Mpc**3)

f_gg = 1 # good for axions
F_g = 0 # sets stimulated decay due to CMB - check different values
f_DM = 0.1 # DM fraction 


def signal_DM(nu_obs, m_DM, decay_DM):

    m_DM *= u.eV
    decay_DM *= u.s**-1

    z = -1 + m_DM / (4*np.pi*cu.hbar) / nu_obs  
    if zmin_gal < z < zmax_gal:
        # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
        eps_DM = (f_gg * f_DM * camb_pars.omegac * rho_crit  * light**2 * decay_DM * (1+F_g)/nu_obs).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
    else:
        eps_DM = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3

    return eps_DM


def dJdz_DM(z, detector, m_DM, decay_DM, run = False,\
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = ''):

    if run:
        if detector == 'GALEX_NUV':
            nu_min = nu_min_gNUV
            nu_max = nu_max_gNUV

        elif detector == 'GALEX_FUV':
            nu_min = nu_min_gFUV
            nu_max = nu_max_gFUV

        elif detector == 'ULTRASAT':
            nu_min = nu_min_US
            nu_max = nu_max_US

        else:
            print('Check detector!')
            return 


        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, decay_DM) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)

        unit = (Response(lambda_from_nu(1.*u.Hz), detector) *sum_signal(1*u.Hz).unit * np.exp(-tau_Lya(lambda_from_nu(1*u.Hz),0.))/(1*u.Hz)).unit


        intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        rest_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.AA).value
            intg_arr[i] = intg(nu_obs_arr[i])

        dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * np.trapz(intg_arr,nu_obs_arr)*(unit*u.Hz/u.steradian) 

    else:
        zval, dJdzval = np.genfromtxt(filename)
        dJdz = interp1d(zval,dJdzval)(z) * u.Jy/u.steradian  

    return (dJdz.to(u.Jy/u.steradian)).value


def bJ_z_DM(z, detector,  m_DM, decay_DM,run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = ''):

    if run:
        if detector == 'GALEX_NUV':
            nu_min = nu_min_gNUV
            nu_max = nu_max_gFUV

        elif detector == 'GALEX_FUV':
            nu_min = nu_min_gFUV
            nu_max = nu_max_gFUV

        elif detector == 'ULTRASAT':
            nu_min = nu_min_US
            nu_max = nu_max_US

        else:
            print('Check detector!')
            return 
#
        
        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, decay_DM) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)

        sum_signal_bias = lambda nu_obs: signal_DM(nu_obs, m_DM, decay_DM) + bias(nu_obs*(1+z),z,val_bias)*signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900) # !! check

        intg_num = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal_bias(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
#
        intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        rest_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        intg_den_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.AA).value
            intg_arr[i] = intg_num(nu_obs_arr[i])
            intg_den_arr[i] = intg_den(nu_obs_arr[i])

        num = np.trapz(intg_arr,nu_obs_arr)
        den = np.trapz(intg_den_arr,nu_obs_arr)

        bJ = num / den

    else:
        zval, bJval = np.genfromtxt(filename)
        bJ = interp1d(zval,bJval)(z) 

    return bJ


def compute_vals_DM(m_DM, decay_DM, reduced_z = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False):

    reduced_label = '_reduced' if reduced_z else ''

    dir_DM = './results/EBL/DM_' +str(f_DM) + '_' + str(m_DM) + ',' + str(np.log10(decay_DM))+'/'

    if not os.path.exists(dir_DM):
        os.makedirs(dir_DM)

    with Pool(6) as pool:
        
        print('Doing NUV')
        dJdz_NUV = partial(dJdz_DM,detector='GALEX_NUV',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900, filename=dir_DM+'dJdz_GALEX_NUV' + reduced_label + '.dat')
        if reduced_z:
            dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
            np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_nuv)))
        else:
            dJ_nuv = pool.map(dJdz_NUV, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_nuv)))

        print('Doing FUV')
        dJdz_FUV = partial(dJdz_DM,detector='GALEX_FUV',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
            np.savetxt(dir_DM+'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_fuv)))
        else:
            dJ_fuv = pool.map(dJdz_FUV, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_fuv)))
    
        print('Doing b NUV')
        bJ_nuv_f =  partial(bJ_z_DM,detector='GALEX_NUV',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
            np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_nuv)))
        else:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_nuv)))

        print('Doing b FUV')
        bJ_fuv_f =  partial(bJ_z_DM,detector='GALEX_FUV',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)
            np.savetxt(dir_DM+  'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_fuv)))
        else:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_fuv)))
     
        print('Doing ULTRASAT')
        dJdz_f = partial(dJdz_DM,detector='ULTRASAT',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat')

        if reduced_z:
            dJ_U = pool.map(dJdz_f, z_gals_interp)
            np.savetxt(dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_U)))
        else:
            dJ_U = pool.map(dJdz_f, z_gals)
            np.savetxt(dir_DM+ 'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(dJ_U)))

        print('Doing b ULTRASAT')
        bJ_f =  partial(bJ_z_DM,detector='ULTRASAT',m_DM=m_DM, decay_DM=decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat')
        if reduced_z:
            bJ_U = pool.map(bJ_f, z_gals_interp)        
            np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_U)))
        else:
            bJ_U = pool.map(bJ_f, z_gals)        
            np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(bJ_U)))
    return




def plot_collective_DM(m_DM, decay_DM):

    fig, ax = plt.subplots(2,1, figsize=(9,11))

    dir_base = './results/EBL/'
    dir_DM = './results/EBL/DM_' +str(f_DM) + '_'+ str(m_DM) + ',' + str(np.log10(decay_DM))+'/'

    fontsize_label = 23
    fontsize_legend = 21
    fontsize_tick = 19
###########################
    wavelenght, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T

    #plt.subplot(311)
    z = [0.,0.5,1.,2]

    rest_wave = np.linspace(700,3500,400)
    for zv in range(len(z)):
        s = np.zeros(len(rest_wave))
        s_tot = np.zeros(len(rest_wave))
        wave = np.zeros((len(rest_wave)))
        for i in range(len(rest_wave)):
            wave[i] = rest_wave[i] * (1.+z[zv])
            if z[zv] <= 1:
                detector = 'GALEX_NUV' # same for galex_fuv in EW
            else:
                detector = 'ULTRASAT' # same for hetdex and spherex
            s[i] = signal(rest_wave[i]*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value

            s_tot[i] = signal(rest_wave[i]*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value + signal_DM(nu_from_lambda(rest_wave[i]*u.AA)/(1+z[zv]), m_DM, decay_DM).value

        s_1500 = signal(1500*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value + signal_DM(nu_from_lambda(1500*u.AA)/(1+z[zv]), m_DM, decay_DM).value

        ax[0].plot(wave, s, color= colors[zv],  zorder= len(z)-zv,linestyle='--')
        ax[0].plot(wave, s_tot, color= colors[zv], label=r'$z = %g$'%z[zv], zorder= len(z)-zv)
        ax[0].plot(1500*(1+z[zv]), s_1500, marker = 'o',color= colors[zv], zorder= len(z)-zv,markersize=7)

    ax[0].set_yscale('log')
    ax[0].set_xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize_label)
    ax[0].set_ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_label)
    ax[0].legend(loc=2,ncol=2,fontsize=fontsize_legend)
    ax[0].set_xlim(700,4900)
    ax[0].set_ylim(1e25,5e28)

    ax[0].set_title(r'$\rm UV\,-\,EBL$',fontsize=fontsize_label,y=1.01)

    ax[0].tick_params(axis='y', labelsize=fontsize_tick) 
    ax[0].tick_params(axis='x', labelsize=fontsize_tick) 
###########################


########################

    dJ_U_tot =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM + 'dJdz_ULTRASAT_reduced.dat')

    dJ_N_tot =  dJdz_DM(z_gals('SDSS'),detector='GALEX_NUV',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM +'dJdz_GALEX_NUV_reduced.dat')

    dJ_F_tot =  dJdz_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U_tot = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM + 'bJ_ULTRASAT_reduced.dat')

    bJ_N_tot = bJ_z_DM(z_gals('SDSS'),detector='GALEX_NUV',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM+ 'bJ_GALEX_NUV_reduced.dat')

    bJ_F_tot = bJ_z_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM,decay_DM=decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM + 'bJ_GALEX_FUV_reduced.dat')

    dJ_U =  dJdz(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base + 'dJdz_ULTRASAT_reduced.dat')

    dJ_N =  dJdz(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base +'dJdz_GALEX_NUV_reduced.dat')

    dJ_F =  dJdz(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U = bJ_z(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base + 'bJ_ULTRASAT_reduced.dat')

    bJ_N = bJ_z(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base+ 'bJ_GALEX_NUV_reduced.dat')

    bJ_F = bJ_z(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base + 'bJ_GALEX_FUV_reduced.dat')

    window_size = 70  # Adjust this value based on your preference
    bJU_smoothed = moving_average(dJ_U*bJ_U, window_size)
    bJGN_smoothed = moving_average(dJ_N*bJ_N, window_size)
    bJGF_smoothed = moving_average(dJ_F*bJ_F, window_size)
# 
    bJU_tot_smoothed = moving_average(dJ_U_tot*bJ_U_tot, window_size)
    bJGN_tot_smoothed = moving_average(dJ_N_tot*bJ_N_tot, window_size)
    bJGF_tot_smoothed = moving_average(dJ_F_tot*bJ_F_tot, window_size)
# 
    ax[1].plot(z_gals('SDSS')[:len(bJGF_smoothed)],bJGF_smoothed,color=color_FUV,linestyle='--')
    ax[1].plot(z_gals('SDSS')[:len(bJGN_smoothed)],bJGN_smoothed,color=color_NUV,linestyle='--')
    ax[1].plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_smoothed,color=color_ULTRASAT,linestyle='--')

    ax[1].plot(z_gals('SDSS')[:len(bJGF_smoothed)],bJGF_tot_smoothed,label=r'$\rm GALEX\, FUV$',color=color_FUV)
    ax[1].plot(z_gals('SDSS')[:len(bJGN_smoothed)],bJGN_tot_smoothed,label=r'$\rm GALEX\, NUV$',color=color_NUV)
    ax[1].plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_tot_smoothed,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)


    ax[1].set_xlabel(r'$z$',fontsize=fontsize_label)
    ax[1].set_ylabel(r'$b_JdJ_{\nu_{\rm obs}}/dz\,[{\rm Jy/sr}]$',fontsize=fontsize_label)
    ax[1].legend(fontsize=fontsize_legend,loc=1)

    ax[1].tick_params(axis='y', labelsize=fontsize_tick) 
    ax[1].tick_params(axis='x', labelsize=fontsize_tick) 

    ax[1].set_title(r'$\rm CBR\,Reconstruction$',fontsize=fontsize_label,y=1.01)

    ax[1].set_xlim(0,2.3)
###############################

    ax[0].yaxis.set_label_coords(-.07, 0.5)
    ax[1].yaxis.set_label_coords(-.07, 0.5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.37)
    #plt.savefig('results/PLOTS/EBL/collective_plot.png',bbox_inches='tight')

    plt.show()

    return 
