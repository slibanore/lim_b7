from EBL_forecast import *
import EBL_plots_paper as p 

rho_crit = (3*H(0)**2/(8*np.pi*cu.G)).to(u.Msun*u.Mpc**-3)


f_gg = 1 # good for axions
F_g = 0 # sets stimulated decay due to CMB - check different values

logaddexp = lambda x,y: x + np.log(1 + np.exp(y - x)) if (x > y) else y + np.log(1 + np.exp(x - y))

sigmoid = lambda lambda_rest: np.exp(-logaddexp(0, -(lambda_rest.value - (1216*u.AA).value) / 100))


def signal_DM(nu_obs, m_DM, f_DM_decay_DM):

    m_DM *= u.eV
    f_DM_decay_DM *= u.s**-1

    nu_rest_DM = (m_DM / (4*np.pi*cu.hbar)).to(u.Hz)

    z = -1 + nu_rest_DM/ nu_obs  
    if zmin_gal <= z <= zmax_gal:
        # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
        eps_DM = fLyC(z,False) * sigmoid(lambda_from_nu(nu_rest_DM)) * (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
    else:
        eps_DM = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3

    return eps_DM


def signal_DM_zdep(z, m_DM, f_DM_decay_DM):

    m_DM *= u.eV
    f_DM_decay_DM *= u.s**-1

    nu_rest_DM = m_DM / (4*np.pi*cu.hbar) #/ (1+z)

    eps_DM = fLyC(z,False) * sigmoid(lambda_from_nu(nu_rest_DM)) * (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions

    return eps_DM


def plot_model_DM(m_DM=10,f_DM_decay = 1e-23):

    plt.figure(figsize=(9,6))

    fontsize_label = 23
    fontsize_legend = 21
    fontsize_tick = 19
###########################
    wavelenght, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T

    #plt.subplot(311)
    z = [1.]

    rest_wave = np.linspace(700,3500,400)
    for zv in range(len(z)):

        lambda_obs_DM = lambda_from_nu(m_DM*u.eV / (4*np.pi*cu.hbar) / (1+z[zv])).value

        s = np.zeros(len(rest_wave)+1)
        s_tot = np.zeros(len(rest_wave)+1)
        wave = np.zeros((len(rest_wave)+1))
        unsort_wave = np.zeros((len(rest_wave)+1))
        for i in range(len(rest_wave)):

            unsort_wave[i] = rest_wave[i] * (1.+z[zv])

        unsort_wave[-1] = lambda_obs_DM

        wave = np.sort(unsort_wave)

        for i in range(len(wave)):

            detector = 'GALEX_NUV' # same for galex_fuv in EW

            s[i] = signal(wave[i]/(1.+z[zv])*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value

            if wave[i] == lambda_obs_DM:
                use_signal_DM = signal_DM_zdep(z[zv], m_DM, f_DM_decay).value
            else:
                use_signal_DM = 0.

            s_tot[i] = signal(wave[i]/(1+z[zv])*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value + use_signal_DM

        s_1500 = signal(1500*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value 

        if z[zv] == 1.:
            plt.plot(wave, s_tot, color= 'k')
            plt.plot(1500*(1+z[zv]), s_1500, marker = 'o',color= 'k', zorder= len(z)-zv,markersize=7)
   
    plt.yscale('log')
    plt.ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_legend)
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize_label)
    plt.xlim(700,5900)
    plt.ylim(5e24,1e28)

    plt.text(1500*(1+1), 3e26, r'$\alpha_{1500}(z)$', rotation=-1.5, fontsize=fontsize_legend, ha='center', va='center')
    plt.annotate('', xy=(0.445, 0.45), xytext=(0.445, 0.01),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='<->'), ha='center', va='center')

    plt.text(2455*(1+1), 6e25, r'$m_{\rm DM}c^2$', rotation=-1.5, fontsize=fontsize_legend, ha='center', va='center')
    plt.annotate('', xy=(0.817, 0.443), xytext=(0.817, 0.35),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='->'), ha='center', va='center')

    plt.text(1570*(1+1), 3e25, r'$\epsilon_{1500}(z)$', rotation=90, fontsize=fontsize_legend, ha='center', va='center')

    plt.text(1030*(1+1), 1.8e26, r'$\alpha_{1100}(z)$', rotation=25, fontsize=fontsize_legend, ha='center', va='center')

    plt.text(750*(1+1), 2.6e25, r'$\alpha_{900}(z)$', rotation=18, fontsize=fontsize_legend, ha='center', va='center')

    plt.text(990*(1+1), 2.5e25, r'$f_{\rm LyC}(z)$', rotation=90, fontsize=fontsize_legend, ha='center', va='center')

    plt.text(1216*(1+1), 5e27, r'$\rm EW(z)$', rotation=0, fontsize=fontsize_legend, ha='center', va='center')


    plt.xticks([]) 
    plt.yticks([]) 
    plt.title(r'$\rm Comoving\,Volume\,Emissivity$',fontsize=fontsize_label,y=1.01)
#    plt.set_title(r'$\rm Comoving\,volume\,emissivity$',fontsize=fontsize_label, y=.8, x=.71)

#    plt.xaxis.set_label_coords(0.5, -.06)
#    plt.yaxis.set_label_coords(-.07, 0.5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.37)
    plt.savefig('results/PLOTS/EBL_DM/collective_plot.png',bbox_inches='tight')

    plt.show()

    return 

def dJdz_DM(z, detector, m_DM, f_DM_decay_DM, run = False,\
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


        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)

        unit = (Response(lambda_from_nu(1.*u.Hz), detector) *sum_signal(1*u.Hz).unit * np.exp(-tau_Lya(lambda_from_nu(1*u.Hz),0.))/(1*u.Hz)).unit


        intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        obs_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(obs_wave))
        intg_arr = np.zeros(len(obs_wave))
        for i in range(len(obs_wave)):
            nu_obs_arr[i] = nu_from_lambda(obs_wave[i]*u.AA).value
            intg_arr[i] = intg(nu_obs_arr[i])

        dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * np.trapz(intg_arr,nu_obs_arr)*(unit*u.Hz/u.steradian) 

    else:
        zval, dJdzval = np.genfromtxt(filename)
        dJdz = interp1d(zval,dJdzval)(z) * u.Jy/u.steradian  

    return (dJdz.to(u.Jy/u.steradian)).value


def bJ_z_DM(z, detector,  m_DM, f_DM_decay_DM,run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = ''):

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
        
        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)

        sum_signal_bias = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM) + bias(nu_obs*(1+z),z,val_bias)*signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900) # !! check

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


def compute_vals_DM(m_DM, f_DM_decay_DM, reduced_z = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False):

    reduced_label = '_reduced' if reduced_z else ''

    dir_DM = './results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    if not os.path.exists(dir_DM):
        os.makedirs(dir_DM)

    with Pool(6) as pool:
        
        print('Doing NUV')
        dJdz_NUV = partial(dJdz_DM,detector='GALEX_NUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900, filename=dir_DM+'dJdz_GALEX_NUV' + reduced_label + '.dat')
        if reduced_z:
            dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
            np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_nuv)))
        else:
            dJ_nuv = pool.map(dJdz_NUV, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_nuv)))

        print('Doing FUV')
        dJdz_FUV = partial(dJdz_DM,detector='GALEX_FUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
            np.savetxt(dir_DM+'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_fuv)))
        else:
            dJ_fuv = pool.map(dJdz_FUV, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_fuv)))
    
        print('Doing b NUV')
        bJ_nuv_f =  partial(bJ_z_DM,detector='GALEX_NUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
            np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_nuv)))
        else:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_nuv)))

        print('Doing b FUV')
        bJ_fuv_f =  partial(bJ_z_DM,detector='GALEX_FUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)
            np.savetxt(dir_DM+  'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_fuv)))
        else:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals('SDSS'))
            np.savetxt(dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_fuv)))
     
        print('Doing ULTRASAT')
        dJdz_f = partial(dJdz_DM,detector='ULTRASAT',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat')

        if reduced_z:
            dJ_U = pool.map(dJdz_f, z_gals_interp)
            np.savetxt(dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_U)))
        else:
            dJ_U = pool.map(dJdz_f, z_gals)
            np.savetxt(dir_DM+ 'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(dJ_U)))

        print('Doing b ULTRASAT')
        bJ_f =  partial(bJ_z_DM,detector='ULTRASAT',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat')
        if reduced_z:
            bJ_U = pool.map(bJ_f, z_gals_interp)        
            np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_U)))
        else:
            bJ_U = pool.map(bJ_f, z_gals)        
            np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(bJ_U)))
    return




def plot_collective_DM():

    m_DM_1 = 10
    m_DM_2 = 20
    f_DM_decay_DM = 1e-23
    f_DM_decay_DM_low = 1e-25
    f_DM_decay_DM_high = 1e-21

    plt.figure(figsize=(10,8))

    dir_base = './results/EBL/'
    dir_DM_1 = './results/EBL/DM_' + str(m_DM_1) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    dir_DM_2 = './results/EBL/DM_' + str(m_DM_2) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    dir_DM_low = './results/EBL/DM_' + str(m_DM_1) + ',' + str(np.log10(f_DM_decay_DM_low))+'/'
    dir_DM_high = './results/EBL/DM_' + str(m_DM_1) + ',' + str(np.log10(f_DM_decay_DM_high))+'/'

    fontsize_label = 23
    fontsize_legend = 21
    fontsize_tick = 19
###########################
########################
    dJ_U_tot =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_1 + 'dJdz_ULTRASAT_reduced.dat')

    dJ_F_tot =  dJdz_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_1 + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U_tot = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_1 + 'bJ_ULTRASAT_reduced.dat')

    bJ_F_tot = bJ_z_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_1 + 'bJ_GALEX_FUV_reduced.dat')

    dJ_U =  dJdz(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base + 'dJdz_ULTRASAT_reduced.dat')

    dJ_F =  dJdz(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U = bJ_z(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base + 'bJ_ULTRASAT_reduced.dat')

    bJ_F = bJ_z(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base + 'bJ_GALEX_FUV_reduced.dat')


    dJ_U_tot_m2 =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_2 + 'dJdz_ULTRASAT_reduced.dat')

    dJ_F_tot_m2 =  dJdz_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_2 + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U_tot_m2 = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_2 + 'bJ_ULTRASAT_reduced.dat')

    bJ_F_tot_m2 = bJ_z_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_2 + 'bJ_GALEX_FUV_reduced.dat')


    dJ_U_tot_low =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_low + 'dJdz_ULTRASAT_reduced.dat')

    dJ_F_tot_low =  dJdz_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_low + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U_tot_low = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_low + 'bJ_ULTRASAT_reduced.dat')

    bJ_F_tot_low = bJ_z_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_low + 'bJ_GALEX_FUV_reduced.dat')


    dJ_U_tot_high =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_high + 'dJdz_ULTRASAT_reduced.dat')

    dJ_F_tot_high =  dJdz_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_high + 'dJdz_GALEX_FUV_reduced.dat')

    bJ_U_tot_high = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_high + 'bJ_ULTRASAT_reduced.dat')

    bJ_F_tot_high = bJ_z_DM(z_gals('SDSS'),detector='GALEX_FUV',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_high + 'bJ_GALEX_FUV_reduced.dat')



    window_size = 70  # Adjust this value based on your preference
    bJU_smoothed = moving_average(dJ_U*bJ_U, window_size)
    bJGF_smoothed = moving_average(dJ_F*bJ_F, window_size)
# 
    bJU_tot_smoothed = moving_average(dJ_U_tot*bJ_U_tot, window_size)
    bJGF_tot_smoothed = moving_average(dJ_F_tot*bJ_F_tot, window_size)
# 
    bJU_tot_smoothed_low = moving_average(dJ_U_tot_low*bJ_U_tot_low, window_size)
    bJGF_tot_smoothed_low = moving_average(dJ_F_tot_low*bJ_F_tot_low, window_size)

    bJU_tot_smoothed_2 = moving_average(dJ_U_tot_m2*bJ_U_tot_m2, window_size)
    bJGF_tot_smoothed_2 = moving_average(dJ_F_tot_m2*bJ_F_tot_m2, window_size)

    bJU_tot_smoothed_high = moving_average(dJ_U_tot_high*bJ_U_tot_high, window_size)
    bJGF_tot_smoothed_high = moving_average(dJ_F_tot_high*bJ_F_tot_high, window_size)

    #plt.plot(z_gals('SDSS')[:len(bJGF_smoothed)],bJGF_smoothed,color=color_FUV,label=r'$\rm GALEX\, FUV$')
    plt.plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_smoothed,color=color_ULTRASAT,label=r'$\rm UV-EBL$')

    #plt.plot(z_gals('SDSS')[:len(bJGF_tot_smoothed)],bJGF_tot_smoothed,color=color_FUV,linestyle='--')
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed)],bJU_tot_smoothed,color=color_ULTRASAT,linestyle='--',label=r'$f_{\rm DM}\Gamma = 10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM))

    #plt.plot(z_gals('SDSS')[:len(bJGF_tot_smoothed_low)],bJGF_tot_smoothed_low,color=color_FUV,linestyle=':')
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_low)],bJU_tot_smoothed_low,color=color_ULTRASAT,linestyle=':',label=r'$f_{\rm DM}\Gamma = 10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM_low))

    #plt.plot(z_gals('SDSS')[:len(bJGF_tot_smoothed_high)],bJGF_tot_smoothed_high,color=color_FUV,linestyle='-.')
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_high)],bJU_tot_smoothed_high,color=color_ULTRASAT,linestyle='-.',label=r'$f_{\rm DM}\Gamma = 10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM_high))
#
    #plt.plot(z_gals('SDSS')[:len(bJGF_tot_smoothed_2)],bJGF_tot_smoothed_2,alpha=0.3,color=color_FUV,)
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_2)],bJU_tot_smoothed_2,alpha=0.3,label=r'$m_{\rm DM}c^2 = %g\,{\rm eV}\,$'%m_DM_2 + r'$(10^{%g}{\rm s}^{-1})$'%np.log10(f_DM_decay_DM))
#
    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.ylabel(r'$b_JdJ_{\nu_{\rm obs}}/dz\,[{\rm Jy/sr}]$',fontsize=fontsize_label)
    plt.legend(fontsize=fontsize_legend,loc=2)

    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    plt.xlim(0.1,2.)
    plt.ylim(10,1e4)
    plt.yscale('log')

    plt.tight_layout()
    plt.subplots_adjust(hspace=.37)
    plt.savefig('results/PLOTS/EBL_DM/DM_summary.png',bbox_inches='tight')

    plt.show()

    return 


def wJgz_DM(z,detector,gal_survey,m_DM,f_DM_decay_DM,
            run=False,
            vals_eps1500=False,
            vals_alpha1500=False,
            vals_alpha1100=False,
            val_EW=False,
            val_flyc=False,
            val_alpha900=False,
            val_bias=False,
            filename = '',
            filename_wm = '',
            filename_dJ = '',
            filename_bJ = ''):

    if not run and os.path.exists(filename):

        z_arr, bar_wJg_all = np.genfromtxt(filename)
        bar_wJg = interp1d(z_arr,bar_wJg_all)(z)

    else:

        bJdJdz = lambda zv: bJ_z_DM(zv, detector,m_DM,f_DM_decay_DM,False,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,val_bias,filename_bJ) * dJdz_DM(zv, detector, m_DM,f_DM_decay_DM,False, vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,filename_dJ)#.value

        bar_wm = lambda zv: wz_m(zv,detector,False,filename_wm)

        bg = lambda zv: b_gal(zv,gal_survey)  

        intg = lambda zv: bar_wm(zv) * bg(zv) * bJdJdz(zv)

        bar_wJg = intg(z) 

    return bar_wJg


def run_wJg_DM(m_DM, f_DM_decay_DM, reduced_z = True,):

    dir_DM = './results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    if reduced_z:
        use_z = z_gals_interp
        use_z_ultrasat = z_gals_interp
        reduced_label = '_reduced'
    else:
        use_z = z_gals('DESI')
        reduced_label = ''
        use_z_ultrasat = z_gals('DESI')
        reduced_label = ''

    wn = np.zeros(len(use_z))
    wf = np.zeros(len(use_z))
    wu = np.zeros(len(use_z_ultrasat))


    for i in (range(len(use_z))):
       print('\nDoing z = ' + str(use_z[i]))
       print('NUV')
       filename_wm = 'results/EBL/wmz_GALEX_NUV' + reduced_label + '.dat'
       wn[i] = wJgz(use_z[i],'GALEX_NUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_GALEX_NUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_GALEX_NUV' + reduced_label + '.dat')
       print('FUV')
       filename_wm = 'results/EBL/wmz_GALEX_FUV' + reduced_label + '.dat'
       wf[i] = wJgz(use_z[i],'GALEX_FUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_GALEX_FUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_GALEX_FUV' + reduced_label + '.dat')
#
    for i in (range(len(use_z_ultrasat))):
        print('ULTRASAT')
        filename_wm = 'results/EBL/wmz_ULTRASAT' + reduced_label + '.dat'
        wu[i] = wJgz(use_z_ultrasat[i],'ULTRASAT','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_ULTRASAT' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_ULTRASAT' + reduced_label + '.dat')
        
    np.savetxt(dir_DM + 'wJg_GALEX_NUV,DESI' + reduced_label + '.dat',(use_z, wn))
    np.savetxt(dir_DM + 'wJg_GALEX_FUV,DESI' + reduced_label + '.dat',(use_z, wf))
    filename_ULT = dir_DM + 'wJg_ULTRASAT,DESI' + reduced_label + '.dat'

    np.savetxt(filename_ULT,(use_z_ultrasat, wu))
    
    return 



pars_fid_DM = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2','log_fLyC_1','log_fLyC_2','bias_1500_0','gamma_bv','gamma_bz','alpha_900','log10_fGDM']
#
fid_vals_det_DM = lambda detector, log10_f_DM_decay_DM: [fiducials['eps1500'][0],
      fiducials['eps1500'][1],
      fiducials['alpha1500'][0],
      fiducials['alpha1500'][1],
      fiducials['alpha1100'][0],
      fiducials['alpha1100'][1],
      fiducials['EW'](detector)[0],
      fiducials['EW'](detector)[1],
      fiducials['fescape'][0],
      fiducials['fescape'][1],
      fiducials['bias'][0],
      fiducials['bias'][1],
      fiducials['bias'][2],
        fiducials['alpha900'],
      log10_f_DM_decay_DM]

def ders_DM(parameter,detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z = False):

    reduced_label = '_reduced' if reduced_z else ''

    print('Doing der: ' + str(parameter))

    dir_DM = './results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    filename = dir_DM + 'der/d' + parameter + '_' + detector + ',' + gal_survey + reduced_label + '.dat'

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT':
        filename_wm = 'results/EBL/wmz_' + detector + reduced_label + '.dat'
    else:
        filename_wm = -1

    filename_dJ = dir_DM + 'der/dJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'
    filename_bJ = dir_DM + 'der/bJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'

    filename_dJ_fid = dir_DM + 'dJdz_' + detector  + reduced_label +  '.dat'
    filename_bJ_fid = dir_DM + 'bJ_' + detector  + reduced_label +  '.dat' 

    fid_vals = fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))
    if os.path.exists(filename):
        use_z = z_gals_interp 

        z_arr, der_wz_arr = np.genfromtxt(filename)
        der_wz = interp1d(z_arr,der_wz_arr)(use_z)

    else:

        if reduced_z:
            use_z = z_gals_interp if detector == 'GALEX_NUV' or detector ==    'GALEX_FUV' or detector == 'ULTRASAT' else -1
        else:
            use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector ==    'GALEX_FUV' or detector == 'ULTRASAT' else -1

        use_pars_up = []
        for i in range(len(pars_fid_DM)):
            
            if pars_fid_DM[i] == parameter:
                par_up = fid_vals[i]*(1+delta_der)
                step = fid_vals[i]*delta_der
            else:
                par_up = fid_vals[i]

            use_pars_up.append(par_up)

        with Pool(6) as pool:

            print('Doing parameter up')
            dJ_up_f = partial(dJdz_DM,detector=detector,
                    m_DM=m_DM,\
                    f_DM_decay_DM=10**use_pars_up[pars_fid_DM.index('log10_fGDM')],\
                    run=True,\
                    vals_eps1500=[use_pars_up[pars_fid_DM.index('eps_1500_0')],use_pars_up[pars_fid_DM.index('gamma_1500')]],\
                    vals_alpha1500=[use_pars_up[pars_fid_DM.index('alpha_1500_0')],use_pars_up[pars_fid_DM.index('C_alpha_1500')]],\
                    vals_alpha1100=[use_pars_up[pars_fid_DM.index('alpha_1100_0')],use_pars_up[pars_fid_DM.index('C_alpha_1100')]],\
                    val_EW=[use_pars_up[pars_fid_DM.index('EW_z1')],use_pars_up[pars_fid_DM.index('EW_z2')]],\
                    val_flyc=[use_pars_up[pars_fid_DM.index('log_fLyC_1')],use_pars_up[pars_fid_DM.index('log_fLyC_2')]],\
                    val_alpha900=use_pars_up[pars_fid_DM.index('alpha_900')],\
                    filename=filename_dJ)
            
            dJ_up = pool.map(dJ_up_f, use_z)

            print('Doing bias up')
            bJ_up_f = partial(bJ_z_DM,detector=detector,\
                    m_DM=m_DM,\
                    f_DM_decay_DM=10**use_pars_up[pars_fid_DM.index('log10_fGDM')],\
                    run=True,\
                    vals_eps1500=[use_pars_up[pars_fid_DM.index('eps_1500_0')],use_pars_up[pars_fid_DM.index('gamma_1500')]],\
                    vals_alpha1500=[use_pars_up[pars_fid_DM.index('alpha_1500_0')],use_pars_up[pars_fid_DM.index('C_alpha_1500')]],\
                    vals_alpha1100=[use_pars_up[pars_fid_DM.index('alpha_1100_0')],use_pars_up[pars_fid_DM.index('C_alpha_1100')]],\
                    val_EW=[use_pars_up[pars_fid_DM.index('EW_z1')],use_pars_up[pars_fid_DM.index('EW_z2')]],\
                    val_flyc=[use_pars_up[pars_fid_DM.index('log_fLyC_1')],use_pars_up[pars_fid_DM.index('log_fLyC_2')]],\
                    val_alpha900=use_pars_up[pars_fid_DM.index('alpha_900')],\
                    val_bias=[use_pars_up[pars_fid_DM.index('bias_1500_0')],use_pars_up[pars_fid_DM.index('gamma_bv')],use_pars_up[pars_fid_DM.index('gamma_bz')]],\
                    filename=filename_bJ)
            
            bJ_up = pool.map(bJ_up_f, use_z)

        np.savetxt(filename_dJ,(use_z,dJ_up))
        np.savetxt(filename_bJ,(use_z,bJ_up))


        filename_fid = dir_DM + 'wJg_' + detector + ',' + gal_survey + reduced_label + '.dat'
       
        if os.path.exists(filename_fid):
                run_fid = False 
        else:
                run_fid = True       
       
        der_wz = np.zeros(len(use_z))
        for z in tqdm(range(len(use_z))):
            wz_up = wJgz_DM(use_z[z],detector,gal_survey, \
                m_DM,\
                10**use_pars_up[pars_fid_DM.index('log10_fGDM')], True,
                vals_eps1500=[use_pars_up[pars_fid_DM.index('eps_1500_0')],use_pars_up[pars_fid_DM.index('gamma_1500')]],
                vals_alpha1500=[use_pars_up[pars_fid_DM.index('alpha_1500_0')],use_pars_up[pars_fid_DM.index('C_alpha_1500')]],
                vals_alpha1100=[use_pars_up[pars_fid_DM.index('alpha_1100_0')],use_pars_up[pars_fid_DM.index('C_alpha_1100')]],
                val_EW=[use_pars_up[pars_fid_DM.index('EW_z1')],use_pars_up[pars_fid_DM.index('EW_z2')]],
                val_flyc=[use_pars_up[pars_fid_DM.index('log_fLyC_1')],use_pars_up[pars_fid_DM.index('log_fLyC_2')]],
                val_alpha900=use_pars_up[pars_fid_DM.index('alpha_900')],
                val_bias=[use_pars_up[pars_fid_DM.index('bias_1500_0')],use_pars_up[pars_fid_DM.index('gamma_bv')],use_pars_up[pars_fid_DM.index('gamma_bz')]],\
                #filename = filename, 
                filename_wm = filename_wm,
                filename_dJ = filename_dJ,
                filename_bJ =filename_bJ)

            wz = wJgz_DM(use_z[z],detector,gal_survey,m_DM,f_DM_decay_DM,
                run=run_fid,
                vals_eps1500=False,
                vals_alpha1500=False,
                vals_alpha1100=False,
                val_EW=False,
                val_flyc=False,
                val_alpha900=False,
                val_bias=False,
                filename =filename_fid, 
                filename_wm = filename_wm,
                filename_dJ = filename_dJ_fid,
                filename_bJ =filename_bJ_fid)

            try:
                der_wz[z] = ((wz_up - wz)/step)
            except:
                der_wz[z] = ((wz_up - wz)/step).value

        np.savetxt(filename,(use_z, der_wz))
        if run_fid:
            np.savetxt(filename_fid,(use_z, (use_z,wz)))

    return der_wz


def run_ders_DM(detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z=False):

    for parameter in pars_fid_DM:
        ders_DM(parameter,detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z)

    return

def Jnu_monopole_DM(detector,m_DM,f_DM_decay_DM):

    dir_DM = './results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    gal_survey = 'SPHEREx' if detector == 'ULTRASAT' else 'SDSS'
    use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

    reduced_label = '_reduced' #if reduced_z else ''

    filename = dir_DM + 'dJdz_' + detector + reduced_label + '.dat'
    
    J = np.zeros(len(use_z))
    for i in range(len(use_z)):

        J[i] = dJdz_DM(use_z[i], detector,m_DM,f_DM_decay_DM,run = False,\
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=filename)

    Jnu = np.trapz(J,use_z)

    return Jnu


def sigma_wz_DM(z, detector, gal_survey,m_DM,f_DM_decay_DM, group_vox):

    delta_zc = 1e-2

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV':
        if not group_vox:
            px_size = 5*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 50.*u.arcsec

    elif detector == 'ULTRASAT':
        if not group_vox:
            px_size = 5.4*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 5.4*u.arcsec

    if gal_survey == 'SDSS':
        Omega_survey = (5500*u.deg**2).to(u.steradian)
    elif gal_survey == 'DESI':
        Omega_survey = (14000*u.deg**2).to(u.steradian)
    else:
        print('\nCheck Omega survey\n')
        return 

    Avox = (((px_size)**2).to(u.steradian))
    dzi = delta_zi(gal_survey)

    angle_observed =  (np.pi*(theta_max(cosmo.angular_diameter_distance(z),detector))**2*u.steradian)

    noise = (dzi) / delta_zc * np.sqrt(Avox  * (Jnu_monopole_DM(detector,m_DM,f_DM_decay_DM)**2 + sigma_N(detector, group_vox).value**2) / (dNgdz(z,gal_survey)*delta_zi_interp )  / angle_observed)


    return noise



def fisher_matrix_DM(pars,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox,run = False):
    
    dir_DM = './results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    if group_vox: 
        filename = dir_DM + 'FISHER_' + detector + ',' + gal_survey + '.dat'
    else:
        filename = dir_DM + 'FISHER_' + detector + ',' + gal_survey + '_singleVOX.dat'

    if not run and os.path.exists(filename):
        Fisher_start = np.genfromtxt(filename)

        Fisher = np.zeros((len(pars),len(pars)))
        for i in range(len(pars)):
            if pars[i] not in pars_fid_DM:
                print('Check the parameter set!!!!')
                return 
            for j in range(len(pars)):
                Fisher[i,j] = Fisher_start[pars_fid_DM.index(pars[i])][pars_fid_DM.index(pars[j])]
        return Fisher

    Fisher = np.zeros((len(pars),len(pars)))
    
    use_z = z_gals_interp #z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

    der_pars = {}
    for p in range(len(pars)):
        #temp = ders(pars[p],detector,gal_survey,reduced_z = True)
        #all_zbin_der = interp1d(use_z,temp)
        all_zbin_der = ders_DM(pars[p],detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z = True)

        der_pars[pars[p]] = all_zbin_der 

    sigma2 = np.zeros(len(use_z))
    print('Doing sigma2')
    for zv in tqdm(range(len(use_z))):

        sigma2[zv] = sigma_wz_DM(use_z[zv],detector,gal_survey,m_DM,f_DM_decay_DM,group_vox).value**2

    Fisher = np.zeros((len(pars),len(pars)))
    print('Doing Fisher')
    for p1 in tqdm(range(len(pars))):
        for p2 in range(len(pars)):
            if p1 <= p2:
                der_p1 = np.asarray(der_pars[pars[p1]])                
                der_p2 = np.asarray(der_pars[pars[p2]])        
                        
                temp = sum(der_p1*der_p2/sigma2)
                    
                Fisher[p1,p2] = temp
            else: 
                Fisher[p1,p2] = Fisher[p2,p1]

    if group_vox: 
        np.savetxt(filename,Fisher)
    else:
        np.savetxt(filename,Fisher)

    return Fisher


def Fisher_change_var_DM(pars,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox,run=False):

    F = fisher_matrix_DM(pars_fid_DM,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox, run)

#    dlogeps_dP = 1 
    deps_dP = np.log(10) * fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index('eps_1500_0')]
    db_dP = np.log(10) * fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index('bias_1500_0')]

    Fisher_prime = np.zeros((len(pars),len(pars)))
    
    for i in range(len(pars)):

        if pars[i] == 'log10_epsbias_1500_0':

            Fisher_prime[i,i] = F[pars_fid_DM.index('eps_1500_0')][pars_fid_DM.index('eps_1500_0')] * deps_dP**2 + 2*F[pars_fid_DM.index('eps_1500_0')][pars_fid_DM.index('bias_1500_0')] * deps_dP*db_dP + F[pars_fid_DM.index('bias_1500_0')][pars_fid_DM.index('bias_1500_0')] * db_dP**2

            for j in range(len(pars)):
                
                if i < j:
                    Fisher_prime[i,j] = F[pars_fid_DM.index('eps_1500_0')][pars_fid_DM.index(pars[j])] * deps_dP + F[pars_fid_DM.index('bias_1500_0')][pars_fid_DM.index(pars[j])] * db_dP

        else:
            for j in range(len(pars)):
                if i <= j:
                    Fisher_prime[i,j] = F[pars_fid_DM.index(pars[i])][pars_fid_DM.index(pars[j])]

    for i in range(len(pars)):
        for j in range(len(pars)):
            if i > j:
                Fisher_prime[i,j] = Fisher_prime[j,i]

    return Fisher_prime


pars_all_DM = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2','log_fLyC_1','log_fLyC_2','gamma_bv','gamma_bz','alpha_900','log10_fGDM']


def compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT', pars =pars_original_c18, prior = False, plot_flag = True, group_vox = True):

    print('Doing contour plots')

    if detector == 'GALEX':

       F_N = Fisher_change_var_DM(pars,'GALEX_NUV','SDSS',m_DM,f_DM_decay_DM,group_vox,False)
       F_F = Fisher_change_var_DM(pars,'GALEX_FUV','SDSS',m_DM,f_DM_decay_DM,group_vox,False)
       temp = F_N + F_F

    elif detector == 'GALEX_DESI':

        F_N = Fisher_change_var_DM(pars,'GALEX_NUV','DESI',m_DM,f_DM_decay_DM,group_vox,False)
        F_F = Fisher_change_var_DM(pars,'GALEX_FUV','DESI',m_DM,f_DM_decay_DM,group_vox,False)
        temp = F_N + F_F

    elif detector == 'ULTRASAT_DESI':

        F = Fisher_change_var_DM(pars,'ULTRASAT','DESI',m_DM,f_DM_decay_DM,group_vox,False)
        temp = F

    elif detector == 'GALEX_ULTRASAT_DESIDESI':

        F_N = Fisher_change_var_DM(pars,'GALEX_NUV','DESI',m_DM,f_DM_decay_DM,group_vox,False)
        F_F = Fisher_change_var_DM(pars,'GALEX_FUV','DESI',m_DM,f_DM_decay_DM,group_vox,False)
        F_GALEX = F_N + F_F
        F_ULTRASAT = Fisher_change_var_DM(pars,'ULTRASAT','DESI',m_DM,f_DM_decay_DM,group_vox,False)

        F_both = np.zeros((len(pars)+1,len(pars)+1))
        for i in range(len(pars)):
            if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
            elif pars[i] == 'EW_z1':
                for j in range(len(pars)):
                    F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
            elif pars[i] == 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
        for j in range(len(pars)):
            if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
            elif pars[j] == 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
        F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]

        temp = F_both

    if prior: 
        for j in range(len(pars)):
            # if pars[j] == 'gamma_1500':
            #     temp[j,j] += 1/.3**2
            if detector == 'GALEX':
                if pars[j] == 'C_alpha_1500':
                    temp[j,j] += 1/1.5**2
                if pars[j] == 'C_alpha_1100':
                    temp[j,j] += 1/1.5**2
            if pars[j] == 'gamma_bz':
                temp[j,j] += 1/0.3**2
            if pars[j] == 'gamma_bv':
               temp[j,j] += 1/1.3**2

    inv_F = np.linalg.inv(temp)


    names = []
    fiducials_pars = []    

    if detector == 'GALEX_ULTRASAT' or detector=='GALEX_ULTRASAT_DESI' or detector=='GALEX_ULTRASAT_DESIDESI':
        for i in pars:
            if i == 'log10_epsbias_1500_0':
                fiducials_pars.append(np.log10(fid_vals_det_DM('GALEX_NUV',np.log10(f_DM_decay_DM))[pars_fid_DM.index('eps_1500_0')]*fid_vals_det_DM('GALEX_NUV',np.log10(f_DM_decay_DM))[pars_fid_DM.index('bias_1500_0')]))
            elif i == 'EW_z1' or i == 'EW_z2':
                fiducials_pars.append(fid_vals_det_DM('GALEX_NUV',np.log10(f_DM_decay_DM))[pars_fid_DM.index(i)].value)
            else:
                fiducials_pars.append(fid_vals_det_DM('GALEX_NUV',np.log10(f_DM_decay_DM))[pars_fid_DM.index(i)])

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else  r'$\alpha_{900}$' if i == 'alpha_900' else r'$f_{\rm DM}\Gamma$' if i == 'log10_fGDM' else -1
            names.append(name)

        fiducials_pars.append(fid_vals_det_DM('ULTRASAT',np.log10(f_DM_decay_DM))[pars_fid_DM.index('EW_z2')].value)
        names.append(r'$EW^{z=2}$')

    else:
        for i in pars:
            if i == 'log10_epsbias_1500_0':
                fiducials_pars.append(np.log10(fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index('eps_1500_0')]*fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index('bias_1500_0')]))
            elif i == 'log10_fGDM':
                fiducials_pars.append(np.log10(f_DM_decay_DM))
            else:
                try:
                    fiducials_pars.append(fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index(i)].value)
                except:
                    fiducials_pars.append(fid_vals_det_DM(detector,np.log10(f_DM_decay_DM))[pars_fid_DM.index(i)])

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else r'$\alpha_{900}$' if i == 'alpha_900' else r'$f_{\rm DM}\Gamma$' if i == 'log10_fGDM' else -1

            names.append(name)

    print('DETECTOR: ' + detector)
    for i in range(len(names)-1):
        print('--- ' + names[i] + ': ' + str(fiducials_pars[i]) + ' +- ' + str(round(np.sqrt(inv_F[i][i]),6)))  

    print('--- ' + names[-1] + ': ' + str(fiducials_pars[-1]) + ' +- ' + str(round(np.sqrt(inv_F[-1][-1]),6)))            

    if plot_flag:

        label_det = r'$\rm GALEX\times SDSS$' if detector == 'GALEX' else r'$\rm NUV\times SDSS$' if detector == 'GALEX_NUV' else r'$\rm FUV\times SDSS$' if detector == 'GALEX_FUV' else  r'$\rm ULTRASAT\times SPHEREx$' if detector == 'ULTRASAT' else  r'$\rm ULTRASAT\times DESI$' if detector == 'ULTRASAT_DESI' else r'$\rm Combined$' if detector == 'GALEX_ULTRASAT' else r'$\rm Combined\, with\,DESI$' if detector == 'GALEX_ULTRASAT_DESI' else -1

        color = color_FUV if detector == 'GALEX' else color_ULTRASAT if (detector == 'ULTRASAT' or detector == 'ULTRASAT_DESI') else color_ULTRASAT if (detector == 'GALEX_ULTRASAT' or detector == 'GALEX_ULTRASAT_DESI') else color_FUV if detector == 'GALEX_FUV' else color_NUV if detector == 'GALEX_NUV' else -1

        line_arg = {'ls':'-', 'color':color} 

        plot_dist_list = []
        plot_dist_list.append(GaussianND(fiducials_pars, inv_F, names=names))

        settings = gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = fontsize*.8
        settings.axes_labelsize = fontsize*.8
        settings.legend_fontsize = fontsize*.9
        settings.fontsize = fontsize
        settings.fig_width_inch = 12

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = False
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14
        g.settings.linewidth= linewidth

        g.triangle_plot(plot_dist_list, names,  
        colors=[color],line_args=line_arg,\
        legend_labels=[label_det],legend_loc='upper right',\
            markers={'x2':0})

        for i in range(len(pars)-1):
            ax = g.subplots[i+1, 0]
            ax.yaxis.set_label_coords(-0.7, 0.5)

        plt.tight_layout()

        if detector == 'GALEX_ULTRASAT':
            name_plot = 'ellipse_combined' 
        elif detector == 'GALEX_ULTRASAT_DESI':
            name_plot = 'ellipse_combined_DESI' 
        else:
            name_plot = 'ellipse_' + detector
        filename = 'results/PLOTS/EBL_DM/' + name_plot + '.png'
        plt.savefig(filename,bbox_inches='tight')

        plt.show()

    return

