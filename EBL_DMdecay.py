from EBL_forecast import *
import EBL_plots_paper as p 

from scipy.special import gamma 

rho_crit = (3*H(0)**2/(8*np.pi*cu.G)).to(u.Msun*u.Mpc**-3)


f_gg = 1 # good for axions
F_g = 0 # sets stimulated decay due to CMB - check different values

logaddexp = lambda x,y: x + np.log(1 + np.exp(y - x)) if (x > y) else y + np.log(1 + np.exp(x - y))

sigmoid = lambda lambda_rest: np.exp(-logaddexp(0, -(lambda_rest.value - (1216*u.AA).value) / 100))

alpha_fid = 1.67
m0_fid = 3.82

mass_spectrum_range = np.linspace(7.5,15,11)

#use_dir_0 = 'results/EBL/'
#use_dir = lambda m_DM, f_DM_decay_DM: use_dir_0 + 'DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM)) + '/'

use_dir_0 = 'test_results/EBL_DM/'
use_dir = lambda m_DM, f_DM_decay_DM: use_dir_0 + 'DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM)) + '/'


def ddm_mass_probability_density(m_ddm, alpha, m_0, plot = False):
    term1 = alpha / m_0
    term2 = 1 / gamma(1 + 1/alpha)
    term3 = lambda use_mass: (use_mass / m_0) ** alpha
    term4 = lambda use_mass: np.exp(-((use_mass / m_0) ** alpha))

    if plot:
        m_ddm_values = np.linspace(0,20)
        probability_values = term1 * term2 * term3(m_ddm_values) * term4(m_ddm_values)

        term1_1 = 1.66 / m_0
        term2_1 = 1 / gamma(1 + 1/1.66)
        term3_1 = lambda use_mass: (use_mass / 3.21) ** 1.66
        term4_1 = lambda use_mass: np.exp(-((use_mass / 3.21) ** 1.66))
        probability_values_1 = term1_1 * term2_1 * term3_1(m_ddm_values) * term4_1(m_ddm_values)

        plt.figure(figsize=(10,9))
        fontsize_label = 30
        fontsize_legend = 27
        fontsize_tick = 27

        plt.plot(m_ddm_values, probability_values, 'k', label=r'$\alpha={}, m_0={}$'.format(alpha, m_0))
        plt.plot(m_ddm_values, probability_values_1, 'k--', label=r'$\alpha={}, m_0={}$'.format(1.66, 3.21))
        plt.xlabel(r'$m_{\mathrm{DDM}}c^2$',fontsize=fontsize_label)
        plt.ylabel(r'$p(m_{\rm DDM})$',fontsize=fontsize_label)

        z = [0.1,0.3,0.5]#1.]
        for i in range(len(z)):
            mass_min_galex= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_min('GALEX_FUV'))).to(u.eV)).value
            mass_max_galex= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_max('GALEX_FUV'))).to(u.eV)).value

            mass_min_galex_nuv= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_min('GALEX_NUV'))).to(u.eV)).value
            mass_max_galex_nuv= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_max('GALEX_NUV'))).to(u.eV)).value

            mass_min_ULTRASAT= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_min('ULTRASAT'))).to(u.eV)).value
            mass_max_ULTRASAT= ((4*np.pi*cu.hbar*(1+z[i])*nu_from_lambda(wavelenght_max('ULTRASAT'))).to(u.eV)).value

            plt.fill_betweenx(y=[-0.005,0.03+i/70],x1=mass_min_galex,x2=mass_max_galex,color=color_FUV,alpha=0.6-i/5,label=r'$\rm GALEX\,FUV$')
    # 
            plt.fill_betweenx(y=[-0.005,0.03+i/70],x1=mass_min_galex_nuv,x2=mass_max_galex_nuv,color=color_NUV,alpha=0.6-i/5,label=r'$\rm GALEX\,NUV$')
    # 
            plt.fill_betweenx(y=[-0.005,0.03+i/70],x1=mass_min_ULTRASAT,x2=mass_max_ULTRASAT,color=color_ULTRASAT,alpha=0.6-i/5,label=r'$\rm ULTRASAT$')

            plt.text(mass_max_ULTRASAT,0.031+i/70,r'$z=%g$'%z[i],fontsize=25)
            if i == 0 :
                plt.axvline(mass_min_galex,ymax=0.19, color=color_FUV,linestyle=':')
                plt.axvline(mass_max_galex,ymax=0.19,color=color_FUV,linestyle=':')

                plt.axvline(mass_min_galex_nuv,ymax=0.19,color=color_NUV,linestyle=':')
                plt.axvline(mass_max_galex_nuv,ymax=0.19,color=color_NUV,linestyle=':')

                plt.axvline(mass_min_ULTRASAT,ymax=0.19,color=color_ULTRASAT,linestyle=':')
                plt.axvline(mass_max_ULTRASAT,ymax=0.19,color=color_ULTRASAT,linestyle=':')

                plt.legend(loc = 1,fontsize=fontsize_legend)

        plt.xlim(0.,20.)
        plt.ylim(-0.005,0.185)
        plt.tight_layout()
        plt.yticks(fontsize=fontsize_tick)
        plt.xticks(fontsize=fontsize_tick)

        plt.savefig('results/PLOTS/EBL_DM/massfunction.png',bbox_inches='tight')

        plt.show()

    return term1 * term2 * term3(m_ddm) * term4(m_ddm)


def ddm_nu_probability_density(m_ddm, alpha, m_0, plot = False):
    
    pm  = lambda use_m: ddm_mass_probability_density(use_m, alpha, m_0, plot = False)
    if plot:
        pm_1  = lambda use_m: ddm_mass_probability_density(use_m, 1.66, 3.21, plot = False)
        nu_rest_DM = lambda use_m: ((use_m*u.eV / (4*np.pi*cu.hbar)).to(u.Hz))

    unit = cu.hbar.unit
    p_nu_unn = lambda use_m: ( pm(use_m) * (4*np.pi*cu.hbar)).value
    
    p_nu = lambda use_m: p_nu_unn(use_m) * unit / (quad(p_nu_unn,1e-5,100)[0] * unit * u.eV)

    if plot:

        p_lambda_unn = lambda use_m: -( pm(use_m) * (4*np.pi*cu.hbar) * (cu.c.to(u.Angstrom/u.s)/lambda_from_nu(nu_rest_DM(use_m))**2)).value
        
        p_lambda_unn_1 = lambda use_m: -( pm_1(use_m) * (4*np.pi*cu.hbar) * (cu.c.to(u.Angstrom/u.s)/lambda_from_nu(nu_rest_DM(use_m))**2)).value
        
        unit_lambda = cu.hbar.unit*u.AA/u.s/u.AA**2

        p_lambda = lambda use_m: p_lambda_unn(use_m) * unit_lambda / (quad(p_lambda_unn,1e-5,100)[0]*unit_lambda * u.eV)

        p_lambda_1 = lambda use_m: p_lambda_unn_1(use_m) * unit_lambda / (quad(p_lambda_unn_1,1e-5,100)[0]*unit_lambda * u.eV)

        m_ddm_values = np.linspace(.5,20)

        plt.figure(figsize=(10,9))
        fontsize_label = 27
        fontsize_legend = 23
        fontsize_tick = 23

        plt.plot(lambda_from_nu(nu_rest_DM(m_ddm_values)), p_lambda(m_ddm_values), 'k', label=r'$\alpha, m_0={},{}$'.format(alpha, m_0))
        
        plt.plot(lambda_from_nu(nu_rest_DM(m_ddm_values)), p_lambda_1(m_ddm_values), 'k--', label=r'$\alpha, m_0={},{}$'.format(1.66,3.21))
        plt.xlabel(r'$\lambda_{\mathrm{DM}}\,[{\rm \AA}]$',fontsize=fontsize_label)
        plt.ylabel(r'$p(\lambda_{\rm DDM})$',fontsize=fontsize_label)
        
        plt.axvline(wavelenght_min('GALEX_FUV').value,0,1.5,color=color_FUV)
        plt.axvline(wavelenght_max('GALEX_FUV').value,0,1.5,color=color_FUV)
        plt.fill_betweenx(y=[-0.1,0.2],x1=wavelenght_min('GALEX_FUV').value,x2=wavelenght_max('GALEX_FUV').value,color=color_FUV,alpha=0.2,label=r'$\rm GALEX\,FUV$')

        plt.axvline(wavelenght_min('GALEX_NUV').value,0,1.5,color=color_NUV)
        plt.axvline(wavelenght_max('GALEX_NUV').value,0,1.5,color=color_NUV)
        plt.fill_betweenx(y=[-0.1,0.2],x1=wavelenght_min('GALEX_NUV').value,x2=wavelenght_max('GALEX_NUV').value,color=color_NUV,alpha=0.2,label=r'$\rm GALEX\,NUV$')

        plt.axvline(wavelenght_min('ULTRASAT').value,0,1.5,color=color_ULTRASAT)
        plt.axvline(wavelenght_max('ULTRASAT').value,0,1.5,color=color_ULTRASAT)
        plt.fill_betweenx(y=[-0.1,0.2],x1=wavelenght_min('ULTRASAT').value,x2=wavelenght_max('ULTRASAT').value,color=color_ULTRASAT,alpha=0.2,label=r'$\rm ULTRASAT$')

        plt.legend(loc=1,fontsize=fontsize_legend)

        plt.xlim(1200,35000)
        plt.ylim(-0.01,0.2)
        plt.tight_layout()
        plt.xscale('log')
        plt.yticks(fontsize=fontsize_tick)
        plt.xticks([1500,3000,5000,10000,20000],[r'$1500$',r'$3000$',r'$5000$',r'$10000$',r'$20000$'],fontsize=fontsize_tick)
        plt.subplots_adjust(hspace=.37)
        plt.savefig('results/PLOTS/EBL_DM/pnufunction.png',bbox_inches='tight')

        plt.show()

    return p_nu(m_ddm)


def signal_DM(nu_obs, m_DM, f_DM_decay_DM, run_compare=False, mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    f_DM_decay_DM *= u.s**-1
    
    if mass_function:

        eps_DM = 0.
        for m in range(len(mass_spectrum_range)):        

            nu_rest_DM = (mass_spectrum_range[m] * u.eV / (4*np.pi*cu.hbar)).to(u.Hz)

            z = -1 + nu_rest_DM/ nu_obs  
            if not run_compare:
                if zmin_gal <= z <= zmax_gal:
                    # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
                    eps_DM_m = fLyC(z,False) * sigmoid(lambda_from_nu(nu_rest_DM)) * (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
                else:
                    eps_DM_m = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3
            else:
                if zmin_gal <= z <= zmax_gal:
                    # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
                    eps_DM_m = (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
                else:
                    eps_DM_m = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3

            eps_DM += eps_DM_m * ddm_nu_probability_density(mass_spectrum_range[m], alpha, m_0, plot = False) * u.eV

    else:
        m_DM *= u.eV

        nu_rest_DM = (m_DM / (4*np.pi*cu.hbar)).to(u.Hz)

        z = -1 + nu_rest_DM/ nu_obs  
        if not run_compare:
            if zmin_gal <= z <= zmax_gal:
                # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
                eps_DM = fLyC(z,False) * sigmoid(lambda_from_nu(nu_rest_DM)) * (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
            else:
                eps_DM = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3 
        else:
            if zmin_gal <= z <= zmax_gal:
                # * (1+z)**3 the factor should not be included, compare eq 2 chiang menard with eq A1 in Jose 
                eps_DM = (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions
            else:
                eps_DM = 0. * u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3

    return eps_DM


def signal_DM_zdep(z, m_DM, f_DM_decay_DM, mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    f_DM_decay_DM *= u.s**-1

    if mass_function:

        eps_DM = 0.

        for m in range(len(mass_spectrum_range)):
        
            nu_rest_DM =( mass_spectrum_range[m] / (4*np.pi*cu.hbar) * u.eV).to(u.Hz) #/ (1+z)

            eps_DM_m = fLyC(z,False) * sigmoid(lambda_from_nu(nu_rest_DM)) * (f_gg * f_DM_decay_DM * camb_pars.omegac * rho_crit  * cu.c**2 * (1+F_g)/nu_rest_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) # need to be weighted by absorptions

            eps_DM += eps_DM_m * ddm_nu_probability_density(mass_spectrum_range[m], alpha, m_0, plot = False) * u.eV

    else:
        m_DM *= u.eV

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
        s_tot_spectrum = np.zeros(len(rest_wave)+1)
        wave = np.zeros((len(rest_wave)+1))
        unsort_wave = np.zeros((len(rest_wave)+1))
        for i in range(len(rest_wave)):

            unsort_wave[i] = rest_wave[i] * (1.+z[zv])

        unsort_wave[-1] = lambda_obs_DM

        wave = np.sort(unsort_wave)

        m_range = np.linspace(10,25,100)
        s_spectrum = np.zeros((len(m_range),len(rest_wave)+1))
        for m in range(len(m_range)):
            lambda_obs_DM_m = lambda_from_nu(m_range[m]*u.eV / (4*np.pi*cu.hbar) / (1+z[zv])).value

            for i in range(len(wave)):

                detector = 'ULTRASAT' # same for galex_fuv in EW

                if wave[i] == lambda_obs_DM_m:
                    use_signal_DM_spectrum = signal_DM_zdep(z[zv], m_range[m], f_DM_decay,True,alpha_fid,m0_fid).value
                else:
                    use_signal_DM_spectrum = 0.

                s_spectrum[m,i] = use_signal_DM_spectrum

        for i in range(len(wave)):

            s_tot_spectrum[i] = signal(wave[i]/(1+z[zv])*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value + np.trapz(s_spectrum[:,i],m_range)
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
            #for m in range(len(m_range)):
            plt.plot(wave, s_tot_spectrum, color= 'k',linestyle='--')
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

    plt.text(2455*(1+1), 6e25, r'$m_{\rm DDM}c^2$', rotation=-1.5, fontsize=fontsize_legend, ha='center', va='center')
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
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = '',with_foreground=False,run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

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

        delta_fg = 0.

        if with_foreground:
            if detector == 'GALEX_NUV' or detector == 'ULTRASAT':
                delta_fg = 2.2
            elif detector == 'GALEX_FUV':
                delta_fg = 1.8

        obs_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(obs_wave))
        for i in range(len(obs_wave)):
            nu_obs_arr[i] = nu_from_lambda(obs_wave[i]*u.AA).value

        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM,run_compare,mass_function = mass_function, alpha=alpha,m_0=m_0) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)*(1+delta_fg)

        unit = (Response(lambda_from_nu(1.*u.Hz), detector) *sum_signal(1*u.Hz).unit * np.exp(-tau_Lya(lambda_from_nu(1*u.Hz),0.))/(1*u.Hz)).unit

        if not run_compare:
            intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
        else:
            intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value /nu_obs 

        intg_arr = np.zeros(len(obs_wave))
        for i in range(len(obs_wave)):
            intg_arr[i] = intg(nu_obs_arr[i])

        dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * np.trapz(intg_arr,nu_obs_arr)*(unit*u.Hz/u.steradian) 

    else:
        zval, dJdzval = np.genfromtxt(filename)
        dJdz = interp1d(zval,dJdzval)(z) * u.Jy/u.steradian  

    return (dJdz.to(u.Jy/u.steradian)).value


def bJ_z_DM(z, detector,  m_DM, f_DM_decay_DM,run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = '',run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

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
        
        sum_signal = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM, run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0) + signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900)

        sum_signal_bias = lambda nu_obs: signal_DM(nu_obs, m_DM, f_DM_decay_DM,run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0) + bias(nu_obs*(1+z),z,val_bias)*signal(lambda_from_nu(nu_obs*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,) 

        if not run_compare:
            intg_num = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal_bias(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
    #
            intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
        else:
            intg_num = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal_bias(nu_obs*u.Hz).value /nu_obs 
    #
            intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * sum_signal(nu_obs*u.Hz).value /nu_obs 

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


def compute_vals_DM(m_DM, f_DM_decay_DM, reduced_z = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    reduced_label = '_reduced' if reduced_z else ''

    dir_DM = use_dir(m_DM,f_DM_decay_DM)
    #'./results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    if not os.path.exists(dir_DM):
        os.makedirs(dir_DM)

    if run_compare:
        dir_DM += 'new_compare_JLB/'
    if not os.path.exists(dir_DM):
        os.makedirs(dir_DM)

    if mass_function:
        dir_DM += 'spectrum_'

    with Pool(6) as pool:
        
        # print('Doing NUV')
        # dJdz_NUV = partial(dJdz_DM,detector='GALEX_NUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900, filename=dir_DM+'dJdz_GALEX_NUV' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)
        # if reduced_z:
        #     dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
        #     np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_nuv)))
        # else:
        #     dJ_nuv = pool.map(dJdz_NUV, z_gals('SDSS'))
        #     np.savetxt(dir_DM+ 'dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_nuv)))

        # print('Doing FUV')
        # dJdz_FUV = partial(dJdz_DM,detector='GALEX_FUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

        # if reduced_z:
        #     dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
        #     np.savetxt(dir_DM+'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_fuv)))
        # else:
        #     dJ_fuv = pool.map(dJdz_FUV, z_gals('SDSS'))
        #     np.savetxt(dir_DM+ 'dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_fuv)))
    
        # print('Doing b NUV')
        # bJ_nuv_f =  partial(bJ_z_DM,detector='GALEX_NUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

        # if reduced_z:
        #     bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
        #     np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_nuv)))
        # else:
        #     bJ_nuv = pool.map(bJ_nuv_f, z_gals('SDSS'))
        #     np.savetxt(dir_DM+ 'bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_nuv)))

        # print('Doing b FUV')
        # bJ_fuv_f =  partial(bJ_z_DM,detector='GALEX_FUV',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

        # if reduced_z:
        #     bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)
        #     np.savetxt(dir_DM+  'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_fuv)))
        # else:
        #     bJ_fuv = pool.map(bJ_fuv_f, z_gals('SDSS'))
        #     np.savetxt(dir_DM+ 'bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_fuv)))
     
        print('Doing ULTRASAT')
        dJdz_f = partial(dJdz_DM,detector='ULTRASAT',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename=dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)
#
        if reduced_z:
            dJ_U = pool.map(dJdz_f, z_gals_interp)
            np.savetxt(dir_DM+'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_U)))
        else:
           dJ_U = pool.map(dJdz_f, z_gals)
           np.savetxt(dir_DM+ 'dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(dJ_U)))

        print('Doing b ULTRASAT')
        bJ_f =  partial(bJ_z_DM,detector='ULTRASAT',m_DM=m_DM, f_DM_decay_DM=f_DM_decay_DM,run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename=dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat',run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)
        if reduced_z:
            bJ_U = pool.map(bJ_f, z_gals_interp)        
            np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_U)))
        else:
           bJ_U = pool.map(bJ_f, z_gals)        
           np.savetxt(dir_DM+ 'bJ_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(bJ_U)))
    return




def plot_collective_DM():

    m_DM_1 = 10.
    m_DM_2 = 20.
    f_DM_decay_DM = 1e-23
    f_DM_decay_DM_low = 1e-25
    f_DM_decay_DM_high = 1e-21

    plt.figure(figsize=(9,8))
    
    dir_DM_1 = use_dir(m_DM_1,f_DM_decay_DM)
    dir_DM_2 = use_dir(m_DM_2,f_DM_decay_DM)

    dir_base = use_results_dir + 'EBL/'
    dir_DM_low = use_results_dir + 'EBL_DM/DM_' + str(m_DM_1) + ',' + str(np.log10(f_DM_decay_DM_low))+'/'
    dir_DM_low_2 = use_results_dir + 'EBL_DM/DM_' + str(m_DM_2) + ',' + str(np.log10(f_DM_decay_DM_low))+'/'
    dir_DM_high = use_results_dir + 'EBL_DM/DM_' + str(m_DM_1) + ',' + str(np.log10(f_DM_decay_DM_high))+'/'

###########################
########################
    dJ_U_tot =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_1 + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U_tot = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_1 + 'bJ_ULTRASAT_reduced.dat')

    dJ_U =  dJdz(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_base + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U = bJ_z(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_base + 'bJ_ULTRASAT_reduced.dat')

    dJ_U_tot_m2 =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_2 + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U_tot_m2 = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_2 + 'bJ_ULTRASAT_reduced.dat')

    dJ_U_tot_low =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_low + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U_tot_low = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_low + 'bJ_ULTRASAT_reduced.dat')

    dJ_U_tot_low_2 =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_low_2 + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U_tot_low_2 = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_2,f_DM_decay_DM=f_DM_decay_DM_low,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_low_2 + 'bJ_ULTRASAT_reduced.dat')


    dJ_U_tot_high =  dJdz_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=dir_DM_high + 'dJdz_ULTRASAT_reduced.dat')

    bJ_U_tot_high = bJ_z_DM(z_gals('DESI'),detector='ULTRASAT',m_DM=m_DM_1,f_DM_decay_DM=f_DM_decay_DM_high,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename=dir_DM_high + 'bJ_ULTRASAT_reduced.dat')


    window_size = 70 # Adjust this value based on your preference
    bJU_smoothed = moving_average(dJ_U*bJ_U, window_size)
# 
    bJU_tot_smoothed = moving_average(dJ_U_tot*bJ_U_tot, window_size)
# 
    bJU_tot_smoothed_low = moving_average(dJ_U_tot_low*bJ_U_tot_low, window_size)
    bJU_tot_smoothed_low_2 = moving_average(dJ_U_tot_low_2*bJ_U_tot_low_2, window_size)

    bJU_tot_smoothed_2 = moving_average(dJ_U_tot_m2*bJ_U_tot_m2, window_size)

    bJU_tot_smoothed_high = moving_average(dJ_U_tot_high*bJ_U_tot_high, window_size)

    fontsize_label = 23
    
    plt.plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_smoothed,color=color_ULTRASAT,label=r'$\rm UV-EBL$')
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed)],bJU_tot_smoothed,color='k',linestyle='--',label=r'$\rm DDM.\, 10\,{\rm eV},\,10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM))
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_2)],bJU_tot_smoothed_2,linestyle='--',color='k',alpha=0.3,label=r'$\rm DDM,\, 20\,{\rm eV},\,10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM))
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_low)],bJU_tot_smoothed_low,color='k',linestyle=':',label=r'$\rm DDM,\, 10\,{\rm eV},\,10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM_low))
    #plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_low_2)],bJU_tot_smoothed_low_2,color='k',linestyle=':',alpha=0.3,label=r'$\rm DDM,\, 20\,{\rm eV},\,10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM_low))
    plt.plot(z_gals('DESI')[:len(bJU_tot_smoothed_high)],bJU_tot_smoothed_high,color='k',linestyle='-.',label=r'$\rm DDM,\, 10\,{\rm eV},\,10^{%g}{\rm s}^{-1}$'%np.log10(f_DM_decay_DM_high))
#
    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.ylabel(r'$b_J(z)dJ_{\nu_{\rm obs}}/dz$',fontsize=fontsize_label)

    plt.xlim(0.1,1.9)
    plt.ylim(5,183)
    #plt.yscale('log')

    plt.legend(loc=3,fontsize=fontsize_label*0.9)
    plt.tight_layout()

    plt.savefig(use_results_dir + 'PLOTS/EBL_DM/DM_summary.png',bbox_inches='tight')

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
            filename_bJ = '',run_compare=False,
            mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    if not run and os.path.exists(filename):

        z_arr, bar_wJg_all = np.genfromtxt(filename)
        bar_wJg = interp1d(z_arr,bar_wJg_all)(z)

    else:

        bJdJdz = lambda zv: bJ_z_DM(zv, detector,m_DM,f_DM_decay_DM,False,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,val_bias,filename_bJ,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0) * dJdz_DM(zv, detector, m_DM,f_DM_decay_DM,False, vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,filename_dJ,False,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)#.value

        bar_wm = lambda zv: wz_m(zv,detector,False,filename_wm)

        bg = lambda zv: b_gal(zv,gal_survey)  

        intg = lambda zv: bar_wm(zv) * bg(zv) * bJdJdz(zv)

        bar_wJg = intg(z) 

    return bar_wJg


def run_wJg_DM(m_DM, f_DM_decay_DM, reduced_z = True, run_compare=False, mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    dir_DM = use_dir(m_DM,f_DM_decay_DM)
    #'./results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    if run_compare:
        dir_DM += 'new_compare_JLB/'

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
       try:
           filename_wm = use_dir_0(m_DM,f_DM_decay_DM) + 'wmz_GALEX_NUV' + reduced_label + '.dat'
       except:
           filename_wm = 'results/EBL/wmz_GALEX_NUV' + reduced_label + '.dat'
       if mass_function:
           wn[i] = wJgz(use_z[i],'GALEX_NUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'spectrum_dJdz_GALEX_NUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'spectrum_bJ_GALEX_NUV' + reduced_label + '.dat')
       else:
           wn[i] = wJgz(use_z[i],'GALEX_NUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_GALEX_NUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_GALEX_NUV' + reduced_label + '.dat')
       print('FUV')
       try:
           filename_wm = use_dir_0(m_DM,f_DM_decay_DM) + 'wmz_GALEX_FUV' + reduced_label + '.dat'
       except:
           filename_wm = 'results/EBL/wmz_GALEX_FUV' + reduced_label + '.dat'
       if mass_function:
           wf[i] = wJgz(use_z[i],'GALEX_FUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'spectrum_dJdz_GALEX_FUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'spectrum_bJ_GALEX_FUV' + reduced_label + '.dat')
       else:
           wf[i] = wJgz(use_z[i],'GALEX_FUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_GALEX_FUV' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_GALEX_FUV' + reduced_label + '.dat')
    for i in (range(len(use_z_ultrasat))):
        print('ULTRASAT')
        try:
           filename_wm = use_dir_0(m_DM,f_DM_decay_DM) + 'wmz_ULTRASAT' + reduced_label + '.dat'
        except:
            filename_wm = 'results/EBL/wmz_ULTRASAT' + reduced_label + '.dat'
        if mass_function:
           wu[i] = wJgz(use_z_ultrasat[i],'ULTRASAT','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'spectrum_dJdz_ULTRASAT' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'spectrum_bJ_ULTRASAT' + reduced_label + '.dat')    
        else:
           wu[i] = wJgz(use_z_ultrasat[i],'ULTRASAT','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = dir_DM + 'dJdz_ULTRASAT' + reduced_label + '.dat',
           filename_bJ = dir_DM + 'bJ_ULTRASAT' + reduced_label + '.dat')

    if mass_function:
        np.savetxt(dir_DM + 'spectrum_wJg_GALEX_NUV,DESI' + reduced_label + '.dat',(use_z, wn))
        np.savetxt(dir_DM + 'spectrum_wJg_GALEX_FUV,DESI' + reduced_label + '.dat',(use_z, wf))
        filename_ULT = dir_DM + 'spectrum_wJg_ULTRASAT,DESI' + reduced_label + '.dat'

        np.savetxt(filename_ULT,(use_z_ultrasat, wu))

    else:
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

def ders_DM(parameter,detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z = False, run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    reduced_label = '_reduced' if reduced_z else ''

    print('Doing der: ' + str(parameter))

    dir_DM = use_dir(m_DM,f_DM_decay_DM)
    #'./results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    if run_compare:
        dir_DM += 'new_compare_JLB/'

    if not os.path.exists(dir_DM + 'der/'):
        os.makedirs(dir_DM+'der/')

    if mass_function:
        filename = dir_DM + 'der/spectrum_d' + parameter + '_' + detector + ',' + gal_survey + reduced_label + '.dat'
    else:
        filename = dir_DM + 'der/d' + parameter + '_' + detector + ',' + gal_survey + reduced_label + '.dat'

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT':
        try:
            filename_wm = use_dir_0 + 'wmz_' + detector + reduced_label + '.dat'
        except:
            filename_wm = 'results/EBL/wmz_' + detector + reduced_label + '.dat'
    else:
        filename_wm = -1

    if mass_function:
        filename_dJ = dir_DM + 'der/spectrum_dJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'
        filename_bJ = dir_DM + 'der/spectrum_bJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'

        filename_dJ_fid = dir_DM + 'spectrum_dJdz_' + detector  + reduced_label +  '.dat'
        filename_bJ_fid = dir_DM + 'spectrum_bJ_' + detector  + reduced_label +  '.dat' 
    else:
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
                    filename=filename_dJ,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)
            
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
                    filename=filename_bJ,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)
            
            bJ_up = pool.map(bJ_up_f, use_z)

        np.savetxt(filename_dJ,(use_z,dJ_up))
        np.savetxt(filename_bJ,(use_z,bJ_up))


        if mass_function:
            filename_fid = dir_DM + 'spectrum_wJg_' + detector + ',' + gal_survey + reduced_label + '.dat'
        else:
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
                filename_bJ =filename_bJ,run_compare=run_compare)

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
                filename_bJ =filename_bJ_fid,run_compare=run_compare)

            try:
                der_wz[z] = ((wz_up - wz)/step)
            except:
                der_wz[z] = ((wz_up - wz)/step).value

        np.savetxt(filename,(use_z, der_wz))
        if run_fid:
            np.savetxt(filename_fid,(use_z, (use_z,wz)))

    return der_wz


def run_ders_DM(detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z=False,run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    for parameter in pars_fid_DM:
        ders_DM(parameter,detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

    return


def Jnu_monopole_DM(detector,m_DM,f_DM_decay_DM, with_foreground = False, run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    dir_DM = use_dir(m_DM,f_DM_decay_DM)
    #'./results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'

    if run_compare:
        dir_DM += 'new_compare_JLB/'

    gal_survey = 'SPHEREx' if detector == 'ULTRASAT' else 'SDSS'
    use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

    reduced_label = '_reduced' #if reduced_z else ''

    if mass_function:
        filename = dir_DM + 'spectrum_dJdz_' + detector + reduced_label + '.dat'
    else:
        filename = dir_DM + 'dJdz_' + detector + reduced_label + '.dat'
    
    J = np.zeros(len(use_z))
    for i in range(len(use_z)):

        J[i] = dJdz_DM(use_z[i], detector,m_DM,f_DM_decay_DM,run = False,\
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=filename,with_foreground=with_foreground,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

    Jnu = np.trapz(J,use_z)

    return Jnu


def plot_monopole():

    nuv_c18 = 259
    fuv_c18 = 90
    
    err_up_nuv_c18 = 62
    err_low_nuv_c18 = 33

    err_up_fuv_c18 = 28
    err_low_fuv_c18 = 16

    nuv_ebl = Jnu_monopole('GALEX_NUV',reduced_z = True)
    fuv_ebl = Jnu_monopole('GALEX_FUV',reduced_z = True)
    ultrasat_ebl = Jnu_monopole('ULTRASAT',reduced_z = True)

    f_decay = [1e-27,1e-26,5e-26,5e-25,1e-24,5e-24,5e-23,1e-22,5e-22,1e-21]
    nuv_10 = np.zeros(len(f_decay))
    fuv_10 = np.zeros(len(f_decay))
    ultrasat_10 = np.zeros(len(f_decay))

    nuv_20 = np.zeros(len(f_decay))
    fuv_20 = np.zeros(len(f_decay))
    ultrasat_20 = np.zeros(len(f_decay))

    for i in range(len(f_decay)):
        nuv_10[i] = Jnu_monopole_DM('GALEX_NUV',10,f_decay[i],False)
        fuv_10[i] = Jnu_monopole_DM('GALEX_FUV',10,f_decay[i],False)
        #ultrasat_10[i] = Jnu_monopole_DM('ULTRASAT',10,f_decay[i],False)

        nuv_20[i] = Jnu_monopole_DM('GALEX_NUV',20,f_decay[i],False)
        fuv_20[i] = Jnu_monopole_DM('GALEX_FUV',20,f_decay[i],False)
        #ultrasat_20[i] = Jnu_monopole_DM('ULTRASAT',20,f_decay[i],False)


    plt.figure(figsize=(8,6))
#    plt.semilogx(f_decay,nuv_c18*np.ones(len(f_decay)),color_NUV,label=r'$\rm GALEX\,NUV\,C18$')
#    plt.semilogx(f_decay,fuv_c18*np.ones(len(f_decay)),color_FUV,label=r'$\rm GALEX\,FUV\,C18$')

    plt.fill_between(f_decay, np.ones(len(f_decay))*(nuv_c18+err_up_nuv_c18),np.ones(len(f_decay))*(nuv_c18-err_low_nuv_c18), color=color_NUV, alpha = 0.4,label=r'$\rm GALEX\,NUV$')

    plt.fill_between(f_decay, np.ones(len(f_decay))*(fuv_c18+err_up_fuv_c18),np.ones(len(f_decay))*(fuv_c18-err_low_fuv_c18), color=color_FUV, alpha = 0.4,label=r'$\rm GALEX\,FUV$')

    plt.semilogx(f_decay,nuv_ebl*np.ones(len(f_decay)),color_NUV,linestyle='-')
    plt.semilogx(f_decay,fuv_ebl*np.ones(len(f_decay)),color_FUV,linestyle='-')
    #plt.loglog(f_decay,ultrasat_ebl*np.ones(len(f_decay)),color_ULTRASAT,linestyle=':',label=r'$\rm ULTRASAT\,LK23$')

    #plt.loglog(f_decay,ultrasat_10,color_ULTRASAT,linestyle='--',label=r'$ULTRASAT,\,m_{\rm DM}c^2=10\,{\rm eV}$')
    plt.semilogx(f_decay,nuv_10,'k',linestyle='--',label=r'$m_{\rm DDM}c^2=10\,{\rm eV}$')
    plt.semilogx(f_decay,fuv_10,'k',linestyle='--')
    plt.semilogx(f_decay,nuv_20,'k',linestyle='-.',label=r'$m_{\rm DDM}c^2=20\,{\rm eV}$')
    plt.semilogx(f_decay,fuv_20,'k',linestyle='-.')

    #plt.loglog(f_decay,ultrasat_20,color_ULTRASAT,linestyle='-.')

    plt.vlines(3e-24,4e1,4e2,'k',linewidth=1)

    plt.xlabel(r'$f_{\rm DDM}\Gamma\,[{\rm s}^{-1}]$',fontsize=22)    
    plt.ylabel(r'$J_{\nu_{\rm obs}}$',fontsize=22)    

    plt.xlim(1e-27,1e-22)
    plt.ylim(4e1,4e2)
    plt.tick_params(axis='x', labelsize=20) 
    plt.tick_params(axis='y', labelsize=20) 

    plt.legend(loc=2,ncol=1,fontsize=20)

    #axs[2].tick_params(axis='x', labelsize=fontsize_tick) 
    #axs[2].tick_params(axis='x', labelsize=fontsize_tick) 
        
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL_DM/monopole.png',bbox_inches='tight')

    plt.show()

    return 
    


def sigma_wz_DM(z, detector, gal_survey,m_DM,f_DM_decay_DM, group_vox,run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):

    delta_zc = max(1e-3,H(z)*scale_physical_min/(cu.c).to(u.km/u.s) * (1+z))

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
    #dzi = delta_zi(gal_survey)
    dzi = max(delta_zi(gal_survey),delta_zc)

    use_delta_zi = delta_zi_interp if gal_survey == 'SDSS' else delta_zi_interp if gal_survey == 'DESI' else -1

    angle_observed =  (np.pi*(theta_max(cosmo.angular_diameter_distance(z),detector))**2*u.steradian)

    ratio_factors = dzi / delta_zc

    noise = ratio_factors * np.sqrt(Avox  * (Jnu_monopole_DM(detector,m_DM,f_DM_decay_DM,True,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)**2 + sigma_N(detector, group_vox).value**2) / (dNgdz(z,gal_survey)*use_delta_zi )  / angle_observed)


    return noise



def run_all_analysis(run_compare=False, m_DM_range = False, f_DM_decay_DM_range = False,):

    if not m_DM_range:
        m_DM_range = np.concatenate((np.linspace(10,30,6),np.array((27.,30.))))
    if not f_DM_decay_DM_range:
        f_DM_decay_DM_range = [1e-23,1e-25,1e-26,1e-27,5e-27,1e-28,5e-28,1e-29]

    for m_DM in m_DM_range:
        print('Doing mass: ' + str(m_DM))
        for f_DM_decay_DM in tqdm(f_DM_decay_DM_range):
            compute_vals_DM(m_DM, f_DM_decay_DM, reduced_z = True, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,run_compare=run_compare,mass_function =False,)
            run_wJg_DM(m_DM, f_DM_decay_DM, reduced_z = True,run_compare=run_compare,mass_function =False, )
            run_ders_DM('GALEX_FUV','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =False, )                
            run_ders_DM('GALEX_NUV','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =False, )                
            run_ders_DM('ULTRASAT','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =False,)

    return 


def run_all_analysis_spectrum(run_compare=False, f_DM_decay_DM_range = False,):

    if not f_DM_decay_DM_range:
        f_DM_decay_DM_range = [1e-27,5e-28,1e-28] #1e-27,5e-28,1e-28 #1e-23,1e-25,1e-27,1e-28,1e-29]

    m_DM = False
    #print('Doing mass: ' + str(m_DM))
    for f_DM_decay_DM in tqdm(f_DM_decay_DM_range):
        compute_vals_DM(m_DM, f_DM_decay_DM, reduced_z = True, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,run_compare=run_compare,mass_function =True,alpha=alpha_fid,m_0=m0_fid)
        run_wJg_DM(m_DM, f_DM_decay_DM, reduced_z = True,run_compare=run_compare,mass_function =True,alpha=alpha_fid,m_0=m0_fid )
        run_ders_DM('GALEX_FUV','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =True,alpha=alpha_fid,m_0=m0_fid )                
        run_ders_DM('GALEX_NUV','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =True,alpha=alpha_fid,m_0=m0_fid )                
        run_ders_DM('ULTRASAT','DESI',m_DM,f_DM_decay_DM,reduced_z=True,run_compare=run_compare,mass_function =True,alpha=alpha_fid,m_0=m0_fid)

    return 


def fisher_matrix_DM(pars,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox,run = False,run_compare=False,mass_function =False, alpha=alpha_fid,m_0=m0_fid):
    
    dir_DM = use_dir(m_DM,f_DM_decay_DM)
    print(dir_DM)
    #'./results/EBL/DM_' + str(m_DM) + ',' + str(np.log10(f_DM_decay_DM))+'/'
    if run_compare:
        dir_DM += 'new_compare_JLB/'

    if mass_function:
        dir_DM += 'spectrum_'

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
        all_zbin_der = ders_DM(pars[p],detector,gal_survey,m_DM,f_DM_decay_DM,reduced_z = True,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0)

        der_pars[pars[p]] = all_zbin_der 

    sigma2 = np.zeros(len(use_z))
    print('Doing sigma2')
    for zv in tqdm(range(len(use_z))):

        sigma2[zv] = sigma_wz_DM(use_z[zv],detector,gal_survey,m_DM,f_DM_decay_DM,group_vox,run_compare=run_compare,mass_function =mass_function, alpha=alpha,m_0=m_0).value**2

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


def Fisher_change_var_DM(pars,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox,run=False,run_compare=False,mass_function=False,alpha_DDM=False,m_0_DDM=False):

    if (m_DM and not mass_function) or (not m_DM and mass_function) :
        F = fisher_matrix_DM(pars_fid_DM,detector,gal_survey,m_DM,f_DM_decay_DM,group_vox, run, run_compare=run_compare,mass_function=mass_function,alpha=alpha_DDM,m_0=m_0_DDM)
    else:
        F = fisher_matrix_DM_massspectrum(pars,detector,gal_survey,f_DM_decay_DM,alpha_DDM,m_0_DDM,group_vox,run = False,run_compare=False)
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


def compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT', pars =pars_original_c18, prior = False, group_vox = True, run_compare=False,mass_function=False,alpha_DDM=False,m_0_DDM=False):

    if detector == 'GALEX':

       F_N = Fisher_change_var_DM(pars,'GALEX_NUV','SDSS',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
       F_F = Fisher_change_var_DM(pars,'GALEX_FUV','SDSS',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
       temp = F_N + F_F

    elif detector == 'GALEX_DESI':

        F_N = Fisher_change_var_DM(pars,'GALEX_NUV','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
        F_F = Fisher_change_var_DM(pars,'GALEX_FUV','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
        temp = F_N + F_F

    elif detector == 'ULTRASAT_DESI':

        F = Fisher_change_var_DM(pars,'ULTRASAT','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
        temp = F

    elif detector == 'GALEX_ULTRASAT_DESIDESI':

        F_N = Fisher_change_var_DM(pars,'GALEX_NUV','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
        F_F = Fisher_change_var_DM(pars,'GALEX_FUV','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
        F_GALEX = F_N + F_F
        F_ULTRASAT = Fisher_change_var_DM(pars,'ULTRASAT','DESI',m_DM,f_DM_decay_DM,group_vox,False,run_compare=run_compare,mass_function=mass_function,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)

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
            if pars[j] == 'gamma_1500':
                temp[j,j] += 1/.3**2
            if detector == 'GALEX':
                if pars[j] == 'C_alpha_1500':
                    temp[j,j] += 1/1.5**2
                if pars[j] == 'C_alpha_1100':
                    temp[j,j] += 1/1.5**2
            if pars[j] == 'bias_1500_0':
                temp[j,j] += 1/0.05**2
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

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else  r'$\alpha_{900}$' if i == 'alpha_900' else r'$f_{\rm DM}\Gamma$' if i == 'log10_fGDM' else r'$b_{1500}^{z=0}$' if i == 'bias_1500_0' else -1
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

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else r'$\alpha_{900}$' if i == 'alpha_900' else r'$f_{\rm DM}\Gamma$' if i == 'log10_fGDM' else r'$b_{1500}^{z=0}$' if i == 'bias_1500_0' else -1

            names.append(name)

    print('DETECTOR: ' + detector)
    val_DM = 0.
    fid_DM = 0.
    for i in range(len(names)-1):
        print('--- ' + names[i] + ': ' + str(fiducials_pars[i]) + ' +- ' + str(round(np.sqrt(inv_F[i][i]),6)))  
        if names[i] == r'$f_{\rm DM}\Gamma$':
            val_DM = np.sqrt(inv_F[i][i])
            fid_DM = fiducials_pars[i]

    return val_DM, fid_DM


#################################################
#################################################
#################################################

from shapely.geometry import LineString

def results_DDM(get_theta_min = False):

    m_DM_range = [10., 13.,16., 18.,20.,25.] #np.linspace(10,25,6)

    if get_theta_min:
        f_DM_decay_DM_range = [1e-27,1e-28]
    else:
        f_DM_decay_DM_range = [1e-23,1e-25,1e-27,1e-28,1e-29]

        colors= ['#9b2226','#ca6702','#eeb831','#7f9c5e','#0a9396','#005c63']
        plt.figure(figsize=(9,8))
        fontsize_label = 27
        fontsize_legend = 23
        fontsize_tick = 23

    intersec = []
    for m_DM in m_DM_range:
        sigma = []
        fid = []
        sigma_ne = []
        fid_ne = []
        for f_DM_decay_DM in f_DM_decay_DM_range:
            temp = compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=get_theta_min,alpha_DDM=False,m_0_DDM=False)
            sigma.append(abs(temp[0]/temp[1]))
            fid.append(temp[1])
            #temp_ne = compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=True,alpha_DDM=False,m_0_DDM=False)
            #sigma_ne.append(abs(temp_ne[0]/temp_ne[1]))
            #fid_ne.append(temp_ne[1])

        if get_theta_min:
            first_line = LineString(np.column_stack((fid, np.log10(sigma))))
            second_line = LineString(np.column_stack((fid, np.ones(len(fid))*np.log10(.5))))
            intersection = first_line.intersection(second_line)
            intersec.append(*intersection.xy[0])

        else:
            plt.plot(fid,np.log10(sigma),'-',label=r'$%g\,{\rm eV}$'%m_DM, zorder=2,color=colors[list(m_DM_range).index(m_DM)])
            #plt.plot(fid_ne,np.log10(sigma_ne),'o--',label=r'$%g\,{\rm eV}$'%m_DM, zorder=2,color=colors[list(m_DM_range).index(m_DM)])

    if get_theta_min:
        return m_DM_range, intersec
    
    plt.axhline(np.log10(.5),-30,color='k',linewidth=1.5)
    plt.fill_between(fid,np.log10(.5),2.9,color='k',alpha=0.25)

    plt.xlim(-29,-23)
    plt.ylim(-4.9,2.9)
    #plt.yscale('log')
    plt.legend(fontsize=fontsize_legend,ncol=2,loc=3)
    plt.xticks(fontsize=fontsize_tick)
    plt.yticks(fontsize=fontsize_tick)
    plt.xlabel(r'$\log_{10}\Theta_{\rm DDM}$',fontsize=fontsize_label)
    plt.ylabel(r'$\log_{10}|\sigma_{\log_{10}\Theta_{\rm DDM}}^{\rm rel}|$',fontsize=fontsize_label)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL_DM/sigma_theta_DDM.png',bbox_inches='tight')

    plt.show()
    return


use_masses_large = [9.5,10., 11., 12.,13.,14.,15., 16., 17.,18.,19.]
def compare_results_DDM(m_DM_range = use_masses_large, plot_flag=True,alpha_DDM=False,m_0_DDM=False, with_run = False, including_ours = False):

    theta_min_arr = []
    if including_ours:
        theta_min_arr_ours = []

    if m_DM_range is not False:
        for m_DM in m_DM_range:
            f_DM_decay_DM = 1e-27
            if with_run:
                run_all_analysis(run_compare=True, m_DM_range = [m_DM], f_DM_decay_DM_range = [f_DM_decay_DM])
            
            temp = compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
            done_run = False
            while temp[0] < abs(temp[1]/2):
                f_DM_decay_DM_old = f_DM_decay_DM
                if m_DM in [9.5, 10., 13., 16., 19.]:
                    f_DM_decay_DM /= 5
                else:
                    f_DM_decay_DM /= 2
                if with_run:
                    run_all_analysis(run_compare=True, m_DM_range = [m_DM], f_DM_decay_DM_range = [f_DM_decay_DM])
                temp = compare_surveys_DM(m_DM,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
                done_run = True
            if not done_run:
                f_DM_decay_DM_old = f_DM_decay_DM

            theta_min_arr.append(f_DM_decay_DM_old)
         
        if including_ours:

            m_DM_range_ours = use_masses_large + [20.,22.,25.]
            for m_DM in m_DM_range_ours:
                if m_DM >= 20.:
                    f_DM_decay_DM_ours = 5e-27
                else:
                    f_DM_decay_DM_ours = 1e-27
                
                if with_run:
                    run_all_analysis(run_compare=False, m_DM_range = [m_DM],    f_DM_decay_DM_range = [f_DM_decay_DM_ours])
                 
                temp_ours = compare_surveys_DM(m_DM,f_DM_decay_DM_ours, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=False,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
                print(m_DM, f_DM_decay_DM)
                done_run_ours = False
                while temp_ours[0] < abs(temp_ours[1]/2):
                    f_DM_decay_DM_old_ours = f_DM_decay_DM_ours
                    # if m_DM in [ 13., 16., 19.]:
                        # f_DM_decay_DM_ours /= 5
                    # else:
                    f_DM_decay_DM_ours /= 2
                    
                    if with_run:
                        run_all_analysis(run_compare=False, m_DM_range = [m_DM], f_DM_decay_DM_range = [f_DM_decay_DM_ours])
                     
                    temp_ours = compare_surveys_DM(m_DM,f_DM_decay_DM_ours,
                     detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=False,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
                    done_run_ours = True
                if not done_run_ours:
                    f_DM_decay_DM_old_ours = f_DM_decay_DM_ours

                theta_min_arr_ours.append(f_DM_decay_DM_old_ours)

    else:
        if not including_ours:
            f_DM_decay_DM = 1e-26
            temp = compare_surveys_DM(False,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=True,mass_function=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
            done_run = False
            while temp[0] < abs(temp[1]/2):
                f_DM_decay_DM_old = f_DM_decay_DM
                f_DM_decay_DM /= 5
                if with_run:
                    run_all_analysis_spectrum(run_compare=True, m_DM_range = [m_DM], f_DM_decay_DM_range = [f_DM_decay_DM])
                temp = compare_surveys_DM(False,f_DM_decay_DM, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=True,mass_function=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
                done_run = True
            if not done_run:
                f_DM_decay_DM_old = f_DM_decay_DM
            theta_min_arr = f_DM_decay_DM_old
            m_DM_range = False
        else:
            f_DM_decay_DM_ours = 1e-26
            temp_ours = compare_surveys_DM(False,f_DM_decay_DM_ours, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=False,mass_function=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
            done_run_ours = False
            while temp_ours[0] < abs(temp_ours[1]/2):
                f_DM_decay_DM_old_ours = f_DM_decay_DM_ours
                f_DM_decay_DM_ours /= 2.
                if with_run:
                    run_all_analysis_spectrum(run_compare=False, m_DM_range = [m_DM], f_DM_decay_DM_range = [f_DM_decay_DM_ours])
                temp_ours = compare_surveys_DM(False,f_DM_decay_DM_ours, detector = 'GALEX_ULTRASAT_DESIDESI', pars =pars_fid_DM, prior = True, group_vox = True,run_compare=False,mass_function=True,alpha_DDM=alpha_DDM,m_0_DDM=m_0_DDM)
                done_run_ours = True

            if not done_run_ours:
                f_DM_decay_DM_old_ours = f_DM_decay_DM_ours
            theta_min_arr_ours = f_DM_decay_DM_old_ours
            m_DM_range_ours = False
           

    if plot_flag:
        plt.figure(figsize=(13,9))
        fontsize_label = 27
        fontsize_legend = 23
        fontsize_tick = 23

        hirax_vid = np.genfromtxt('DDM_compare/hirax_VID.dat').T
        comap2_vid = np.genfromtxt('DDM_compare/comap2_VID.dat').T
        ccat_vid = np.genfromtxt('DDM_compare/ccat_VID.dat').T
        atlast_vid = np.genfromtxt('DDM_compare/atlast_VID.dat').T
        spherex_vid = np.genfromtxt('DDM_compare/spherex_VID.dat').T
        hetdex_vid = np.genfromtxt('DDM_compare/hetdex_VID.dat').T

        plt.plot(hirax_vid[0],hirax_vid[1],color='navy',linestyle=':')
        plt.plot(comap2_vid[0],comap2_vid[1],color='maroon',linestyle=':')
        plt.plot(ccat_vid[0],ccat_vid[1],color='goldenrod',linestyle=':')
        plt.plot(atlast_vid[0],atlast_vid[1],color='teal',linestyle=':')
        plt.plot(spherex_vid[0],spherex_vid[1],color='yellowgreen',linestyle=':')
        plt.plot(hetdex_vid[0],hetdex_vid[1],color='gray',linestyle=':')

        hera_pk = np.genfromtxt('DDM_compare/HERA_pk.dat').T
        hirax_pk = np.genfromtxt('DDM_compare/hirax_pk.dat').T
        hirax_pk_2 = np.genfromtxt('DDM_compare/hirax_pk_2.dat').T
        comap2_pk = np.genfromtxt('DDM_compare/comap2_pk.dat').T
        comap2_pk_2 = np.genfromtxt('DDM_compare/comap2_pk_2.dat').T
        ccat_pk = np.genfromtxt('DDM_compare/ccat_pk.dat').T
        atlast_pk = np.genfromtxt('DDM_compare/atlast_pk.dat').T
        atlast_pk_2 = np.genfromtxt('DDM_compare/atlast_pk_2.dat').T
        spherex_pk = np.genfromtxt('DDM_compare/spherex_pk.dat').T
        hetdex_pk = np.genfromtxt('DDM_compare/hetdex_pk.dat').T

        plt.plot(hera_pk[0],hera_pk[1],color='plum',linestyle='-',label=r'$\rm HERA$')
        plt.plot(hirax_pk[0],hirax_pk[1],color='navy',linestyle='-',label=r'$\rm HIRAX+CHIME$')
        plt.plot(hirax_pk_2[0],hirax_pk_2[1],color='navy',linestyle='-')
        plt.plot(comap2_pk[0],comap2_pk[1],color='maroon',linestyle='-',label=r'$\rm COMAP2$')
        plt.plot(comap2_pk_2[0],comap2_pk_2[1],color='maroon',linestyle='-')
        plt.plot(ccat_pk[0],ccat_pk[1],color='goldenrod',linestyle='-',label=r'$\rm CCAT-prime$')
        plt.plot(atlast_pk[0],atlast_pk[1],color='teal',linestyle='-',label=r'$\rm AtLAST$')
        plt.plot(atlast_pk_2[0],atlast_pk_2[1],color='teal',linestyle='-')
        plt.plot(spherex_pk[0],spherex_pk[1],color='yellowgreen',linestyle='-',label=r'$\rm SPHEREx$')
        plt.plot(hetdex_pk[0],hetdex_pk[1],color='gray',linestyle='-',label=r'$\rm HETDEX$')

        chime_cross = np.genfromtxt('DDM_compare/chime_cross.dat').T
        comap_cross = np.genfromtxt('DDM_compare/comap_cross.dat').T
        ccat_cross = np.genfromtxt('DDM_compare/ccat_cross.dat').T
        starfire_cross = np.genfromtxt('DDM_compare/starfire_cross.dat').T
        spherex_cross = np.genfromtxt('DDM_compare/spherex_cross.dat').T

        plt.plot(chime_cross[0],chime_cross[1]**-1,color='#fdb515',linestyle='--',label=r'$\rm CHIME\,(cross)$',alpha=0.5,zorder=0)
        plt.plot(comap_cross[0],comap_cross[1]**-1,color='maroon',linestyle='--',label=r'$\rm COMAP\,(cross)$',alpha=0.5,zorder=0)
        plt.plot(ccat_cross[0],ccat_cross[1]**-1,color='k',linestyle='--',label=r'$\rm CCAT\,(cross)$',alpha=0.5,zorder=0)
        plt.plot(starfire_cross[0],starfire_cross[1]**-1,color='mediumseagreen',linestyle='--',label=r'$\rm STARFIRE\,(cross)$',alpha=0.5,zorder=0)
        plt.plot(spherex_cross[0],spherex_cross[1]**-1,color='yellowgreen',linestyle='--',label=r'$\rm SPHEREx\,(cross)$',zorder=0)


        window_size = 7  # Adjust this value based on your preference

        this_work = interp1d(m_DM_range,theta_min_arr)
        this_work_ours = interp1d(m_DM_range_ours,theta_min_arr_ours)
        long_masses = np.linspace(m_DM_range[0],m_DM_range[-1],100)
        long_masses_ours = np.linspace(m_DM_range_ours[0],m_DM_range_ours[-1],100)
        theta_smoothed = moving_average(this_work(long_masses), window_size)
        theta_smoothed_ours = moving_average(this_work_ours(long_masses_ours), window_size)

        plt.fill_between(long_masses[:len(theta_smoothed)],theta_smoothed,1e-23,label=r'$\rm This\, work$',color=color_ULTRASAT,alpha=0.5)
        
        if including_ours:
            plt.plot(m_DM_range_ours,theta_min_arr_ours,linestyle = '--', color=color_ULTRASAT,)

        # nu_CO = 26*u.GHz
        # nu_CO_1 = 34*u.GHz
        # m_DM = (nu_CO * (4*np.pi*cu.hbar)).to(u.eV)
        # m_DM_1 = (nu_CO_1 * (4*np.pi*cu.hbar)).to(u.eV)
        # plt.axvline(m_DM.value,color='k')
        # plt.axvline(m_DM_1.value,color='k')

        plt.xlim(7e-7,30)
        plt.ylim(1e-40,5e-24)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(ncol=2,fontsize=fontsize_legend)
        plt.xticks(fontsize=fontsize_tick)
        plt.yticks(fontsize=fontsize_tick)
        plt.xlabel(r'$m_{\rm DDM}c^2\,[{\rm eV}]$',fontsize=fontsize_label)
        plt.ylabel(r'$\rm min[\tilde{\Theta}_{\rm DDM}]$',fontsize=fontsize_label)

        plt.tight_layout()
        plt.savefig('results/PLOTS/EBL_DM/theta_DDM_compare_full.png',bbox_inches='tight')

#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
        plt.figure(figsize=(9,9))
        fontsize_label = 27
        fontsize_legend = 23
        fontsize_tick = 23

        
        plt.plot((spherex_pk[0]),np.log10(spherex_pk[1]),color='yellowgreen',linestyle='-',label=r'${\rm SPHEREx},\,P(k)$')
        plt.plot((hetdex_pk[0]),np.log10(hetdex_pk[1]),color='teal',linestyle='-',label=r'${\rm HETDEX},\,P(k)$')

        plt.plot((spherex_vid[0]),np.log10(spherex_vid[1]),color='yellowgreen',linestyle=':',label=r'$\rm SPHEREx,\,VID$')
        plt.plot((hetdex_vid[0]),np.log10(hetdex_vid[1]),color='teal',linestyle=':',label=r'$\rm HETDEX,\,VID$')

        #plt.fill_between(m_DM_range,np.log10(theta_min_arr),np.log10(1e-23),label=r'$\rm This\, work,\,\tilde{ \Theta}_{\rm DDM}$',color=color_ULTRASAT,alpha=0.5)
        plt.fill_between((long_masses[:len(theta_smoothed)]),np.log10(theta_smoothed),np.log10(1e-23),label=r'$\rm This\, work,\,\tilde{ \Theta}_{\rm DDM}$',color=color_ULTRASAT,alpha=0.5)

        if including_ours:
            plt.fill_between((long_masses_ours[:len(theta_smoothed_ours)]),np.log10(theta_smoothed_ours),np.log10(1e-23),linestyle = '--', facecolor='none',edgecolor=color_ULTRASAT,label=r'$\rm This\, work,\, \Theta_{\rm DDM}$', linewidth=2, zorder=100)
            #plt.fill_between((m_DM_range_ours),np.log10(theta_min_arr_ours),np.log10(1e-23),linestyle = '--', facecolor='none',edgecolor=color_ULTRASAT,label=r'$\rm This\, work,\, \Theta_{\rm DDM}$', linewidth=2, zorder=100)
    
        plt.xlim(1,30)
        plt.ylim(-33,np.log10(5e-24))
        #plt.xscale('log')
        #plt.yscale('log')
        plt.legend(ncol=2,fontsize=fontsize_legend,loc=8)
        plt.xticks(fontsize=fontsize_tick)
        plt.yticks(fontsize=fontsize_tick)
        plt.xlabel(r'$m_{\rm DDM}c^2\,[{\rm eV}]$',fontsize=fontsize_label)
        plt.ylabel(r'$\log_{10}{\rm min}[\tilde{\Theta}_{\rm DDM}]$',fontsize=fontsize_label)

        plt.tight_layout()
        plt.savefig('results/PLOTS/EBL_DM/theta_DDM_compare.png',bbox_inches='tight')

        plt.show()

    if not including_ours:

        return theta_min_arr, m_DM_range
    else:
        return theta_min_arr_ours, m_DM_range_ours



####################################################
####################################################
####################################################

