from LIM_b7 import *
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

figsize = (12, 8)
fontsize = 30
linewidth = 2.5

plt.rcParams['lines.linewidth'] = linewidth
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize*.8
plt.rcParams['ytick.labelsize'] = fontsize*.8
plt.rcParams['legend.fontsize'] = fontsize*.9
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['legend.columnspacing'] = 0.8

color_ULTRASAT = '#ec4067'
color_NUV = '#468B97'
color_FUV = '#243763'
color_CASTOR = '#7E1717'
colors = [color_FUV, color_NUV, color_ULTRASAT, color_CASTOR]

############################################################

cosmo_input = dict(f_NL=0, H0=67.67,cosmomc_theta=None,
    ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
    num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
    YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
    TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
    theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
    dark_energy_model='ppf',As=2.105e-09, ns=0.967, nrun=0., nrunrun=0.0, r=0.0, nt=None, ntrun=0.0, 
    pivot_scalar=0.05, pivot_tensor=0.05,
    parameterization=2,halofit_version='mead')

nonlinear=True

camb_pars = camb.set_params(\
    H0=cosmo_input['H0'], cosmomc_theta=cosmo_input['cosmomc_theta'],
    ombh2=cosmo_input['ombh2'], omch2=cosmo_input['omch2'], omk=cosmo_input['omk'],
    neutrino_hierarchy=cosmo_input['neutrino_hierarchy'], 
    num_massive_neutrinos=cosmo_input['num_massive_neutrinos'],
    mnu=cosmo_input['mnu'], nnu=cosmo_input['nnu'], YHe=cosmo_input['YHe'], 
    meffsterile=cosmo_input['meffsterile'], 
    standard_neutrino_neff=cosmo_input['standard_neutrino_neff'], 
    TCMB=cosmo_input['TCMB'], tau=cosmo_input['tau'], 
    deltazrei=cosmo_input['deltazrei'], 
    bbn_predictor=cosmo_input['bbn_predictor'], 
    theta_H0_range=cosmo_input['theta_H0_range'],
    w=cosmo_input['w'],  wa=cosmo_input['wa'], cs2=cosmo_input['cs2'], 
    dark_energy_model=cosmo_input['dark_energy_model'],
    As=cosmo_input['As'], ns=cosmo_input['ns'], 
    nrun=cosmo_input['nrun'], nrunrun=cosmo_input['nrunrun'], 
    r=cosmo_input['r'], nt=cosmo_input['nt'], ntrun=cosmo_input['ntrun'], 
    pivot_scalar=cosmo_input['pivot_scalar'], 
    pivot_tensor=cosmo_input['pivot_tensor'],
    parameterization=cosmo_input['parameterization'],
    halofit_version=cosmo_input['halofit_version'])

camb_pars.WantTransfer=True    
camb_pars.Transfer.accurate_massive_neutrinos = True

zmax = 15.  
Nz = 150
zcosmo = np.linspace(0.,zmax,Nz)

camb_pars.set_matter_power(redshifts=list(zcosmo))

cosmo = camb.get_results(camb_pars)

if camb_pars.num_nu_massive != 0:
   var = 8
else:
   var = 7

s8lin = cosmo.get_sigma8()
Dgrowth = interp1d(zcosmo[::-1],s8lin/s8lin[-1],kind='cubic', bounds_error=False,fill_value='extrapolate')

kmin = 1e-3*u.Mpc**-1
kmax = 1e3*u.Mpc**-1
nk = int(1e3)
k_edge = np.linspace(kmin.value,kmax.value,nk)*kmin.unit
Nedge = k_edge.size
k_arr =  (k_edge[0:Nedge-1]+k_edge[1:Nedge])/2.



zmin_gal = 0.1
zmax_gal = 2
R = 350
delta_zi = lambda gal_survey: 1/R if gal_survey == 'SDSS' else 0.01 if  gal_survey == 'SPHEREx' else -1
z_edge = np.arange(zmin_gal,zmax_gal,delta_zi('SDSS'))
Nzedge = z_edge.size
z_SDSS =  (z_edge[0:Nzedge-1]+z_edge[1:Nzedge])/2.

z_edge = np.arange(zmin_gal,zmax_gal,delta_zi('SPHEREx'))
Nzedge = z_edge.size
z_SPHEREx =  (z_edge[0:Nzedge-1]+z_edge[1:Nzedge])/2.

z_edge_small =  np.arange(zmin_gal,zmax_gal,0.05)
Nzedge_small = z_edge_small.size
z_small =  (z_edge_small[0:Nzedge_small-1]+z_edge_small[1:Nzedge_small])/2.

z_edge_castor =  np.arange(zmin_gal,3,0.05)
Nzedge_castor = z_edge_castor.size
z_small_castor =  (z_edge_castor[0:Nzedge_castor-1]+z_edge_castor[1:Nzedge_castor])/2.


PK = camb.get_matter_power_interpolator(camb_pars, zmin=0, zmax=zmax, nz_step=64, 
zs=None, kmax=kmax.value, nonlinear=nonlinear,
var1=var, var2=var, hubble_units=False, 
k_hunit=False, return_z_k=False, k_per_logint=None, log_interp=False,  extrap_kmax=True)



############################################################


wavelenght_min = lambda dtc: 175*u.nm if dtc == 'GALEX_NUV' else 135*u.nm if dtc == 'GALEX_FUV' else 230*u.nm if dtc ==  'ULTRASAT' else 105*u.nm if dtc == 'CASTOR_UV'  else 300*u.nm if dtc == 'CASTOR_U'  else 400*u.nm if dtc == 'CASTOR_G'else -1

wavelenght_max = lambda dtc: 280*u.nm if dtc == 'GALEX_NUV' else 175*u.nm if dtc == 'GALEX_FUV' else 290*u.nm if dtc ==  'ULTRASAT' else 300*u.nm if dtc == 'CASTOR_UV'  else 400*u.nm if dtc == 'CASTOR_U'  else 550*u.nm if dtc == 'CASTOR_G' else -1

lambda_dtc = lambda dtc: (wavelenght_max(dtc)+wavelenght_min(dtc))/2.


# posterior values from tab 1 ref_GAL
fiducials = {
    'bias': [0.32,-0.86,0.79],
    'eps150': [pow(10,25.13)/0.32,2.06],
    'alpha150': [-0.08,1.85],
    'alpha110': [-3.71,0.50],
    'alpha90': -1.5,
    'fescape': [-0.53,-0.84],
    'EW': [-0.617*u.nm,8.802*u.nm]
}

def nu_from_lambda(wavelenght):

    # Wavelenght in nm, frequency in Hz
    nu = (cu.c / (wavelenght.to(u.m))).to(u.Hz)

    return nu

nu_min_gNUV = nu_from_lambda(wavelenght_max('GALEX_NUV'))
nu_max_gNUV = nu_from_lambda(wavelenght_min('GALEX_NUV'))

nu_min_gFUV = nu_from_lambda(wavelenght_max('GALEX_FUV'))
nu_max_gFUV = nu_from_lambda(wavelenght_min('GALEX_FUV'))

nu_min_US = nu_from_lambda(wavelenght_max('ULTRASAT'))
nu_max_US = nu_from_lambda(wavelenght_min('ULTRASAT'))

nu_min_CUV = nu_from_lambda(wavelenght_max('CASTOR_UV'))
nu_max_CUV = nu_from_lambda(wavelenght_min('CASTOR_UV'))

nu_min_CU = nu_from_lambda(wavelenght_max('CASTOR_U'))
nu_max_CU = nu_from_lambda(wavelenght_min('CASTOR_U'))

nu_min_CG = nu_from_lambda(wavelenght_max('CASTOR_G'))
nu_max_CG = nu_from_lambda(wavelenght_min('CASTOR_G'))


def lambda_from_nu(nu):

    # Wavelenght in nm, frequency in Hz
    wavelenght = (cu.c / nu).to(u.nm)

    return wavelenght


def bias(nu, z, vals):

    # eq 15 ref_GAL + eq 24
    # normalization at wavelenght = 150nm, z = 0
    if vals is False:
        b_150_z0, gamma_bnu, gamma_bz = fiducials['bias']
        nu_scale = nu_from_lambda(150*u.nm)
        b_nu_z = b_150_z0 * ((nu/nu_scale)**gamma_bnu) * ((z+1.)**gamma_bz)

    elif type(vals) == float or type(vals) == int:
        b_nu_z = vals
    else:
        b_150_z0, gamma_bnu, gamma_bz = vals    
        nu_scale = nu_from_lambda(150*u.nm)
        b_nu_z = b_150_z0 * ((nu/nu_scale)**gamma_bnu) * ((z+1.)**gamma_bz)

    return b_nu_z


def plot_bias():

    # reproduce fig 11

    waves = [150*u.nm, 300*u.nm, 600*u.nm]

    for i in waves:
        nu_i = nu_from_lambda(i)
        b = bias(nu_i,z_small,False)
        plt.plot(z_small, b, label=r'$\lambda = %g\,{\rm nm}$'%i.value)
    
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b(\nu,z)$',fontsize=fontsize)
    plt.legend(loc=2)

    plt.xlim(zmin_gal,zmax_gal)
    plt.ylim(0,3.5)

    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bias_nuz.png',bbox_inches='tight')

    plt.show()

    return 

################################


def gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def eps150(z, vals = False):

    if vals is False:
        eps_150_z0, gamma_e150 = fiducials['eps150']
        eps_150 = eps_150_z0 * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * (1.+z)**gamma_e150 
    elif type(vals) == float or type(vals) == int:
        eps_150 = vals * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3)
    else:
        eps_150_z0, gamma_e150 = vals
        eps_150 = eps_150_z0 * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * (1.+z)**gamma_e150 

    return eps_150

def alpha150(z, vals = False):

    if vals is False:
        alpha_150_z0, Calpha150 = fiducials['alpha150']
        alpha_150 = alpha_150_z0 + Calpha150*np.log10(1.+z)
    elif type(vals) == float or type(vals) == int:
        alpha_150 = vals    
    else:
        alpha_150_z0, Calpha150 = vals
        alpha_150 = alpha_150_z0 + Calpha150*np.log10(1.+z)

    return alpha_150

def non_ionizing_continuum(nu,z, vals_eps150 = False, vals_alpha150 = False):

    # above the Lya line - 150 A 

    eps_150 = eps150(z, vals_eps150)

    alpha_150 = alpha150(z, vals_alpha150)

    nu_scale = nu_from_lambda(150*u.nm)

    eps = eps_150 * (nu / nu_scale)**alpha_150

    return eps


def alpha110(z, vals = False):

    if vals is False:
        alpha_110_z0, Calpha110 = fiducials['alpha110']
        alpha_110 = alpha_110_z0 + Calpha110*np.log10(1.+z)
    elif type(vals) == float or type(vals) == int:
        alpha_110 = vals
    else:
        alpha_110_z0, Calpha110 = vals
        alpha_110 = alpha_110_z0 + Calpha110*np.log10(1.+z)

    return alpha_110 


def EW_val(z, vals = False):

    if vals is False:
        EW03, EW1 = fiducials['EW']
        CLya = (EW1 - EW03) / np.log10((1+1.)/(1+0.3))
        EW = (CLya * np.log10((1.+z)/(1+0.3)) + EW03) 
    elif type(vals) == float or type(vals) == int:
        EW = vals*u.nm
    else:
        EW03, EW1 = vals
        CLya = (EW1 - EW03) / np.log10((1+1.)/(1+0.3))
        EW = (CLya * np.log10((1.+z)/(1+0.3)) + EW03) 

    return EW

def Lya_line(nu,z, vals_alpha110 = False, val_EW = False, to_plot = False):

    nu_scale = nu_from_lambda(121.6*u.nm)

    # non ionizing continuum 

    alpha_110 = alpha110(z,vals_alpha110)
    eps_Nion_cont = (nu / nu_scale)**alpha_110

    EW = EW_val(z, val_EW) 
    # Lya emission line with equivalent width
    eps_line = nu**2 / (cu.c.to(u.nm/u.s)) * EW / u.Hz

    # line width for illustration purposes
    if to_plot:
        FWHM = 0.5*u.nm
        # delta dirac is a sharp gaussian
        line =  eps_line * gaussian(lambda_from_nu(nu).value, (121.6*u.nm).value, (FWHM/(2 * np.sqrt(2 * np.log(2)))).value)
         #/ u.Hz if (nu_from_lambda(121.6*u.nm + FWHM) <= nu <= nu_from_lambda(121.6*u.nm - FWHM)) else 0. 

    else:
        line =  eps_line  if nu == nu_from_lambda(121.6*u.nm) else 0. 

    eps = eps_Nion_cont + line 
    
    return eps

def fLyC(z, vals):

    # Lyman continuum escape fraction
    if vals is False:
        logf_LyC1, logf_LyC2 = fiducials['fescape']
        CLyC = (logf_LyC2 - logf_LyC1) / (np.log10((1+2.)/(1+1.)))
        f_LyC = pow(10,CLyC * np.log10((1.+z)/(1.+1)) + logf_LyC2)
    elif type(vals) == float or type(vals) == int:
        f_LyC = vals
    else:
        logf_LyC1, logf_LyC2 = vals
        CLyC = (logf_LyC2 - logf_LyC1) / (np.log10((1+2.)/(1+1.)))
        f_LyC = pow(10,CLyC * np.log10((1.+z)/(1.+1)) + logf_LyC2)

    return f_LyC

def Lya_break(nu,z,vals_alpha110 = False, val_flyc = False,val_alpha90 = False):

    # ioinizing continuum at 90 A - not resolved 
    # related with Lya escape fraction 
    nu_scale_121 = nu_from_lambda(121.6*u.nm)
    nu_scale_91 = nu_from_lambda(91.2*u.nm)

    alpha_110 = alpha110(z,vals_alpha110)
    if val_alpha90 is False:
        alpha_90 = fiducials['alpha90']
    else:
        alpha_90 = val_alpha90

    # Lyman continuum escape fraction
    f_LyC = fLyC(z,val_flyc) 

    eps_Nion_cont = ( nu_scale_91 / nu_scale_121 )**alpha_110

    eps_ion_cont = ( nu / nu_scale_91 )**alpha_90

    eps = f_LyC * eps_Nion_cont * eps_ion_cont

    return eps


def signal(wavelenght,z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90, to_plot = False):

    if wavelenght.value <= 91.2:
        eps = non_ionizing_continuum(nu_from_lambda(121.6*u.nm),z,vals_eps150,vals_alpha150) * Lya_break(nu_from_lambda(wavelenght),z,vals_alpha110,val_flyc,val_alpha90)

    elif 91.2 < wavelenght.value <= 121.6:
        eps = non_ionizing_continuum(nu_from_lambda(121.6*u.nm),z,vals_eps150,vals_alpha150) * Lya_line(nu_from_lambda(wavelenght),z,vals_alpha110,val_EW,to_plot)

    else:
        eps = non_ionizing_continuum(nu_from_lambda(wavelenght),z,vals_eps150,vals_alpha150)

    return eps 

def signal_with_DM(wavelenght,z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90, to_plot = False):

    mass_DM_eV = 5 # eV 
    nu0 = 1.21e14 * mass_DM_eV * u.Hz # ---> 495 eV , in ultrasat band from z=1
    gamma = 1e-18 * u.s**-1

    rho_crit = (2.77536627e11*(u.Msun.to(u.kg)*u.kg*u.Mpc**-3) * cu.c**2 ).to(u.erg/u.Mpc**3)     
    rhoM = rho_crit*(camb_pars.omegam-camb_pars.omeganu)

    f_DM = 1
    rho_DM_0 = f_DM * rhoM
    rho_DM = (rho_DM_0*(1+z)**3) / u.Hz #* np.exp(-gamma*(t-t0)) --> exponent around 1 

    eps_DM = (gamma * rho_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * gaussian(nu_from_lambda(wavelenght), nu_from_lambda(121.6*u.nm), nu_from_lambda(wavelenght)**2 / (cu.c.to(u.nm/u.s)) * EW/1e10) #if int(wavelenght.value) == int(lambda_from_nu(nu0).value) else 0.

    if wavelenght.value <= 91.2:
        eps = non_ionizing_continuum(nu_from_lambda(121.6*u.nm),z,vals_eps150,vals_alpha150) * Lya_break(nu_from_lambda(wavelenght),z,vals_alpha110,val_flyc,val_alpha90)

    elif 91.2 < wavelenght.value <= 121.6:
        eps = non_ionizing_continuum(nu_from_lambda(121.6*u.nm),z,vals_eps150,vals_alpha150) * Lya_line(nu_from_lambda(wavelenght),z,vals_alpha110,val_EW,to_plot)

    else:
        eps = non_ionizing_continuum(nu_from_lambda(wavelenght),z,vals_eps150,vals_alpha150)

    return eps + eps_DM


def plot_parameters():

    # reproduce fig 7
 
    z = z_small
    eps_150 = eps150(z)
    #plt.subplot(511)
    plt.plot(z, eps_150, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.ylim(1e25,1e27)
    plt.ylabel(r'$\epsilon_{150}$',fontsize=fontsize)
    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/eps_150.png',bbox_inches='tight')
    plt.show()

    alpha_150 = alpha150(z)
    #plt.subplot(512)
    plt.plot(z, alpha_150, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{150}$',fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/alpha_150.png',bbox_inches='tight')

    plt.show()

    alpha_110 = alpha110(z)
    #plt.subplot(513)
    plt.plot(z, alpha_110, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{110}$',fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/alpha_110.png',bbox_inches='tight')

    plt.show()

    EW = EW_val(z)*10
    #plt.subplot(514)
    plt.plot(z, EW, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.ylim(-100,300)
    plt.ylabel(r'$\rm EW_{Ly\alpha} [A]$',fontsize=fontsize)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/EW.png',bbox_inches='tight')

    plt.show()

    f_LyC = fLyC(z,False)
    #plt.subplot(515)
    plt.plot(z, f_LyC, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.ylim(0.,1)
    plt.ylabel(r'$\rm \rm f_{LyC}$',fontsize=fontsize)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/fLyC.png',bbox_inches='tight')

    plt.show()

    return 



def plot_rest_signal(vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False):

    # reproduce fig 2, top panel and fig 8

    z = 0 
    wave = np.concatenate((np.linspace(70,121.6,2000), np.linspace(121.6,495,2000), np.linspace(495,300,200)))
    s = np.zeros(len(wave))
    for i in range(len(wave)):
        s[i] = signal(wave[i]*u.nm,z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value
        #s[i] = signal(wave[i]*u.nm,z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value

    plt.loglog(wave, s, label=r'$z = %g$'%z)
    #plt.yscale('log')
    plt.xlabel(r'$\lambda_{\rm rest}\,[{\rm nm}]$',fontsize=fontsize)
    plt.ylabel(r'$\epsilon_\nu\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.xlim(wave[0],wave[-1])
    plt.ylim(1e24,1e29)

    plt.show()

    return 


def plot_z_signal(vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False):

    # reproduce fig 2, mid panel and fig 8

    #z = [0,0.3,0.6,1,1.5,2.1]
    z = [0,0.5,1.,2]

    rest_wave = np.concatenate((np.linspace(70,121.6,2000), np.linspace(121.6,350,2000)))

    for zv in range(len(z)):
        s = np.zeros(len(rest_wave))
        wave = np.zeros((len(rest_wave)))
        for i in range(len(rest_wave)):
            wave[i] = rest_wave[i] * (1.+z[zv])
            s[i] = signal(rest_wave[i]*u.nm,z[zv],vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90,True).value 

#        plt.subplot(121)
        plt.plot(wave, s, color= colors[zv], label=r'$z = %g$'%z[zv])
## 
#        plt.subplot(122)
#        plt.loglog(rest_wave, s, color= colors[zv], label=r'$z = %g$'%z[zv])

    # plt.subplot(121)
    plt.yscale('log')
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm nm}]$',fontsize=fontsize)
    plt.ylabel(r'$\epsilon_\nu\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize)
    plt.legend(loc=1)
# 
    #plt.subplot(122)
    # plt.xlabel(r'$\lambda\,[{\rm nm}]$')
    # plt.ylabel(r'$\epsilon_\nu\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$')
    # plt.legend(loc=1)
# 
    plt.xlim(80,800)
    plt.ylim(3e24,3e28)

    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/emissivity.png',bbox_inches='tight')
    plt.show()

    return 


def Response(wavelenght,detector):

    # response function from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=GALEX/GALEX.NUV
    if detector == 'GALEX_NUV':
        wave_dat, R_dat = np.genfromtxt('dat/GALEX.NUV.dat').T
        wave_dat *= 0.1
        area_NUV = np.trapz(R_dat/wave_dat,wave_dat) # !!!! integration with change A -> nm but check final value
        interp = interp1d(wave_dat,R_dat / area_NUV,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.
    elif detector == 'GALEX_FUV':
        wave_dat, R_dat = np.genfromtxt('dat/GALEX.FUV.dat').T
        wave_dat *= 0.1
        area_FUV = np.trapz(R_dat/wave_dat,wave_dat) # !!!! integration with change A -> nm but check final value
        interp = interp1d(wave_dat,R_dat / area_FUV,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.

    elif detector == 'ULTRASAT':
        wave_dat, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T
        detector_area = (4.5*4.5) # cm**2
        Ndetector = 4
        R_dat = QE_dat * detector_area * Ndetector
        area_NUV = np.trapz(R_dat/wave_dat,wave_dat) 
        interp = interp1d(wave_dat,R_dat / area_NUV,fill_value=0.)    

        try:
            R = interp(wavelenght)
        except:
            R = 0.

    elif detector == 'CASTOR_UV':
        wave_dat, R_dat = np.genfromtxt('dat/CASTOR_uv.dat').T
        wave_dat *= 0.1
        area_CUV = np.trapz(R_dat/wave_dat,wave_dat) # !!!! integration with change A -> nm but check final value
        interp = interp1d(wave_dat,R_dat / area_CUV,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.

    elif detector == 'CASTOR_U':
        wave_dat, R_dat = np.genfromtxt('dat/CASTOR_u.dat').T
        wave_dat *= 0.1
        area_CU = np.trapz(R_dat/wave_dat,wave_dat) # !!!! integration with change A -> nm but check final value
        interp = interp1d(wave_dat,R_dat / area_CU,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.

    elif detector == 'CASTOR_G':
        wave_dat, R_dat = np.genfromtxt('dat/CASTOR_g.dat').T
        wave_dat *= 0.1
        area_CG = np.trapz(R_dat/wave_dat,wave_dat) # !!!! integration with change A -> nm but check final value
        interp = interp1d(wave_dat,R_dat / area_CG,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.
    else:
        R = 1.

    return R



def plot_response():

    wavelenght, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T
    
    plt.figure()
    plt.plot(wavelenght,QE_dat,label=r'$ULTRASAT$',color=color_ULTRASAT)
    plt.legend(loc=1)
    plt.xlabel(r'$\lambda\,[{\rm nm}]$',fontsize=fontsize)
    plt.ylabel(r'$\rm QE$',fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/QE_ultrasat.png',bbox_inches='tight')

    wavelenght = np.linspace(wavelenght_min('ULTRASAT'),wavelenght_max('ULTRASAT'))

    wavelenght_N = np.linspace(wavelenght_min('GALEX_NUV'),wavelenght_max('GALEX_NUV'))

    wavelenght_F = np.linspace(wavelenght_min('GALEX_FUV'),wavelenght_max('GALEX_FUV'))

    wavelenght_CUV = np.linspace(wavelenght_min('CASTOR_UV'),wavelenght_max('CASTOR_UV'))

    wavelenght_CU = np.linspace(wavelenght_min('CASTOR_U'),wavelenght_max('CASTOR_U'))

    wavelenght_CG = np.linspace(wavelenght_min('CASTOR_G'),wavelenght_max('CASTOR_G'))

 
    R = np.zeros(len(wavelenght))
    RN = np.zeros(len(wavelenght))
    RF = np.zeros(len(wavelenght))
    RCUV = np.zeros(len(wavelenght))
    RCU = np.zeros(len(wavelenght))
    RCG = np.zeros(len(wavelenght))

    for i in range(len(wavelenght)):
        R[i] = Response(wavelenght[i],'ULTRASAT')
        RN[i] = Response(wavelenght_N[i],'GALEX_NUV')
        RF[i] = Response(wavelenght_F[i],'GALEX_FUV')
        RCUV[i] = Response(wavelenght_CUV[i],'CASTOR_UV')
        RCU[i] = Response(wavelenght_CU[i],'CASTOR_U')
        RCG[i] = Response(wavelenght_CG[i],'CASTOR_G')

    plt.figure()
    plt.plot(wavelenght,R,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(wavelenght_N,RN,label=r'$\rm NUV$',color=color_NUV)
    plt.plot(wavelenght_F,RF,label=r'$\rm FUV$',color=color_FUV)
    plt.plot(wavelenght_CUV,RCUV,label=r'$uv$',color=color_CASTOR,linestyle=':')
    plt.plot(wavelenght_CU,RCU,label=r'$u$',color=color_CASTOR,linestyle='--')
    plt.plot(wavelenght_CG,RCG,label=r'$g$',color=color_CASTOR,linestyle='-.')

    plt.legend(loc=1,ncol=2)
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm nm}]$',fontsize=fontsize)
    plt.ylabel(r'$R(\lambda_{\rm obs})$',fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/response.png',bbox_inches='tight')

    plt.show()

    return 




def H(z):
        
    # Hubble parameter from camb
    
    H = cosmo.hubble_parameter(z)*(u.km/u.Mpc/u.s)

    return H


def tau_LS(wavelengh,z):

    # Lyman series abs - sec 4.1 : lyA forest + damped lyA systems

    # Ly series absorption
    lines_j, lambda_j, A_LAF_j1, A_LAF_j2, A_LAF_j3, A_DLA_j1, A_DLA_j2 = np.genfromtxt('dat/Ly_series_absorption.txt').T

    tau_LAF_LS = 0.
    tau_DLA_LS = 0.
    for j in range(len(lines_j)):

        use_lambda_j = lambda_j[j] * 0.1 * u.nm

        val = 0.
        if use_lambda_j < wavelengh < use_lambda_j*(1+z):
            if wavelengh < 2.2*use_lambda_j:     
                val = A_LAF_j1[j]*(wavelengh/use_lambda_j)**1.2
            elif 2.2*use_lambda_j <= wavelengh < 5.7*use_lambda_j:         
                val = A_LAF_j2[j]*(wavelengh/use_lambda_j)**3.7
            else:
                val = A_LAF_j3[j]*(wavelengh/use_lambda_j)**5.5
        else:
            val = 0.

        tau_LAF_LS += val

    for j in range(len(lines_j)):

        use_lambda_j = lambda_j[j] * 0.1 * u.nm

        val = 0.
        if use_lambda_j < wavelengh < use_lambda_j*(1+z):
            if wavelengh < 3*use_lambda_j:         
                val = A_DLA_j1[j]*(wavelengh/use_lambda_j)**2.
            else:
                val = A_DLA_j2[j]*(wavelengh/use_lambda_j)**3.
        else:
            val = 0.

        tau_DLA_LS += val

    tau_LS = tau_LAF_LS + tau_DLA_LS

    return tau_LS


def tau_LAF_LC(wavelengh, z):

    tau_LAF_LC = 0.
    wave_L = 91.18*u.nm # Ly limit

    if wavelengh > wave_L:
        if z < 1.2:
            if wavelengh < wave_L*(1+z):
                tau_LAF_LC = 0.325 * ((wavelengh/wave_L)**1.2 - (1+z)**-0.9 * (wavelengh / wave_L)**2.1)
            else:
                tau_LAF_LC = 0.

        elif 1.2 <= z < 4.7:
            if wavelengh < 2.2*wave_L:
                tau_LAF_LC = 2.55e-2*(1+z)**1.6*(wavelengh / wave_L)**2.1 + 0.325*(wavelengh / wave_L)**1.2 -0.250*(wavelengh / wave_L)**2.1
            elif 2.2*wave_L <= wavelengh < (1+z)*wave_L:
                tau_LAF_LC = 2.55e-2*((1+z)**1.6*(wavelengh / wave_L)**2.1 -(wavelengh / wave_L)**3.7)
            else:
                tau_LAF_LC = 0.

        else:
            if wavelengh < 2.2*wave_L:
                tau_LAF_LC = 5.22e-4*(1+z)**3.4*(wavelengh/wave_L)**2.1 +  0.325 * (wavelengh/wave_L)**1.2 - 3.14e-2 * (wavelengh / wave_L)**2.1
            elif 2.2*wave_L <= wavelengh < 5.7*wave_L:
                tau_LAF_LC = 5.22e-4*(1+z)**3.4*(wavelengh/wave_L)**2.1 +  0.218 * (wavelengh/wave_L)**2.1 - 2.55e-2 * (wavelengh / wave_L)**3.7
            elif 5.7*wave_L <= wavelengh < (1+z)*wave_L:
                tau_LAF_LC = 5.22e-4*((1+z)**3.4*(wavelengh/wave_L)**2.1 - (wavelengh / wave_L)**5.5)
            else:
                tau_LAF_LC = 0.
    else:
        tau_LAF_LC = 0.

    return tau_LAF_LC

def tau_DLA_LC(wavelengh, z):

    tau_DLA_LC = 0.
    wave_L = 91.18*u.nm # Ly limit

    if wavelengh > wave_L:
        if z < 2.:
            if wavelengh < wave_L*(1+z):
                tau_DLA_LC = 0.211*(1+z)**2 - 7.66e-2*(1+z)**2.3*(wavelengh/wave_L)**-0.3 -0.135*(wavelengh/wave_L)**2
            else:
                tau_DLA_LC = 0.
        else:
            if wavelengh < 3*wave_L:
                tau_DLA_LC = 0.634 + 4.7e-2*(1+z)**3 - 1.78e-2*(1+z)**3.3*(wavelengh/wave_L)**-0.3-0.135*(wavelengh/wave_L)**2-0.291*(wavelengh/wave_L)**-0.3
            elif 3*wave_L <= wavelengh < (1+z)*wave_L:
                tau_DLA_LC = 4.7e-2*(1+z)**3 - 1.78e-2*(1+z)**3.3*(wavelengh/wave_L)**-0.3
            else:
                tau_DLA_LC = 0.
    else:
        tau_DLA_LC = 0.

    return tau_DLA_LC

def tau_LC(wavelengh, z):

    # Lyman continuum - sec 4.2 : lyA forest + damped lyA systems

    tau_LC = tau_LAF_LC(wavelengh, z) + tau_DLA_LC(wavelengh, z)

    return tau_LC


def tau_Lya(wavelengh, z):

    # from sec 4 in 1402.0677

    tau = tau_LS(wavelengh,z) + tau_LC(wavelengh, z)

    return tau

def plot_tau():

    wavelenghtN, R_dat = np.genfromtxt('dat/GALEX.NUV.dat').T
    wavelenghtF, R_datF = np.genfromtxt('dat/GALEX.FUV.dat').T

    wavelenghtN *= 0.1
    wavelenghtF *= 0.1

    waves = np.linspace(wavelenght_min('GALEX_FUV'),wavelenght_max('GALEX_NUV'))


    z = [0.,0.5,1.,2]
    lines = ['-','--',':']
    for i in z:

        tau_ls = np.zeros(len(waves))
        tau_lc = np.zeros(len(waves))
        tau_tot = np.zeros(len(waves))
        for l in range(len(waves)):
            tau_ls[l] = tau_LS(waves[l],i)
            tau_lc[l] = tau_LC(waves[l],i)
            tau_tot[l] = tau_Lya(waves[l],i)

        #plt.plot(waves,tau_tot,color=colors[z.index(i)],linestyle=lines[0],label=r'$z=%g$'%i)
        if z.index(i) == 0:
          plt.plot(waves,tau_lc,color=colors[z.index(i)],linestyle='-',#lines[1],
          label=r'${\rm Ly\alpha\,\, continuum}$')
          plt.plot(waves,tau_ls,color=colors[z.index(i)],linestyle=lines[1],label=r'${\rm Ly\alpha\,\, series}$')
        else:
            plt.plot(waves,tau_lc,color=colors[z.index(i)],linestyle='-',label=r'$z=%g$'%i)
            plt.plot(waves,tau_ls,color=colors[z.index(i)],linestyle=lines[1])
# 
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm nm}]$',fontsize=fontsize)
    plt.ylabel(r'$\tau(\lambda_{\rm obs},z)$',fontsize=fontsize)

    plt.legend(ncol=2,loc=1)
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/tau.png',bbox_inches='tight')

    plt.show()

    return 


def dJdz(z, detector,run = False,\
         vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False, filename = ''):

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

        elif detector == 'CASTOR_UV':
            nu_min = nu_min_CUV
            nu_max = nu_max_CUV
        elif detector == 'CASTOR_U':
            nu_min = nu_min_CU
            nu_max = nu_max_CU
        elif detector == 'CASTOR_G':
            nu_min = nu_min_CG
            nu_max = nu_max_CG
        else:
            print('Check detector!')
            return 

        unit = (Response(lambda_from_nu(1.*u.Hz), detector) *signal(lambda_from_nu(1*u.Hz),0.,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90) * np.exp(-tau_Lya(lambda_from_nu(1*u.Hz),0.))/(1*u.Hz)).unit

        intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs # !!!!

        #intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        #rest_wave = np.concatenate((np.linspace(lambda_from_nu(nu_min).value,121.6,2000), np.linspace(121.6,lambda_from_nu(nu_max).value,2000)))
        rest_wave = np.concatenate((np.linspace(lambda_from_nu(nu_min).value,121.6,200), np.linspace(121.6,495,200), np.linspace(495,lambda_from_nu(nu_max).value,200)))

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.nm).value
            intg_arr[i] = intg(nu_obs_arr[i])

        #dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * quad(intg,nu_min.value,nu_max.value)[0]*(unit*u.Hz/u.steradian) 
        dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * np.trapz(intg_arr,nu_obs_arr)*(unit*u.Hz/u.steradian) 

    else:
        zval, dJdzval = np.genfromtxt(filename)
        dJdz = interp1d(zval,dJdzval)(z) * u.Jy/u.steradian  

    return (dJdz.to(u.Jy/u.steradian)).value



def bJ_z(z, detector, run = False, vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias=False, filename = ''):

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
        elif detector == 'CASTOR_UV':
            nu_min = nu_min_CUV
            nu_max = nu_max_CUV
        elif detector == 'CASTOR_U':
            nu_min = nu_min_CU
            nu_max = nu_max_CU
        elif detector == 'CASTOR_G':
            nu_min = nu_min_CG
            nu_max = nu_max_CG
        else:
            print('Check detector!')
            return 

        #intg_num = lambda nu_obs: bias(nu_obs*u.Hz*(1+z),z,val_bias) * Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        #intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
#
        intg_num = lambda nu_obs: bias(nu_obs*u.Hz*(1+z),z,val_bias) * Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs # !!! 
#
        intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs # !!!
#
  #      rest_wave = np.concatenate((np.linspace(lambda_from_nu(nu_min).value,121.6,2000), np.linspace(121.6,lambda_from_nu(nu_max).value,2000)))
        rest_wave = np.concatenate((np.linspace(lambda_from_nu(nu_min).value,121.6,200), np.linspace(121.6,495,200), np.linspace(495,lambda_from_nu(nu_max).value,200)))

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        intg_den_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.nm).value
            intg_arr[i] = intg_num(nu_obs_arr[i])
            intg_den_arr[i] = intg_den(nu_obs_arr[i])

        # num = quad(intg_num,nu_min.value,nu_max.value)[0]
        # den = quad(intg_den,nu_min.value,nu_max.value)[0]
        num = np.trapz(intg_arr,nu_obs_arr)
        den = np.trapz(intg_den_arr,nu_obs_arr)

        bJ = num / den

    else:
        zval, bJval = np.genfromtxt(filename)
        bJ = interp1d(zval,bJval)(z) 

    return bJ


def compute_vals(vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False):


    with Pool(6) as pool:

        #print('Doing NUV')
        #dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
        #dJ_nuv = pool.map(dJdz_NUV, z_small)
        #print('Doing FUV')
        #dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
        #dJ_fuv = pool.map(dJdz_FUV, z_small)
 
       # print('Doing b NUV')
       # bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_NUV.dat')
       # bJ_nuv = pool.map(bJ_nuv_f, z_small)
 #
       # print('Doing b FUV')
       # bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_FUV.dat')
       # bJ_fuv = pool.map(bJ_fuv_f, z_small)


        print('Doing ULTRASAT')
        dJdz_f = partial(dJdz,detector='ULTRASAT',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
        dJ_U = pool.map(dJdz_f, z_small)

        print('Doing b ULTRASAT')
        bJ_f =  partial(bJ_z,detector='ULTRASAT',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False)
        bJ_U = pool.map(bJ_f, z_small)
    
      #  print('Doing CASTOR UV')
      #  dJdz_f = partial(dJdz,detector='CASTOR_UV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
      #  dJ_CUV = pool.map(dJdz_f, z_small_castor)
 #
      #  print('Doing b CASTOR UV')
      #  bJ_f =  partial(bJ_z,detector='CASTOR_UV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False)
      #  bJ_CUV = pool.map(bJ_f, z_small_castor)
    #
      #  print('Doing CASTOR U')
      #  dJdz_f = partial(dJdz,detector='CASTOR_U',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
      #  dJ_CU = pool.map(dJdz_f, z_small_castor)
 #
      #  print('Doing b CASTOR U')
      #  bJ_f =  partial(bJ_z,detector='CASTOR_U',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False)
      #  bJ_CU = pool.map(bJ_f, z_small_castor)
 #
      #  print('Doing CASTOR G')
      #  dJdz_f = partial(dJdz,detector='CASTOR_G',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90)
      #  dJ_CG = pool.map(dJdz_f, z_small_castor)
 #
      #  print('Doing b CASTOR G')
      #  bJ_f =  partial(bJ_z,detector='CASTOR_G',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False)
      #  bJ_CG = pool.map(bJ_f, z_small_castor)
    #
#
#    np.savetxt('dat/dJdz_GALEX_NUV.dat',(z_small,np.asarray(dJ_nuv)))
#    np.savetxt('dat/bJ_GALEX_NUV.dat',(z_small,np.asarray(bJ_nuv)))
#
#    np.savetxt('dat/dJdz_GALEX_FUV.dat',(z_small,np.asarray(dJ_fuv)))
  #  np.savetxt('dat/bJ_GALEX_FUV.dat',(z_small,np.asarray(bJ_fuv)))

    np.savetxt('dat/dJdz_ULTRASAT_TESTGAUS.dat',(z_small,np.asarray(dJ_U)))
    np.savetxt('dat/bJ_ULTRASAT_TESTGAUS.dat',(z_small,np.asarray(bJ_U)))

 #   np.savetxt('dat/dJdz_CASTOR_UV.dat',(z_small_castor,np.asarray(dJ_CUV)))
 #   np.savetxt('dat/bJ_CASTOR_UV.dat',(z_small_castor,np.asarray(bJ_CUV)))
#
 #   np.savetxt('dat/dJdz_CASTOR_U.dat',(z_small_castor,np.asarray(dJ_CU)))
 #   np.savetxt('dat/bJ_CASTOR_U.dat',(z_small_castor,np.asarray(bJ_CU)))
#
 #   np.savetxt('dat/dJdz_CASTOR_G.dat',(z_small_castor,np.asarray(dJ_CG)))
 #   np.savetxt('dat/bJ_CASTOR_G.dat',(z_small_castor,np.asarray(bJ_CG)))
#
    return



def plot_signal(vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias=False):

    dJ_U =  dJdz(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_ULTRASAT.dat')

    dJ_N =  dJdz(z_small,detector='GALEX_NUV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_GALEX_NUV.dat')

    dJ_F =  dJdz(z_small,detector='GALEX_FUV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_GALEX_FUV.dat')

    dJ_CUV =  dJdz(z_small_castor,detector='CASTOR_UV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_CASTOR_UV.dat')

    dJ_CU =  dJdz(z_small_castor,detector='CASTOR_U',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_CASTOR_U.dat')

    dJ_CG =  dJdz(z_small_castor,detector='CASTOR_G',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_CASTOR_G.dat')


    bJ_U = bJ_z(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_ULTRASAT.dat')

    bJ_N = bJ_z(z_small,detector='GALEX_NUV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_GALEX_NUV.dat')

    bJ_F = bJ_z(z_small,detector='GALEX_FUV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_GALEX_FUV.dat')

    bJ_CUV = bJ_z(z_small_castor,detector='CASTOR_UV',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_CASTOR_UV.dat')

    bJ_CU = bJ_z(z_small_castor,detector='CASTOR_U',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_CASTOR_U.dat')

    bJ_CG = bJ_z(z_small_castor,detector='CASTOR_G',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_CASTOR_G.dat')

    # figure 10 
    plt.figure()
    plt.plot(z_small,dJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(z_small,dJ_N,label=r'$\rm NUV$',color=color_NUV)
    plt.plot(z_small,dJ_F,label=r'$\rm FUV$',color=color_FUV)
    plt.plot(z_small_castor,dJ_CUV,label=r'$uv$',color=color_CASTOR,linestyle=':')
    plt.plot(z_small_castor,dJ_CU,label=r'$u$',color=color_CASTOR,linestyle='--')
    plt.plot(z_small_castor,dJ_CG,label=r'$g$',color=color_CASTOR,linestyle='-.')
    plt.xlim(z_small[0],z_small_castor[-1])
    #plt.ylim(-10,200)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$dJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/dJdz.png',bbox_inches='tight')

    plt.figure()
    plt.plot(z_small,bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(z_small,bJ_N,label=r'$\rm NUV$',color=color_NUV)
    plt.plot(z_small,bJ_F,label=r'$\rm FUV$',color=color_FUV)
    plt.plot(z_small_castor,bJ_CUV,label=r'$uv$',color=color_CASTOR,linestyle=':')
    plt.plot(z_small_castor,bJ_CU,label=r'$u$',color=color_CASTOR,linestyle='--')
    plt.plot(z_small_castor,bJ_CG,label=r'$g$',color=color_CASTOR,linestyle='-.')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_J$',fontsize=fontsize)
    plt.xlim(0,3)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bJz.png',bbox_inches='tight')

    plt.figure()
    plt.plot(z_small,dJ_U*bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(z_small,dJ_N*bJ_N,label=r'$\rm NUV$',color=color_NUV)
    plt.plot(z_small,dJ_F*bJ_F,label=r'$\rm FUV$',color=color_FUV)
    plt.plot(z_small_castor,dJ_CUV*bJ_CUV,label=r'$\rm uv$',color=color_CASTOR,linestyle=':')
    plt.plot(z_small_castor,dJ_CU*bJ_CU,label=r'$\rm u$',color=color_CASTOR,linestyle='--')
    plt.plot(z_small_castor,dJ_CG*bJ_CG,label=r'$\rm g$',color=color_CASTOR,linestyle='-.')
    plt.xlim(z_small[0],z_small_castor[-1])
    #plt.ylim(-10,200)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_JdJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bJdJdz.png',bbox_inches='tight')


    plt.show()
    return 


def plot_bdJdz_multi():

    # reproduce fig 3
    #alpha_1500 = [-1.,-0.5,0,0.5,1.,False]
    alpha_1500 = [-1.,0,1.,False]

    #EW = [-5,-2.5,0,2.5,5,False]
    EW = [-5,0,5,False]
    
    #f_Lyc = [0.,0.03,0.1,0.3,1,False]
    f_Lyc = [0.,0.1,1,False]
    #color= ['b','c','k','orange','r','k']
    color= ['b','k','r','k']

    vals_eps150= False
    vals_alpha110= False
    val_alpha90= False
    
    plt.figure()
    val_EW= False
    val_flyc= False
 
    print('Doing alpha150')
    for A in tqdm(alpha_1500):
 
        vals_alpha150= A
 
        with Pool(6) as pool:

            print('Doing NUV')
            dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_nuv = pool.map(dJdz_NUV, z_small)
            print('Doing FUV')
            dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_fuv = pool.map(dJdz_FUV, z_small)
    
            bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_NUV.dat')
            bJ_nuv = pool.map(bJ_nuv_f, z_small)
 
            print('Doing b FUV')
            bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_FUV.dat')
            bJ_fuv = pool.map(bJ_fuv_f, z_small)

            bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
            bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

        if not vals_alpha150:
            linestyle = '--'
        else: 
            linestyle = '-'

        plt.subplot(121)
        # plt.subplot(231)
        plt.plot(z_small,bdJ_fuv,color= color[alpha_1500.index(A)],label=r'$\alpha_{150} = %g$'%A,linestyle=linestyle)
 
        # plt.subplot(234)
        plt.subplot(122)
        plt.plot(z_small,bdJ_nuv,color= color[alpha_1500.index(A)],label=r'$\alpha_{150} = %g$'%A,linestyle=linestyle)
 
    # plt.subplot(231)
    plt.subplot(121)
    plt.xlim(z_small[0],z_small[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)
 
    plt.subplot(122)
    # plt.subplot(234)
    plt.xlim(z_small[0],z_small[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)


    plt.figure()

    vals_alpha150= False
    val_flyc= False

    print('Doing EW')
    for A in tqdm(EW):
        val_EW = A

        with Pool(6) as pool:

            print('Doing NUV')
            dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_nuv = pool.map(dJdz_NUV, z_small)
            print('Doing FUV')
            dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_fuv = pool.map(dJdz_FUV, z_small)
    
            bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_NUV.dat')
            bJ_nuv = pool.map(bJ_nuv_f, z_small)
 
            print('Doing b FUV')
            bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_FUV.dat')
            bJ_fuv = pool.map(bJ_fuv_f, z_small)

            bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
            bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

        if not val_EW:
            linestyle = '--'
        else: 
            linestyle = '-'
        # plt.subplot(232)
        plt.subplot(121)
        plt.plot(z_small,bdJ_fuv,color= color[EW.index(A)],label=r'$\rm EW = %g$'%A,linestyle=linestyle)

        plt.subplot(122)
        # plt.subplot(235)
        plt.plot(z_small,bdJ_nuv,color= color[EW.index(A)],label=r'$\rm EW = %g$'%A,linestyle=linestyle)

    plt.subplot(121)
    # plt.subplot(232)
    plt.xlim(z_small[0],z_small[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.subplot(122)
    # plt.subplot(235)
    plt.xlim(z_small[0],z_small[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)


    plt.figure()
    val_EW= False
    vals_alpha150= False

    print('Doing fLyC')
    for A in tqdm(f_Lyc):

        val_flyc= A
        with Pool(6) as pool:

            print('Doing NUV')
            dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_nuv = pool.map(dJdz_NUV, z_small)
            print('Doing FUV')
            dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,)
            dJ_fuv = pool.map(dJdz_FUV, z_small)
    
            bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_NUV.dat')
            bJ_nuv = pool.map(bJ_nuv_f, z_small)
 
            print('Doing b FUV')
            bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = False,filename='dat/dJdz_GALEX_FUV.dat')
            bJ_fuv = pool.map(bJ_fuv_f, z_small)

            bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
            bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

        if not val_flyc:
            linestyle = '--'
        else: 
            linestyle = '-'

        plt.subplot(121)
        # plt.subplot(233)
        plt.plot(z_small,bdJ_fuv,color= color[f_Lyc.index(A)],label=r'$f_{\rm LyC} = %g$'%A,linestyle=linestyle)

        plt.subplot(122)
        # plt.subplot(236)
        plt.plot(z_small,bdJ_nuv,color= color[f_Lyc.index(A)],label=r'$f_{\rm LyC} = %g$'%A,linestyle=linestyle)

    plt.subplot(121)
    # plt.subplot(233)
    plt.xlim(z_small[0],z_small[-1])
#    plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)

    # plt.subplot(236)
    plt.subplot(122)
    plt.xlim(z_small[0],z_small[-1])
#    plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.show()
    return 





def plot_signal_withDM(vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias=False):

    dJ_U =  dJdz(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_ULTRASAT_TESTGAUS.dat')

    
    dJ_U_DM =  dJdz(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,filename='dat/dJdz_ULTRASAT_withDM.dat')


    bJ_U = bJ_z(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_ULTRASAT_TESTGAUS.dat')

    bJ_U_DM = bJ_z(z_small,detector='ULTRASAT',run=False,vals_eps150=vals_eps150,vals_alpha150=vals_alpha150,vals_alpha110=vals_alpha110,val_EW=val_EW,val_flyc=val_flyc,val_alpha90=val_alpha90,val_bias = val_bias,filename='dat/bJ_ULTRASAT_withDM.dat')

    # figure 10 
    plt.figure()
    plt.plot(z_small,dJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    #plt.plot(z_small,dJ_U_DM,label=r'$\rm with\, DM$',color='k')
    plt.xlim(z_small[0],z_small_castor[-1])
    #plt.ylim(-10,200)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$dJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/dJdz_withDM.png',bbox_inches='tight')

    plt.figure()
    plt.plot(z_small,bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(z_small,bJ_U_DM,label=r'$\rm with\,DM$',color='k')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_J$',fontsize=fontsize)
    plt.xlim(0,3)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bJz_withDM.png',bbox_inches='tight')

    plt.figure()
    plt.plot(z_small,dJ_U*bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(z_small,dJ_U_DM*bJ_U_DM,label=r'$\rm with\, DM$',color='k')
    plt.xlim(z_small[0],z_small_castor[-1])
    #plt.ylim(-10,200)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_JdJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
    plt.legend(loc=4, ncol=2)
    
    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bJdJdz_withDM.png',bbox_inches='tight')


    plt.show()
    return 
