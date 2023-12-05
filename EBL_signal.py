from LIM_b7 import *
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial


plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica"
})

figsize = (10, 9)
fontsize = 33
linewidth = 2.5

plt.rcParams['lines.linewidth'] = linewidth
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize*.5
plt.rcParams['ytick.labelsize'] = fontsize*.5
plt.rcParams['legend.fontsize'] = fontsize*.7
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['legend.columnspacing'] = 0.8
plt.rcParams['legend.frameon'] = False

color_FUV = '#1e6091'
color_NUV = '#89c2d9'
color_ULTRASAT = '#eb4464'
color_hetdex = '#ffb703'
color_spherex = '#da2c38'
colors = ['k',color_FUV, color_NUV, color_ULTRASAT, color_hetdex, color_spherex]

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

R_SDSS = 350

R_DESI = 2000 # 2306.06307 sec 2.2.2

delta_zi = lambda gal_survey: 1/R_SDSS if gal_survey == 'SDSS' else 0.01 if  gal_survey == 'SPHEREx' else 1/R_DESI if gal_survey == 'DESI' else -1

z_edge_SDSS = np.arange(zmin_gal,zmax_gal,delta_zi('SDSS'))
Nzedge_SDSS = z_edge_SDSS.size

z_edge_SPHEREx = np.arange(zmin_gal,zmax_gal,delta_zi('SPHEREx'))
Nzedge_SPHEREx = z_edge_SPHEREx.size

z_edge_DESI = np.arange(zmin_gal,zmax_gal,delta_zi('DESI'))
# ALTERNATIVE: np.concatenate((np.arange(zmin_gal,1.6,delta_zi('DESI')),np.arange(1.6,zmax_gal,7*delta_zi('DESI')))) # see tab 5 in 2306.06307 - desi z ucertainty is 0.0005(1+z) below z = 1.6, 0.0025 up to 2.1. Since we have 0.0005*2.6 = 0.001 and 0.0025*3 = 0.007, we take the conservative choice dz = 0.001 for z < 1.6 - which is same order as SDSS - and 0.007 above 
Nzedge_DESI = z_edge_DESI.size

z_gals =  lambda gal_survey: (z_edge_SDSS[0:Nzedge_SDSS-1]+z_edge_SDSS[1:Nzedge_SDSS])/2. if gal_survey == 'SDSS' else (z_edge_SPHEREx[0:Nzedge_SPHEREx-1]+z_edge_SPHEREx[1:Nzedge_SPHEREx])/2. if gal_survey == 'SPHEREx' else (z_edge_DESI[0:Nzedge_DESI-1]+z_edge_DESI[1:Nzedge_DESI])/2. if gal_survey == 'DESI' else -1

delta_zi_interp = 0.1
z_gals_interp = np.arange(zmin_gal,zmax_gal+delta_zi_interp,delta_zi_interp)


#z_edge = np.arange(zmin_gal,zmax_gal,delta_zi('SPHEREx'))
#Nzedge = z_edge.size
#z_SPHEREx =  (z_edge[0:Nzedge-1]+z_edge[1:Nzedge])/2.


PK = camb.get_matter_power_interpolator(camb_pars, zmin=0, zmax=zmax, nz_step=64, 
zs=None, kmax=kmax.value, nonlinear=nonlinear,
var1=var, var2=var, hubble_units=False, 
k_hunit=False, return_z_k=False, k_per_logint=None, log_interp=False,  extrap_kmax=True)



def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')


############################################################


wavelenght_min = lambda dtc: 1750*u.AA if dtc == 'GALEX_NUV' else 1350*u.AA if dtc == 'GALEX_FUV' else 2300*u.AA if dtc ==  'ULTRASAT' else -1
wavelenght_max = lambda dtc: 2800*u.AA if dtc == 'GALEX_NUV' else 1750*u.AA if dtc == 'GALEX_FUV' else 2900*u.AA if dtc ==  'ULTRASAT' else -1

lambda_dtc = lambda dtc: (wavelenght_max(dtc)+wavelenght_min(dtc))/2.


# posterior values from tab 1 ref_GAL
fiducials = {
    'bias': [0.32,-0.86,0.79],
    'eps1500': [pow(10,25.13)/0.32,2.06],
    'alpha1500': [-0.08,1.85],
    'alpha1100': [-3.71,0.50],
    'alpha900': -1.5,
    'fescape': [-0.53,-0.84],
    'EW': lambda detector: [-6.17*u.AA,88.02*u.AA] if detector != 'ULTRASAT' else [88.02*u.AA,176.67*u.AA]
}

def nu_from_lambda(wavelenght):

    # Wavelenght in Angstrom, frequency in Hz
    nu = (cu.c / (wavelenght.to(u.m))).to(u.Hz)

    return nu

nu_min_gNUV = nu_from_lambda(wavelenght_max('GALEX_NUV'))
nu_max_gNUV = nu_from_lambda(wavelenght_min('GALEX_NUV'))

nu_min_gFUV = nu_from_lambda(wavelenght_max('GALEX_FUV'))
nu_max_gFUV = nu_from_lambda(wavelenght_min('GALEX_FUV'))

nu_min_US = nu_from_lambda(wavelenght_max('ULTRASAT'))
nu_max_US = nu_from_lambda(wavelenght_min('ULTRASAT'))


def lambda_from_nu(nu):

    # Wavelenght in angstrom, frequency in Hz
    wavelenght = (cu.c / nu).to(u.AA)

    return wavelenght


def bias(nu, z, vals):

    # eq 15 ref_GAL + eq 24
    # normalization at wavelenght = 1500A, z = 0
    if vals is False:
        b_1500_z0, gamma_bnu, gamma_bz = fiducials['bias']
        nu_scale = nu_from_lambda(1500*u.AA)
        b_nu_z = b_1500_z0 * ((nu/nu_scale)**gamma_bnu) * ((z+1.)**gamma_bz)

    elif type(vals) == float or type(vals) == int:
        b_nu_z = vals
    else:
        b_1500_z0, gamma_bnu, gamma_bz = vals    
        nu_scale = nu_from_lambda(1500*u.AA)
        b_nu_z = b_1500_z0 * ((nu/nu_scale)**gamma_bnu) * ((z+1.)**gamma_bz)

    return b_nu_z


def plot_bias():

    # reproduce fig 11

    waves = [1500*u.AA, 3000*u.AA, 6000*u.AA]

    for i in waves:
        nu_i = nu_from_lambda(i)
        b = bias(nu_i,z_gals_interp,False)
        plt.plot(z_gals_interp, b, label=r'$\lambda = %g\,{\rm A}$'%i.value)
    
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b(\nu,z)$',fontsize=fontsize)
    plt.legend(loc=2)

    plt.xlim(zmin_gal,zmax_gal)
    plt.ylim(0,3.5)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/bias_nuz.png',bbox_inches='tight')

    plt.show()

    return 

################################


def gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def eps1500(z, vals = False):

    if vals is False:
        eps_1500_z0, gamma_e1500 = fiducials['eps1500']
        eps_1500 = eps_1500_z0 * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * (1.+z)**gamma_e1500 
    elif type(vals) == float or type(vals) == int:
        eps_1500 = vals * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3)
    else:
        eps_1500_z0, gamma_e1500 = vals
        eps_1500 = eps_1500_z0 * (u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * (1.+z)**gamma_e1500 

    return eps_1500

def alpha1500(z, vals = False):

    if vals is False:
        alpha_1500_z0, Calpha1500 = fiducials['alpha1500']
        alpha_1500 = alpha_1500_z0 + Calpha1500*np.log10(1.+z)
    elif type(vals) == float or type(vals) == int:
        alpha_1500 = vals    
    else:
        alpha_1500_z0, Calpha1500 = vals
        alpha_1500 = alpha_1500_z0 + Calpha1500*np.log10(1.+z)

    return alpha_1500

def non_ionizing_continuum(nu,z, vals_eps1500 = False, vals_alpha1500 = False):

    # above the Lya line - 1500 A 

    eps_1500 = eps1500(z, vals_eps1500)

    alpha_1500 = alpha1500(z, vals_alpha1500)

    nu_scale = nu_from_lambda(1500*u.AA)

    eps = eps_1500 * (nu / nu_scale)**alpha_1500

    return eps


def alpha1100(z, vals = False):

    if vals is False:
        alpha_1100_z0, Calpha1100 = fiducials['alpha1100']
        alpha_1100 = alpha_1100_z0 + Calpha1100*np.log10(1.+z)
    elif type(vals) == float or type(vals) == int:
        alpha_1100 = vals
    else:
        alpha_1100_z0, Calpha1100 = vals
        alpha_1100 = alpha_1100_z0 + Calpha1100*np.log10(1.+z)

    return alpha_1100 


def EW_val(z, detector, vals = False):

    if vals is False:
        EWA, EWB = fiducials['EW'](detector)
        if detector == 'ULTRASAT':
            zA = 1.
            zB = 2.
        else:
            zA = 0.3
            zB = 1.
        CLya = (EWB - EWA) / np.log10((1+zB)/(1+zA))
        EW = (CLya * np.log10((1.+z)/(1+zA)) + EWA) 

    elif type(vals) == float or type(vals) == int:
        EW = vals*u.AA
    else:
        EWA, EWB = vals
        if detector == 'ULTRASAT':
            zA = 1.
            zB = 2.
        else:
            zA = 0.3
            zB = 1.
        CLya = (EWB - EWA) / np.log10((1+zB)/(1+zA))
        EW = (CLya * np.log10((1.+z)/(1+zA)) + EWA) 

    return EW

def Lya_line(nu,z, detector, vals_alpha1100 = False, val_EW = False, to_plot = False):

    nu_scale = nu_from_lambda(1216*u.AA)

    # non ionizing continuum 

    alpha_1100 = alpha1100(z,vals_alpha1100)
    eps_Nion_cont = (nu / nu_scale)**alpha_1100

    delta_val = 0.005
    EW = EW_val(z, detector, val_EW) 
    # Lya emission line with equivalent width
    eps_line = nu**2 / (cu.c.to(u.AA/u.s)) * (EW / (nu_from_lambda(1216*u.AA)*delta_val) )

    line =  eps_line  if nu_from_lambda(1216*u.AA)*(1-delta_val)  < nu < nu_from_lambda(1216*u.AA)*(1+delta_val) else 0. 

    eps = eps_Nion_cont + line 
    
    return eps

def fLyC(z, vals):

    # Lyman continuum escape fraction
    if vals is False:
        logf_LyC1, logf_LyC2 = fiducials['fescape']
        CLyC = (logf_LyC2 - logf_LyC1) / (np.log10((1+2.)/(1+1.)))
        f_LyC = pow(10,CLyC * np.log10((1.+z)/(1.+1)) + logf_LyC1)
    elif type(vals) == float or type(vals) == int:
        f_LyC = vals
    else:
        logf_LyC1, logf_LyC2 = vals
        CLyC = (logf_LyC2 - logf_LyC1) / (np.log10((1+2.)/(1+1.)))
        f_LyC = pow(10,CLyC * np.log10((1.+z)/(1.+1)) + logf_LyC1)

    return f_LyC

def Lya_break(nu,z,vals_alpha1100 = False, val_flyc = False,val_alpha900 = False):

    # ioinizing continuum at 900 A - not resolved 
    # related with Lya escape fraction 
    nu_scale_1216 = nu_from_lambda(1216*u.AA)
    nu_scale_912 = nu_from_lambda(912*u.AA)

    alpha_1100 = alpha1100(z,vals_alpha1100)
    if val_alpha900 is False:
        alpha_900 = fiducials['alpha900']
    else:
        alpha_900 = val_alpha900

    # Lyman continuum escape fraction
    f_LyC = fLyC(z,val_flyc) 

    eps_Nion_cont = ( nu_scale_912 / nu_scale_1216 )**alpha_1100

    eps_ion_cont = ( nu / nu_scale_912 )**alpha_900

    eps = f_LyC * eps_Nion_cont * eps_ion_cont

    return eps


def signal(wavelenght,z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900, to_plot = False):

    if wavelenght.value <= 912:
        eps = non_ionizing_continuum(nu_from_lambda(1216*u.AA),z,vals_eps1500,vals_alpha1500) * Lya_break(nu_from_lambda(wavelenght),z,vals_alpha1100,val_flyc,val_alpha900)

    elif 912 < wavelenght.value <= 1216:
        eps = non_ionizing_continuum(nu_from_lambda(1216*u.AA),z,vals_eps1500,vals_alpha1500) * Lya_line(nu_from_lambda(wavelenght),z,detector,vals_alpha1100,val_EW,to_plot)

    else:
        eps = non_ionizing_continuum(nu_from_lambda(wavelenght),z,vals_eps1500,vals_alpha1500)

    return eps 

# def signal_with_DM(wavelenght,z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900, to_plot = False):
# 
    # mass_DM_eV = 5 # eV 
    # nu0 = 1.21e14 * mass_DM_eV * u.Hz # ---> 495 eV , in ultrasat band from z=1
    # gamma = 1e-18 * u.s**-1
# 
    # rho_crit = (2.77536627e11*(u.Msun.to(u.kg)*u.kg*u.Mpc**-3) * cu.c**2 ).to(u.erg/u.Mpc**3)     
    # rhoM = rho_crit*(camb_pars.omegam-camb_pars.omeganu)
# 
    # f_DM = 1
    # rho_DM_0 = f_DM * rhoM
    # rho_DM = (rho_DM_0*(1+z)**3) / u.Hz #* np.exp(-gamma*(t-t0)) --> exponent around 1 
# 
    # eps_DM = (gamma * rho_DM).to(u.erg*u.s**-1*u.Hz**-1*u.Mpc**-3) * gaussian(nu_from_lambda(wavelenght), nu_from_lambda(1216*u.AA), nu_from_lambda(wavelenght)**2 / (cu.c.to(u.AA/u.s)) * EW/1e10) #if int(wavelenght.value) == int(lambda_from_nu(nu0).value) else 0.
# 
    # if wavelenght.value <= 912:
        # eps = non_ionizing_continuum(nu_from_lambda(1216*u.AA),z,vals_eps1500,vals_alpha1500) * Lya_break(nu_from_lambda(wavelenght),z,vals_alpha1100,val_flyc,val_alpha900)
# 
    # elif 912 < wavelenght.value <= 1216:
        # eps = non_ionizing_continuum(nu_from_lambda(1216*u.AA),z,vals_eps1500,vals_alpha1500) * Lya_line(nu_from_lambda(wavelenght),z,detector,vals_alpha1100,val_EW,to_plot)
# 
    # else:
        # eps = non_ionizing_continuum(nu_from_lambda(wavelenght),z,vals_eps1500,vals_alpha1500)
# 
    # return eps + eps_DM


def plot_parameters():

    # reproduce fig 7
 
    z = z_gals_interp
    eps_1500 = eps1500(z)
    #plt.subplot(511)
    plt.plot(z, eps_1500, label=r'$\lambda = %g\,{\rm A}$'%1500)
    plt.ylim(1e25,1e27)
    plt.ylabel(r'$\epsilon_{1500}$',fontsize=fontsize)
    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/eps_1500.png',bbox_inches='tight')
    plt.show()

    alpha_1500 = alpha1500(z)
    #plt.subplot(512)
    plt.plot(z, alpha_1500, label=r'$\lambda = %g\,{\rm A}$'%1500)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{1500}$',fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/alpha_1500.png',bbox_inches='tight')

    plt.show()

    alpha_1100 = alpha1100(z)
    #plt.subplot(513)
    plt.plot(z, alpha_1100, label=r'$\lambda = %g\,{\rm A}$'%1500)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{1100}$',fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/alpha_1100.png',bbox_inches='tight')

    plt.show()

    detector = 'GALEX_NUV' # same for galex_fuv
    EW = EW_val(z,detector)
    #plt.subplot(514)
    plt.plot(z, EW, label=r'$\lambda = %g\,{\rm A}$'%1500)
    plt.ylim(-100,300)
    plt.ylabel(r'$\rm EW_{Ly\alpha} [A]$',fontsize=fontsize)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/EW.png',bbox_inches='tight')

    plt.show()

    f_LyC = fLyC(z,False)
    #plt.subplot(515)
    plt.plot(z, f_LyC, label=r'$\lambda = %g\,{\rm A}$'%1500)
    plt.ylim(0.,1)
    plt.ylabel(r'$\rm \rm f_{LyC}$',fontsize=fontsize)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/fLyC.png',bbox_inches='tight')

    plt.show()

    return 



def plot_rest_signal(detector='GALEX_NUV',vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False):

    # reproduce fig 2, top panel and fig 8

    z = 0 
    wave = np.linspace(700,30000)
    s = np.zeros(len(wave))
    for i in range(len(wave)):
        s[i] = signal(wave[i]*u.AA,z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900).value
        #s[i] = signal(wave[i]*u.AA,z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900).value

    plt.loglog(wave, s, label=r'$z = %g$'%z, color='k')
    #plt.yscale('log')
    plt.xlabel(r'$\lambda_{\rm rest}\,[{\rm A}]$',fontsize=fontsize)
    plt.ylabel(r'$\epsilon_\nu\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.xlim(wave[0],wave[-1])
    plt.ylim(1e24,1e29)

    plt.show()

    return 




def Response(wavelenght,detector):

    # response function from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=GALEX/GALEX.NUV
    if detector == 'GALEX_NUV':
        wave_dat, R_dat = np.genfromtxt('dat/GALEX.NUV.dat').T
        area_NUV = np.trapz(R_dat/wave_dat,wave_dat) 
        interp = interp1d(wave_dat,R_dat / area_NUV,fill_value=0.)
        try:
            R = interp(wavelenght)
        except:
            R = 0.
    elif detector == 'GALEX_FUV':
        wave_dat, R_dat = np.genfromtxt('dat/GALEX.FUV.dat').T
        area_FUV = np.trapz(R_dat/wave_dat,wave_dat) 
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

    else:
        R = 1.

    return R



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

        use_lambda_j = lambda_j[j] * u.AA

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

        use_lambda_j = lambda_j[j] * u.AA

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
    wave_L = 911.8*u.AA # Ly limit

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
    wave_L = 911.8*u.AA # Ly limit

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

    waves = np.linspace(700,20000)*u.AA


    z = [0.,0.3,1.,2,5,10]
    lines = ['-','--',':']
    for i in z:

        tau_ls = np.zeros(len(waves))
        tau_lc = np.zeros(len(waves))
        tau_tot = np.zeros(len(waves))
        for l in range(len(waves)):
            tau_ls[l] = tau_LS(waves[l],i)
            tau_lc[l] = tau_LC(waves[l],i)
            tau_tot[l] = tau_Lya(waves[l],i)

        plt.loglog(waves,tau_tot,color=colors[z.index(i)],label=r'$z = %g$'%i)
# 
        #plt.plot(waves,tau_tot,color=colors[z.index(i)],linestyle=lines[0],label=r'$z=%g$'%i)
        #if z.index(i) == 0:
        #  plt.plot(waves,tau_lc,color=colors[z.index(i)],linestyle='-',#lines[1],
        #  label=r'${\rm Ly\alpha\,\, continuum}$')
        #  plt.plot(waves,tau_ls,color=colors[z.index(i)],linestyle=lines[1],label=r'${\rm Ly\alpha\,\, series}$')
        #else:
        #    plt.plot(waves,tau_lc,color=colors[z.index(i)],linestyle='-',label=r'$z=%g$'%i)
        #    plt.plot(waves,tau_ls,color=colors[z.index(i)],linestyle=lines[1])
# 
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm A}]$',fontsize=fontsize)
    plt.ylabel(r'$\tau(\lambda_{\rm obs},z)$',fontsize=fontsize)

    plt.legend(ncol=2,loc=1)
    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/tau.png',bbox_inches='tight')

    plt.show()

    return 


def dJdz(z, detector,run = False,\
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

        unit = (Response(lambda_from_nu(1.*u.Hz), detector) *signal(lambda_from_nu(1*u.Hz),0.,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900) * np.exp(-tau_Lya(lambda_from_nu(1*u.Hz),0.))/(1*u.Hz)).unit

        intg = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        rest_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.AA).value
            intg_arr[i] = intg(nu_obs_arr[i])

        #dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * quad(intg,nu_min.value,nu_max.value)[0]*(unit*u.Hz/u.steradian) 
        dJdz = cu.c.to(u.km/u.s) / (4*np.pi*H(z)*(1+z)) * np.trapz(intg_arr,nu_obs_arr)*(unit*u.Hz/u.steradian) 

    else:
        zval, dJdzval = np.genfromtxt(filename)
        dJdz = interp1d(zval,dJdzval)(z) * u.Jy/u.steradian  

    return (dJdz.to(u.Jy/u.steradian)).value



def bJ_z(z, detector, run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = ''):

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
        intg_num = lambda nu_obs: bias(nu_obs*u.Hz*(1+z),z,val_bias) * Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 
#
        intg_den = lambda nu_obs: Response(lambda_from_nu(nu_obs*u.Hz), detector) * signal(lambda_from_nu(nu_obs*u.Hz*(1+z)),z,detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs*u.Hz),z))/nu_obs 

        rest_wave = np.linspace(lambda_from_nu(nu_min).value,lambda_from_nu(nu_max).value,500)

        nu_obs_arr = np.zeros(len(rest_wave))
        intg_arr = np.zeros(len(rest_wave))
        intg_den_arr = np.zeros(len(rest_wave))
        for i in range(len(rest_wave)):
            nu_obs_arr[i] = nu_from_lambda(rest_wave[i]*u.AA).value
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


def compute_vals(reduced_z = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False):

    reduced_label = '_reduced' if reduced_z else ''

    with Pool(6) as pool:
        
        print('Doing NUV')
        dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900, filename='./results/EBL/dJdz_GALEX_NUV' + reduced_label + '.dat')
        if reduced_z:
            dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
            np.savetxt('./results/EBL/dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_nuv)))
        else:
            dJ_nuv = pool.map(dJdz_NUV, z_gals('SDSS'))
            np.savetxt('./results/EBL/dJdz_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_nuv)))

        print('Doing FUV')
        dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
            np.savetxt('./results/EBL/dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_fuv)))
        else:
            dJ_fuv = pool.map(dJdz_FUV, z_gals('SDSS'))
            np.savetxt('./results/EBL/dJdz_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(dJ_fuv)))
    
        print('Doing b NUV')
        bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/bJ_GALEX_NUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
            np.savetxt('./results/EBL/bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_nuv)))
        else:
            bJ_nuv = pool.map(bJ_nuv_f, z_gals('SDSS'))
            np.savetxt('./results/EBL/bJ_GALEX_NUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_nuv)))

        print('Doing b FUV')
        bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/bJ_GALEX_FUV' + reduced_label + '.dat')

        if reduced_z:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)
            np.savetxt('./results/EBL/bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_fuv)))
        else:
            bJ_fuv = pool.map(bJ_fuv_f, z_gals('SDSS'))
            np.savetxt('./results/EBL/bJ_GALEX_FUV' + reduced_label + '.dat', (z_gals('SDSS'),np.asarray(bJ_fuv)))
     
        print('Doing ULTRASAT')
        dJdz_f = partial(dJdz,detector='ULTRASAT',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_ULTRASAT' + reduced_label + '.dat')

        if reduced_z:
            dJ_U = pool.map(dJdz_f, z_gals_interp)
            np.savetxt('./results/EBL/dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(dJ_U)))
        else:
            dJ_U = pool.map(dJdz_f, z_gals)
            np.savetxt('./results/EBL/dJdz_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(dJ_U)))

        print('Doing b ULTRASAT')
        bJ_f =  partial(bJ_z,detector='ULTRASAT',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/bJ_ULTRASAT' + reduced_label + '.dat')
        if reduced_z:
            bJ_U = pool.map(bJ_f, z_gals_interp)        
            np.savetxt('./results/EBL/bJ_ULTRASAT' + reduced_label + '.dat', (z_gals_interp,np.asarray(bJ_U)))
        else:
            bJ_U = pool.map(bJ_f, z_gals)        
            np.savetxt('./results/EBL/bJ_ULTRASAT' + reduced_label + '.dat', (z_gals,np.asarray(bJ_U)))
    return




def plot_bdJdz_multi():

    # reproduce fig 3
    #alpha_1500 = [-1.,-0.5,0,0.5,1.,False]
    alpha_1500 = [-1.,0,1.,False]

    #EW = [-5,-2.5,0,2.5,5,False]
    EW = [-50,0,50,False]
    
    #f_Lyc = [0.,0.03,0.1,0.3,1,False]
    f_Lyc = [0.,0.1,1,False]
    #color= ['b','c','k','orange','r','k']
    color= ['b','k','r','k']

    vals_eps1500= False
    vals_alpha1100= False
    val_alpha900= False
    
    # plt.figure()
    # val_EW= False
    # val_flyc= False
 
    # print('Doing alpha1500')
    # for A in tqdm(alpha_1500):
 
    #     vals_alpha1500= A
 
    #     with Pool(6) as pool:

    #         print('Doing NUV')
    #         dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
    #         dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)

    #         print('Doing FUV')
    #         dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
    #         dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
    
    #         print('Doing bias NUV')
    #         bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_NUV.dat')
    #         bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
 
    #         print('Doing bias FUV')
    #         bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_FUV.dat')
    #         bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)

    #         bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
    #         bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

    #     if not vals_alpha1500:
    #         linestyle = '--'
    #     else: 
    #         linestyle = '-'

    #     plt.subplot(121)
    #     # plt.subplot(231)
    #     plt.plot(z_gals('SDSS'),bdJ_fuv,color= color[alpha_1500.index(A)],label=r'$\alpha_{1500} = %g$'%A,linestyle=linestyle)
 
    #     # plt.subplot(234)
    #     plt.subplot(122)
    #     plt.plot(z_gals('SDSS'),bdJ_nuv,color= color[alpha_1500.index(A)],label=r'$\alpha_{1500} = %g$'%A,linestyle=linestyle)
 
    # # plt.subplot(231)
    # plt.subplot(121)
    # plt.xlim(z_gals_interp,z_gals('SDSS')[-1])
    # #plt.ylim(-10,65)
    # plt.xlabel(r'$z$',fontsize=fontsize)
    # plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    # plt.legend(loc=1)
 
    # plt.subplot(122)
    # # plt.subplot(234)
    # plt.xlim(z_gals_interp,z_gals('SDSS')[-1])
    # #plt.ylim(-10,65)
    # plt.xlabel(r'$z$',fontsize=fontsize)
    # plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    # plt.legend(loc=1)


    plt.figure()

    vals_alpha1500= False
    val_flyc= False

    print('Doing EW')
    for A in tqdm(EW):
        val_EW = A

        with Pool(6) as pool:

            print('Doing NUV')
            dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
            dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
            print('Doing FUV')
            dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
            dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
    
            bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_NUV.dat')
            bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
 
            print('Doing b FUV')
            bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_FUV.dat')
            bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)

            bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
            bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

        if not val_EW:
            linestyle = '--'
        else: 
            linestyle = '-'
        # plt.subplot(232)
        plt.subplot(121)
        plt.plot(z_gals_interp,bdJ_fuv,color= color[EW.index(A)],label=r'$\rm EW = %g$'%A,linestyle=linestyle)

        plt.subplot(122)
        # plt.subplot(235)
        plt.plot(z_gals_interp,bdJ_nuv,color= color[EW.index(A)],label=r'$\rm EW = %g$'%A,linestyle=linestyle)

    plt.subplot(121)
    # plt.subplot(232)
    plt.xlim(z_gals('SDSS')[0],z_gals('SDSS')[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.subplot(122)
    # plt.subplot(235)
    plt.xlim(z_gals('SDSS')[0],z_gals('SDSS')[-1])
    #plt.ylim(-10,65)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)


#     plt.figure()
#     val_EW= False
#     vals_alpha1500= False

#     print('Doing fLyC')
#     for A in tqdm(f_Lyc):

#         val_flyc= A
#         with Pool(6) as pool:

#             print('Doing NUV')
#             dJdz_NUV = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
#             dJ_nuv = pool.map(dJdz_NUV, z_gals_interp)
#             print('Doing FUV')
#             dJdz_FUV = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,)
#             dJ_fuv = pool.map(dJdz_FUV, z_gals_interp)
    
#             bJ_nuv_f =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_NUV.dat')
#             bJ_nuv = pool.map(bJ_nuv_f, z_gals_interp)
 
#             print('Doing b FUV')
#             bJ_fuv_f =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = False,filename='results/EBL/dJdz_GALEX_FUV.dat')
#             bJ_fuv = pool.map(bJ_fuv_f, z_gals_interp)

#             bdJ_nuv = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
#             bdJ_fuv = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 

#         if not val_flyc:
#             linestyle = '--'
#         else: 
#             linestyle = '-'

#         plt.subplot(121)
#         # plt.subplot(233)
#         plt.plot(z_gals_interp,bdJ_fuv,color= color[f_Lyc.index(A)],label=r'$f_{\rm LyC} = %g$'%A,linestyle=linestyle)

#         plt.subplot(122)
#         # plt.subplot(236)
#         plt.plot(z_gals_interp,bdJ_nuv,color= color[f_Lyc.index(A)],label=r'$f_{\rm LyC} = %g$'%A,linestyle=linestyle)

#     plt.subplot(121)
#     # plt.subplot(233)
#     plt.xlim(z_gals('SDSS')[0],z_gals('SDSS')[-1])
# #    plt.ylim(-10,65)
#     plt.xlabel(r'$z$',fontsize=fontsize)
#     plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
#     plt.legend(loc=1)

#     # plt.subplot(236)
#     plt.subplot(122)
#     plt.xlim(z_gals('SDSS')[0],z_gals('SDSS')[-1])
# #    plt.ylim(-10,65)
#     plt.xlabel(r'$z$',fontsize=fontsize)
#     plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
#     plt.legend(loc=1)

#     plt.show()
    return 


def plot_bdJdz_epsbias():

    color= ['b','k','r','k']

    vals_eps1500_mod = [fiducials['eps1500'][0]*2,fiducials['eps1500'][1]]
    vals_bias_mod = [fiducials['bias'][0]*2,fiducials['bias'][1],fiducials['bias'][2]]

    plt.figure()
 
    with Pool(6) as pool:

        print('Doing NUV')
        #dJdz_NUV_fid = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        #dJ_nuv = pool.map(dJdz_NUV_fid, z_gals_interp)
#
        #dJdz_NUV_2eps = partial(dJdz,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        #dJ_nuv_2eps = pool.map(dJdz_NUV_2eps, z_gals_interp)

        print('Doing FUV')
        #dJdz_FUV_fid = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        #dJ_fuv = pool.map(dJdz_FUV_fid, z_gals_interp)
#
        #dJdz_FUV_2eps = partial(dJdz,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        #dJ_fuv_2eps = pool.map(dJdz_FUV_2eps, z_gals_interp)

        print('Doing ULTRASAT')
        dJdz_ULTRASAT_fid = partial(dJdz,detector='ULTRASAT',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        dJ_ULTRASAT = pool.map(dJdz_ULTRASAT_fid, z_gals_interp)

        dJdz_ULTRASAT_2eps = partial(dJdz,detector='ULTRASAT',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,)
        dJ_ULTRASAT_2eps = pool.map(dJdz_ULTRASAT_2eps,z_gals_interp)


        print('Doing bias NUV')
        #bJ_nuv_f_fid =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/dJdz_GALEX_NUV.dat')
        #bJ_nuv = pool.map(bJ_nuv_f_fid,z_gals_interp)
#
        #bJ_nuv_f_2eps =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/dJdz_GALEX_NUV.dat')
        #bJ_nuv_2eps = pool.map(bJ_nuv_f_2eps,z_gals_interp)

        #bJ_nuv_f_2bias =  partial(bJ_z,detector='GALEX_NUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = vals_bias_mod,filename='results/EBL/dJdz_GALEX_NUV.dat')
        #bJ_nuv_2bias = pool.map(bJ_nuv_f_2bias, z_gals_interp)
#
        print('Doing bias FUV')
        #bJ_fuv_f_fid =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/dJdz_GALEX_FUV.dat')
        #bJ_fuv = pool.map(bJ_fuv_f_fid, z_gals_interp)
#
        #bJ_fuv_f_2eps =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='...')
        #bJ_fuv_2eps = pool.map(bJ_fuv_f_2eps, z_gals_interp)
#
        #bJ_fuv_f_2bias =  partial(bJ_z,detector='GALEX_FUV',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = vals_bias_mod,filename='...')
        #bJ_fuv_2bias = pool.map(bJ_fuv_f_2bias, z_gals_interp)

        print('Doing bias ULTRASAT')
        bJ_ULTRASAT_f_fid =  partial(bJ_z,detector='ULTRASAT',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='...')
        bJ_ULTRASAT = pool.map(bJ_ULTRASAT_f_fid, z_gals_interp)

        bJ_ULTRASAT_f_2eps =  partial(bJ_z,detector='ULTRASAT',run=True,vals_eps1500=vals_eps1500_mod,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='...')
        bJ_ULTRASAT_2eps = pool.map(bJ_ULTRASAT_f_2eps, z_gals_interp)

        bJ_ULTRASAT_f_2bias =  partial(bJ_z,detector='ULTRASAT',run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = vals_bias_mod,filename='...')
        bJ_ULTRASAT_2bias = pool.map(bJ_ULTRASAT_f_2bias,z_gals_interp)



        #bdJ_nuv_fid = np.asarray(dJ_nuv) * np.asarray(bJ_nuv)
        #bdJ_fuv_fid = np.asarray(dJ_fuv) * np.asarray(bJ_fuv) 
        bdJ_ULTRASAT_fid = np.asarray(dJ_ULTRASAT) * np.asarray(bJ_ULTRASAT) 

        #bdJ_nuv_2eps = np.asarray(dJ_nuv_2eps) * np.asarray(bJ_nuv_2eps)
        #bdJ_fuv_2eps = np.asarray(dJ_fuv_2eps) * np.asarray(bJ_fuv_2eps) 
        bdJ_ULTRASAT_2eps = np.asarray(dJ_ULTRASAT_2eps) * np.asarray(bJ_ULTRASAT_2eps) 

        #bdJ_nuv_2bias = np.asarray(dJ_nuv) * np.asarray(bJ_nuv_2bias)
        #bdJ_fuv_2bias = np.asarray(dJ_fuv) * np.asarray(bJ_fuv_2bias) 
        bdJ_ULTRASAT_2bias = np.asarray(dJ_ULTRASAT) * np.asarray(bJ_ULTRASAT_2bias) 

    print(bdJ_ULTRASAT_2eps)
    print(bdJ_ULTRASAT_2bias)

    plt.subplot(121)
    #plt.plot(z_gals_interp,bdJ_fuv_fid,color= color_FUV,label=r'$\rm Fiducial\, GALEX\, FUV$')
    #plt.plot(z_gals_interp,bdJ_nuv_fid,color= color_NUV,label=r'$\rm Fiducial\, GALEX\, NUV$')
 
    #plt.plot(z_gals_interp,bdJ_fuv_2eps,color= color_FUV,linestyle='--',label=r'$\rm 2\epsilon_{1500}^{z=0}$')
    #plt.plot(z_gals_interp,bdJ_nuv_2eps,color= color_NUV,linestyle='--',label=r'$\rm 2\epsilon_{1500}^{z=0}$')
 
    #plt.plot(z_gals_interp,bdJ_fuv_2bias,color= color_FUV,linestyle=':',label=r'$\rm 2b_{1500}^{z=0}$')
    #plt.plot(z_gals_interp,bdJ_nuv_2bias,color= color_NUV,linestyle=':',label=r'$\rm 2b_{1500}^{z=0}$')
 
    plt.xlim(z_gals_interp[0],z_gals_interp[-1])
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)
 
    plt.subplot(122)
    plt.plot(z_gals_interp,bdJ_ULTRASAT_fid,color= color_ULTRASAT,label=r'$\rm Fiducial\, ULTRASAT$') 
    plt.plot(z_gals_interp,bdJ_ULTRASAT_2eps,color= color_ULTRASAT,linestyle='--',label=r'$\rm 2\epsilon_{1500}^{z=0}$') 
    plt.plot(z_gals_interp,bdJ_ULTRASAT_2bias,color= 'k',linestyle=':',label=r'$\rm 2b_{1500}^{z=0}$')
 
    plt.xlim(z_gals_interp[0],z_gals_interp[-1])
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\rm b\,dJ/dz\,[Jy/sr]$',fontsize=fontsize)
    plt.legend(loc=1)
 


    return 


# def plot_signal_withDM(vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False):

#     dJ_U =  dJdz(z_gals,detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_ULTRASAT_TESTGAUS.dat')

    
#     dJ_U_DM =  dJdz(z_gals,detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_ULTRASAT_withDM.dat')


#     bJ_U = bJ_z(z_gals,detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = val_bias,filename='results/EBL/bJ_ULTRASAT_TESTGAUS.dat')

#     bJ_U_DM = bJ_z(z_gals,detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha11000=vals_alpha11000,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = val_bias,filename='results/EBL/bJ_ULTRASAT_withDM.dat')

#     # figure 10 
#     plt.figure()
#     plt.plot(z_gals,dJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
#     #plt.plot(z_gals,dJ_U_DM,label=r'$\rm with\, DM$',color='k')
#     #plt.ylim(-10,200)
#     plt.xlabel(r'$z$',fontsize=fontsize)
#     plt.ylabel(r'$dJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
#     plt.legend(loc=4, ncol=2)
    
#     plt.tight_layout()
#     plt.savefig('results/PLOTS/EBL/dJdz_withDM.png',bbox_inches='tight')

#     plt.figure()
#     plt.plot(z_gals,bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
#     plt.plot(z_gals,bJ_U_DM,label=r'$\rm with\,DM$',color='k')
#     plt.xlabel(r'$z$',fontsize=fontsize)
#     plt.ylabel(r'$b_J$',fontsize=fontsize)
#     plt.xlim(0,3)
#     plt.legend(loc=4, ncol=2)
    
#     plt.tight_layout()
#     plt.savefig('results/PLOTS/EBL/bJz_withDM.png',bbox_inches='tight')

#     plt.figure()
#     plt.plot(z_gals,dJ_U*bJ_U,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
#     plt.plot(z_gals,dJ_U_DM*bJ_U_DM,label=r'$\rm with\, DM$',color='k')
#     #plt.ylim(-10,200)
#     plt.xlabel(r'$z$',fontsize=fontsize)
#     plt.ylabel(r'$b_JdJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
#     plt.legend(loc=4, ncol=2)
    
#     plt.tight_layout()
#     plt.savefig('results/PLOTS/EBL/bJdJdz_withDM.png',bbox_inches='tight')


#     plt.show()
#     return 
