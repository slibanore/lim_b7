# SL: last update 01/18/2023

from LIM_b7.modelling.overall_setting import * 


###########################################
# SPECIFIC ANALYSIS RELATED PARAMETERS
###########################################

# line rest frame
rest_nu = 1.17e6 

# speed of light
light = 2.998e8 # m / s 

# kmax for Pk computation
kmax_fid = 100.*u.Mpc**-1

#developer = 'CDS' 
CDS_alpha_1_fid = 1.
CDS_alpha_2_fid = 1.
CDS_alpha_3_fid = 1. 
CDS_alpha_4_fid = 1.
k_12_fid = 0.5
k_23_fid = 2.
k_34_fid = 10.

#developer = 'axions' 
m_axion_fid = 1e-24
f_axion_fid = 0.1

# developer = PMF
n_B_fid = -2.9
sigma_B_0_fid = 0.1
smooth_scale_fid = 1

###########################################
# STATISTICAL ANALYSIS PARAMETERS
###########################################

delta = 1.e-3 # step for numerical derivatives
error_type = 'poisson' # VID stat error
Npt_i = 200 # Ti bins
min_Ti = 3. # min Ti bin
max_Ti = 100000. # max Ti bin
min_Bi = 1e-1 # min value of the VID

do_Jysr_fid = True 

###########################################
# DEFAULT HALO MODEL
###########################################

hmf_fid = 'Tinker'
delta_tinker = 200.
alpha_tinker_fid=10**(-(0.75/np.log10(delta_tinker/75.0))**1.2)
A_tinker_fid= lambda z: 0.186*(1.0+z)**(-0.14)
a_tinker_fid = lambda z: 1.47*(1.0+z)**(-0.06)
b_tinker_fid = lambda z: 2.57*(1.0+z)**(-alpha_tinker_fid)
c_tinker_fid = lambda z: 1.19

fNL_fid = 0.0

hmf_pars_fid = lambda z: dict(A_tinker = A_tinker_fid(z),a_tinker = a_tinker_fid(z),b_tinker = b_tinker_fid(z),c_tinker = c_tinker_fid(z))

###########################################
# DEFAULT ASTRO MODEL
###########################################

astro_fid = 'Lyalpha'
# astro_alpha_fid = 1.37
# astro_beta_fid = -1.74
# astro_dMF_fid = 1.0
# astro_sig_SFR_fid = 0.3

astro_sig_scatter_fid = 0.3
astro_Lcut_fid = 2.6e7*u.Lsun

SFR_file_fid = os.getcwd() + '/LIM_b7/modelling/SFR_tables/sfr_release.dat'

csi_Lyalpha_fid = 1.6
z0_Lyalpha_fid = 0.87
zeta_Lyalpha_fid = 1/4.
f0_Lyalpha_fid = 0.35
SFR0_Lyalpha_fid = 1.29 # Msun yr-1
t_Lyalpha_fid = 1.25

astro_par_fid = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)


Mmin_fid = 1e10*u.Msun
Mmax_fid = 1e13*u.Msun

Lmin_fid = 1e-5*u.Lsun
Lmax_fid = 1e10*u.Lsun

###########################################
# COSMO_ASTRO PARAMETERS
###########################################

# The full set of parameters is set in
# cosmo_astro prescriptions, here you only
# have the ones we need to change for 

astrocosmo_dict = lambda developer, z: \
	dict(\

	cosmo_input_camb = dict(
	f_NL=0, H0=67.67, cosmomc_theta=None,
	ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
	num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
	YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
	TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
	theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
	dark_energy_model='ppf',As=2.105e-09, ns=0.967, nrun=0., nrunrun=0.0, r=0.0, nt=None, ntrun=0.0, 
	pivot_scalar=0.05, pivot_tensor=0.05,
	parameterization=2,halofit_version='mead'),

	# HALO PROPERTIES
	Mmin = Mmin_fid,
	Mmax = Mmax_fid,
	hmf_model = hmf_fid,
	hmf_pars = hmf_pars_fid(z),

	# ASTROPHYSICAL MODEL FOR LINE EMISSION
	model_name = astro_fid,
	model_par = astro_par_fid,

	sigma_scatter = astro_sig_scatter_fid,
	Lmin = Lmin_fid,
	Lmax = Lmax_fid,

	# POWER SPECTRUM - RELATED QUANTITIES
	#kmin = 1e-2*u.Mpc**-1,
	kmax = kmax_fid,

	# VID - RELATED QUANTITES
	#Tmin_VID = 1e-5*u.uK,
	#Tmax_VID = 80.*u.uK,
	do_fast_VID = True,
	do_Jysr = do_Jysr_fid,
	#Ngal_max = 500,

	# model to be analysed and model related parameters
	developer = developer,

	#developer = 'CDS' 
	CDS_alpha_1 = CDS_alpha_1_fid, 
	CDS_alpha_2 = CDS_alpha_2_fid, 
	CDS_alpha_3 = CDS_alpha_3_fid, 
	CDS_alpha_4 = CDS_alpha_4_fid,
	k_12 = k_12_fid,
	k_23 = k_23_fid,
	k_34 = k_34_fid,

	#developer = 'axions' 
	m_axion = m_axion_fid,
	f_axion = f_axion_fid,

	#developer = 'PMF'
	n_B = n_B_fid,         
	sigma_B_0 = sigma_B_0_fid,    
	smooth_scale = smooth_scale_fid,   
	)


###########################################
# SURVEY PARAMETERS
###########################################

# possible default experiments 
# see arXiv:2208.01658

# High z experiment (z = 5.3, 6.25, 7.25, 8)
# Possible configurations:
# COMAP-EoR , COMAP-EoR Deep, COS3 
obs_params_highz = lambda z, which: \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 1*u.GHz,
			dnu = 2*u.MHz,
			Tsys_NEFD = 20*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 38,
			beam_FWHM = (4.2)*u.arcmin,
			tobs =  7000.*u.hr,
			Omega_field =  0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4*u.deg**2,
			Nfield = 1) if z == 8.2 else \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 2*u.GHz,
			dnu = 2*u.MHz,
			Tsys_NEFD = 20*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 38,
			beam_FWHM = (4)*u.arcmin,
			tobs =  7000.*u.hr,
			Omega_field =  0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4.*u.deg**2,
			Nfield = 1) if z == 7.25 else \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 2*u.GHz,
			dnu = 2*u.MHz,
			Tsys_NEFD = 22*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 38,
			beam_FWHM = (3.7)*u.arcmin,
			tobs = 7000.*u.hr,
			Omega_field =  0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4.*u.deg**2,
			Nfield = 1) if z == 6.25 else \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 3*u.GHz,
			dnu = 2*u.MHz,
			Tsys_NEFD = 27*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 38,
			beam_FWHM = (3.3)*u.arcmin,
			tobs =  7000.*u.hr,
			Omega_field =  0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4.*u.deg**2,
			Nfield = 1) if z == 5.3 else -1


# Low z experiment (z = 2.60, 3)
# Possible configurations:
# COMAP-EoR , COMAP-EoR Deep, COS3 
obs_params_lowz = lambda z, which: \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 4*u.GHz,
			dnu = 2*u.MHz, 
			Tsys_NEFD = 34*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 19,
			beam_FWHM = (4.5)*u.arcmin,
			tobs = 5000.*u.hr,
			Omega_field = 0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4.*u.deg**2,
			Nfield = 1) if z == 2.60 else \
		dict(nuObs = (rest_nu/(1.+z))*u.GHz,
			Delta_nu = 4*u.GHz,
			dnu = 2*u.MHz, 
			Tsys_NEFD = 34*u.K,
			Nfeeds = 1000 if which == 'COS3' else 1000 if which == 'wide' else 19,
			beam_FWHM = (3.9)*u.arcmin,
			tobs = 5000.*u.hr,
			Omega_field =  0.5*u.deg**2 if which == 'deep' else 10*u.deg**2 if which == 'wide' else 4.*u.deg**2,
			Nfield = 1) if z == 3. else -1


# COMAP 1 (see arXiv:1907.10065)
obs_params_COMAP1 = lambda z: dict(\
        nuObs = (rest_nu / (1.+ z))*u.GHz, 
        Delta_nu = 8*u.GHz, 
        dnu = 15.6*u.MHz,
        Tsys_NEFD = 44*u.K, 
        Nfeeds = 19,
        beam_FWHM = 4*u.arcmin, 
        tobs = 5000.*u.hr, 
        Omega_field = 4*u.deg**2, 
        Nfield = 1)  

# IMS3 (1907.10065)
obs_params_IMS3 = lambda z: dict(\
        nuObs = (rest_nu / (1.+ z))*u.GHz,
			Delta_nu = 24*u.GHz,
			dnu = (rest_nu / (1.+ z))*1e3/9150*u.MHz,
			Tsys_NEFD = max(20,(rest_nu / (1.+ z)))*u.K,
			Nfeeds = 1000,
			beam_FWHM = 4*u.arcmin,
			tobs = 1e4*u.hr,
			Omega_field = 1000.*u.deg**2,
			Nfield = 1)

# ULTRASAT (Libanore, Kovetz (2023) - ULTRASAT white paper)
obs_params_ULTRASAT = lambda z: dict(\
        nuObs = (rest_nu / (1.+ z))*u.GHz, 
        Delta_nu = 2.70e8*u.GHz, 
        dnu = (2.70e8*u.GHz).to(u.MHz),
		mag_limit = 23.,
        nsigma = 5,
        px_side = 5.45*u.arcsec,
        combined_pixels = 3,
        Nfeeds = 1,
        beam_FWHM = 8.3*u.arcsec, 
        tobs = 540*u.hr, 
        Omega_field = (4*np.pi*(180/np.pi)**2)*u.deg**2, # full sky
        Nfield = 1)  

