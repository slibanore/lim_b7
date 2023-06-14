# SL: last update 01/18/2023

from .VID_computation import *


default_par = dict(

    # COSMOLOGICAL MODEL
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
    Mmin = 1e9*u.Msun,
    Mmax = 1e15*u.Msun,
    nM = 1000,
    hmf_model = 'Tinker',
    tinker_pars = dict(\
    A_tinker = 0.186, a_tinker = 1.47,
    b_tinker = 2.57, c_tinker = 1.19),
    bias_model = 'Tinker10',
    bias_par = {},

    # ASTROPHYSICAL MODEL FOR LINE EMISSION
    nu = 115.271*u.GHz,
    nuObs = 30*u.GHz,
    model_name = 'TonyLi',
    model_par = {\
    'alpha':1.37, 'beta':-1.74, \
    'dMF':1.,'sig_SFR':0.3, \
    'SFR_file':os.getcwd() + '/LIM_b7/modelling/SFR_tables/sfr_release.dat'},
    sigma_scatter=0.3,
    Lmin = 1e-5*u.Lsun,
    Lmax = 1e10*u.Lsun,
    nL = 5000,
    dndL_Lcut=50.*u.Lsun,

    # POWER SPECTRUM - RELATED QUANTITIES
    kmin = 1e-2*u.Mpc**-1,
    kmax = 100.*u.Mpc**-1,
    nk = 500,
    fduty=1.,
    nonlinear=False,
    do_onehalo=False,
    do_Jysr=False,
    do_RSD=False,
    sigma_NL=7*u.Mpc,
    nmu=1000,
    FoG_damp='Lorentzian',
    smooth=False,
    do_conv_Wkmin = False,

    # VID - RELATED QUANTITES
    Tmin_VID=1e-3*u.uK,
    Tmax_VID=80.*u.uK,
    nT=2**22,
    do_fast_VID=True,
    Ngal_max=500,
    subtract_VID_mean=True,
    linear_VID_bin=False,
    do_sigma_G = True,
    sigma_G_input = 1.6,

    # set which kind of analysis you want to use
    developer = 'CDS', 
    
    ###########################################
    # PROPERTIES RELATED WITH SPECIFIC ANALYSIS
    ###########################################

    # CDS : arXiv:2208.01658
    # LCDM deviations and related scales
    CDS_alpha_1=1., 
    CDS_alpha_2=1., 
    CDS_alpha_3=1., 
    CDS_alpha_4=1., 
    k_12 = 0.5, 
    k_23 = 2., 
    k_34 = 10.,

    # axions
    # mass and damping scale for fuzzy DM 
    m_axion = 1e-24, 
    f_axion = 0.1, 

    # PMF
    n_B = -2.9,         # PMF spectral index
    sigma_B_0 = 0.1,    # PMF variance in nGauss
    smooth_scale = 1,   # PMF smooth scale in Mpc/h

    )


# create lim objects
def lim_VID(params=default_par):
    
    if type(params) is not dict:
        raise ValueError('model_params must be a dict')
        
    #Check whether there is any parameter wrong
    defpars = get_default_params(LineModel.__init__)
    x1 = get_default_params(LineObs.__init__)
    defpars.update(x1)

    x2 = get_default_params(VID.__init__)
    defpars.update(x2)

    for key in params:
        if key not in defpars:
            raise ValueError(key+" is not a valid input parameter for lim_VID. Check your input, please.")
    
    return VID(**params)


###########################################
# OTHER FUNCTIONS
###########################################

def cut_T(model, Ti, Bi, max_Ti, min_Ti):

    if type(max_Ti) == list:
        mask_voxels = max(int(model.Nvox * max_Ti[0]),1)
        print('Masked voxels = ' + str(mask_voxels))
        remove_vox_id = len(Bi) - 1 
        removed_vox = 0 
        while removed_vox < mask_voxels:
            for i in range(int(Bi[remove_vox_id])):
                removed_vox += 1
                Bi[remove_vox_id] -= 1
                if removed_vox == mask_voxels:
                    break
            remove_vox_id -= 1
        use_max_Ti = Ti[remove_vox_id+1] 
    else:			
        use_max_Ti = max_Ti
        
    got_Ti = False
    min_Ti_id = 0 
    got_Ti_max = False
    use_max_Ti_id = -1 
    for i in range(len(Ti)):
        try:
            if Ti[i].value >= min_Ti:
                min_Ti_id = i
                got_Ti = True
            if got_Ti:
                break
        except:
            if Ti[i] >= min_Ti:
                min_Ti_id = i
                got_Ti = True
            if got_Ti:
                break

    for i in range(len(Ti)):
        try:
            if Ti[i].value > use_max_Ti:
                use_max_Ti_id = i
                got_Ti_max = True
        except:
            if Ti[i] > use_max_Ti:
                use_max_Ti_id = i
                got_Ti_max = True
        if got_Ti_max:
            break        

    Ti = np.asarray(Ti[min_Ti_id:use_max_Ti_id])*u.uK
    Bi = np.asarray(Bi[min_Ti_id:use_max_Ti_id])*u.dimensionless_unscaled

    return Ti, Bi


def save_VID(filename, Ti, Bi):

    f = open(filename,'w')
    f.write('Ti\tBi\n')
    for i in range(len(Ti)):
        f.write(str(Ti[i]) + '\t')
        f.write(str(Bi[i]) + '\n')
    f.close()

    return 


def import_VID(filename, label, model, max_Ti, min_Ti, cut_TB):

    print(label)
    f = open(filename,'r')
    lines = f.readlines()

    Ti_full = []
    Bi_full = []
    for x in range(len(lines)-1):
        words = lines[x+1].split('\t')
        Ti_full.append(float(words[0][:-3]))
        Bi_full.append(float(words[1]))
        
    if cut_TB:
        Ti, Bi = cut_T(model, Ti_full, Bi_full, max_Ti, min_Ti)
    else:
        Ti = np.asarray(Ti_full)*u.uK
        Bi = np.asarray(Bi_full)*u.dimensionless_unscaled

    return Ti, Bi 


def save_fisher(filename, derive_pars,Fisher):
    f = open(filename,'w')
    f.write(str(derive_pars))
    f.write('\n')
        
    f.write('\n\n')
    for i in range(len(Fisher)):
        f.write('[')
        for j in range(len(Fisher[i])):
            f.write(str(Fisher[i][j]))
            if j != len(Fisher[i])-1:
                f.write(',')
            else:
                f.write(']')
        f.write('\n')
    f.close()

    return 