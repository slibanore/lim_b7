# SL: last update 01/28/2023

import get_Pk as pk 
import Pkline_fisher as pf
from LIM_b7.fiducial_pars import *



###########################################
# CHOOSE MODULE
###########################################

get_matter_pk = True
get_line_pk = True

do_RSD = True
smooth = True
do_conv_Wkmin = False
nonlinear = False

save_flag = True
plot_flag = True

Pk_fisher = False

###########################################
# SETUP ANALYSIS
###########################################

model_data = dict(\
    developer = 'CDS',               

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
    hmf_model = 'Tinker',
    bias_model = 'Tinker10',

    # POWER SPECTRUM - RELATED QUANTITIES
    kmin = kmin_fid,
    kmax = kmax_fid,
    do_onehalo=True,
    do_RSD=do_RSD,
    smooth= smooth,
    do_conv_Wkmin = do_conv_Wkmin,
    nonlinear=nonlinear,
    nk = 500,
    fduty=1.,
    nmu=1000,

    )

survey = 'COS3' # EOR, deep, wide
redshift = 6.25 # 2.6, 3, 5.3, 6.25, 7.25

def main():

    ###########################################
    # COMPUTE AND PLOT POWER SPECTRUM
    ###########################################

    pk.run_pk(get_matter_pk,
           get_line_pk,
           model_data,
           survey,
           redshift,
           save_flag,
           plot_flag)

    ###########################################
    # COMPUTE FISHER
    ###########################################

    if Pk_fisher:

        model_id = 'COS3_hz'
        developer = 'CDS'
        derive_model = []
        derive_Pkline = ['Tfs8','Tbs8','Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0,2]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True
        a1_prior = 0.03
        a_t_prior = 0.2

        pf.fisher_Pkline(model_id, 
                developer, 
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed, \
                a1_prior,
                a_t_prior)


if __name__ == "__main__":
    main()