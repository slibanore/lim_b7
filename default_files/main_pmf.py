# SL: last update 01/28/2023

import plot_lib as vp 
import VID_fisher_lib as vf
import Pkline_fisher as pf
from LIM_b7.fiducial_pars import *
import LIM_b7
def main():

    ###########################################
    # FIDUCIALS
    ###########################################

    # 1) set up the fiducial parameters you want to change wrt the default file
    # NOTE: before changing the hmf or astro models check the syntax in fiducial_pars.py / default_pars.py

    # !!! To fix : if this function updates the fiducial pars list but then the following ones are run with the previous version - how to reload it?? 

    from LIM_b7 import update_fiducials

    par_list = []
    par_value = []

    #par_list = ['n_B_fid', 
    #            'sigma_B_0_fid',
    #            'smooth_scale_fid']
    #par_value = [-2.9,
    #            0.,
    #            1
    #            ]

    file_pars = 'LIM_b7/fiducial_pars.py'

    if len(par_list) > 0:
        update_fiducials(file_pars,par_list,par_value)

    ###########################################
    # CHOOSE MODULE
    ###########################################

    # 2) choose which module to run 

    check_plots = False
    VID_fisher = True
    Pk_fisher = True

    ###########################################
    # SETUP ANALYSIS
    ###########################################

    # 3) set up specific analysis parameters
    # NOTE that each module has inside more setting that can be tuned 

    if check_plots:

        mod_par = 'n_B'
        mod_val = [-2,-2.9]
        save_figure = True
        if type(mod_val) == list:
            get_SNR = True
            vp.compare_multi_vals(mod_par, mod_val, save_figure, get_SNR)
        else:
            vp.compare_models(mod_par, mod_val, save_figure)


    if VID_fisher:

        model_id = 'COS3_hz'
        developer = 'PMF'
        derive_pars = ['n_B','sigma_B_0','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut','a_tinker']
        save_VID_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True
        a1_prior = 0.03
        a_t_prior = 0.2

        vf.fisher_VID(model_id, 
            developer, 
            derive_pars,
            save_VID_flag, 
            save_fisher_flag, 
            plot_sigma_flag, 
            import_allowed, 
            a1_prior, 
            a_t_prior)


    if Pk_fisher:

        model_id = 'COS3_hz'
        developer = 'PMF'
        derive_model = ['n_B','sigma_B_0','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','a_tinker']
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