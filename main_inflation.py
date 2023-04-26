# SL: last update 01/28/2023

import plot_lib as vp 
#!!!
import ns_VID_fisher_lib as vf
import ns_Pkline_fisher as pf
from LIM_b7.ns_fiducial_pars import *


def main():

    ###########################################
    # FIDUCIALS
    ###########################################

    # 1) set up the fiducial parameters you want to change wrt the default file
    # NOTE: before changing the hmf or astro models check the syntax in fiducial_pars.py / default_pars.py

    from LIM_b7 import update_fiducials

    par_list = []#'hmf_fid', 
               # 'hmf_pars_fid',
               # 'fNL_fid']
    par_value = []#str("'NG_Riotto'"),
                #    'lambda z: dict(A_tinker = A_tinker_fid(z),a_tinker = a_tinker_fid(z),b_tinker = b_tinker_fid(z),c_tinker = c_tinker_fid(z),fNL = fNL_fid)',
                #    1.
                #]

    file_pars = 'LIM_b7/ns_fiducial_pars.py'

    if len(par_list) > 0:
        update_fiducials(file_pars,par_list,par_value)

    ###########################################
    # CHOOSE MODULE
    ###########################################

    # 2) choose which module to run 

    check_plots = False
    VID_fisher = False
    Pk_fisher = True

    ###########################################
    # SETUP ANALYSIS
    ###########################################

    # 3) set up specific analysis parameters
    # NOTE that each module has inside more setting that can be tuned 

    if check_plots:

        mod_par = 'nrunrun'
        # mod_val = [lambda z: dict(A_tinker = A_tinker_fid(z),a_tinker = a_tinker_fid(z),b_tinker = b_tinker_fid(z),c_tinker = c_tinker_fid(z),fNL = 100.),
        #            lambda z: dict(A_tinker = A_tinker_fid(z),a_tinker = a_tinker_fid(z),b_tinker = b_tinker_fid(z),c_tinker = c_tinker_fid(z),fNL = 10.),
        #            lambda z: dict(A_tinker = A_tinker_fid(z),a_tinker = a_tinker_fid(z),b_tinker = b_tinker_fid(z),c_tinker = c_tinker_fid(z),fNL = 5.)]

        mod_val = [0.005, 0.010, 0.020]
        save_figure = True
        plot_bias = False
        if type(mod_val) == list:
            get_SNR = True
            vp.compare_multi_vals(mod_par, mod_val, save_figure, get_SNR)
        else:
            vp.compare_models(mod_par, mod_val, save_figure)
        if plot_bias:
            vp.compare_scaledependence(mod_par,mod_val,save_figure)

    if VID_fisher:

        model_id = 'COS3_hz'
        developer = 'CDS'
        derive_pars = ['ns','nrun','nrunrun']#,'astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut']
        #derive_pars = ['fNL','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut']
        save_VID_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True
        ns_prior = 0.0056
        a_t_prior = False #0.2

        vf.fisher_VID(model_id, 
            developer, 
            derive_pars,
            save_VID_flag, 
            save_fisher_flag, 
            plot_sigma_flag, 
            import_allowed, 
            ns_prior, 
            a_t_prior)


    if Pk_fisher:

        model_id = 'COS3_hz'
        developer = 'CDS'
        #derive_model = ['fNL']#,'a1','a2','a3']
        derive_model = ['ns','nrun','nrunrun']
        derive_Pkline = []#'Tfs8','Tbs8','Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0,2]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = False
        import_allowed = True
        ns_prior = False
        a_t_prior = False #0.2

        pf.fisher_Pkline(model_id, 
                developer, 
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed, \
                ns_prior,
                a_t_prior)


if __name__ == "__main__":
    main()