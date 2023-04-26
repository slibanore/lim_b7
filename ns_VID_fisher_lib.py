# SL: last update 01/22/2023

from LIM_b7 import *
from LIM_b7.ns_fiducial_pars import astrocosmo_dict

save_fig_dir = './results_ns/'#deep/'

create_dir(save_fig_dir)
create_dir(save_fig_dir + 'VID/')
create_dir(save_fig_dir + 'fisher_VID/')

direct_id = lambda developer, model_id, z: save_fig_dir + 'VID/' + developer + '/' + model_id + '/z' + str(z)

direct_fisher_id = lambda developer, model_id, z: save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + '/z' + str(z)


def set_der_par(derive_par,z,developer,obs_pars):

    print('\nDoing derivative wrt: ' + str(derive_par))
    # arxiv: 1807.06211 tab 2 TE+lowE
    if derive_par == 'ns': 	
        ns_low = ns_fid*(1. - delta)
        ns_up = ns_fid*(1. + delta)
    else:
        ns_low = ns_fid
        ns_up = ns_fid

    if derive_par == 'nrun': 	
        nrun_low = nrun_fid*(1. - delta)
        nrun_up = nrun_fid*(1. + delta)
    else:
        nrun_low = nrun_fid
        nrun_up = nrun_fid

    if derive_par == 'nrunrun': 	
        nrunrun_low = nrunrun_fid*(1. - delta)
        nrunrun_up = nrunrun_fid*(1. + delta)
    else:
        nrunrun_low = nrunrun_fid
        nrunrun_up = nrunrun_fid

    if derive_par == 'astro_alpha': 	
        use_astro_alpha_low = astro_alpha_fid*(1. - delta)
        use_astro_alpha_up = astro_alpha_fid*(1. + delta)
    else:
        use_astro_alpha_low = astro_alpha_fid 
        use_astro_alpha_up = astro_alpha_fid

    if derive_par == 'astro_beta': 
        use_astro_beta_low = astro_beta_fid*(1. - delta)
        use_astro_beta_up = astro_beta_fid*(1. + delta)
    else:
        use_astro_beta_low = astro_beta_fid 
        use_astro_beta_up = astro_beta_fid

    if derive_par == 'astro_sig_SFR': 	
        use_astro_sig_SFR_low = astro_sig_SFR_fid*(1. - delta)
        use_astro_sig_SFR_up = astro_sig_SFR_fid*(1. + delta)
    else:
        use_astro_sig_SFR_low = astro_sig_SFR_fid 
        use_astro_sig_SFR_up = astro_sig_SFR_fid

    if derive_par == 'astro_sig_scatter': 	
        use_astro_sig_scatter_low = astro_sig_scatter_fid*(1. - delta)
        use_astro_sig_scatter_up = astro_sig_scatter_fid*(1. + delta)
    else:
        use_astro_sig_scatter_low = astro_sig_scatter_fid 
        use_astro_sig_scatter_up = astro_sig_scatter_fid

    if derive_par == 'astro_Lcut': 	
        use_astro_Lcut_low = astro_Lcut_fid*(1. - delta)
        use_astro_Lcut_up = astro_Lcut_fid*(1. + delta)
    else:
        use_astro_Lcut_low = astro_Lcut_fid 
        use_astro_Lcut_up = astro_Lcut_fid

    if derive_par == 'A_tinker':
        use_Atink_low = A_tinker_fid(z)*(1. - delta)
        use_Atink_up = A_tinker_fid(z)*(1. + delta)
    else:
        use_Atink_low = A_tinker_fid(z)
        use_Atink_up = A_tinker_fid(z)

    if derive_par == 'a_tinker':
        use_atink_low = a_tinker_fid(z)*(1. - delta)
        use_atink_up = a_tinker_fid(z)*(1. + delta)
    else:
        use_atink_low = a_tinker_fid(z)
        use_atink_up = a_tinker_fid(z)

    if derive_par == 'b_tinker':
        use_btink_low = b_tinker_fid(z)*(1. - delta)
        use_btink_up = b_tinker_fid(z)*(1. + delta)
    else:
        use_btink_low = b_tinker_fid(z)
        use_btink_up = b_tinker_fid(z)

    if derive_par == 'c_tinker':
        use_ctink_low = c_tinker_fid(z)*(1. - delta)
        use_ctink_up = c_tinker_fid(z)*(1. + delta)
    else:
        use_ctink_low = c_tinker_fid(z)
        use_ctink_up = c_tinker_fid(z)

    low_cosmo = dict(
	f_NL=0, H0=67.67, cosmomc_theta=None,
	ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
	num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
	YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
	TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
	theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
	dark_energy_model='ppf',As=2.105e-09, 
    # !!! Planck 2018 eq. 16/17/18 arXiv:1807.06211 
    ns=ns_low, nrun=nrun_low, nrunrun=nrunrun_low, 
    r=0.0, nt=None, ntrun=0.0, 
	pivot_scalar=0.05, pivot_tensor=0.05,
	parameterization=2,halofit_version='mead')

    low_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_low, \
            beta=use_astro_beta_low, dMF = astro_dMF_fid, \
            sig_SFR=use_astro_sig_SFR_low,SFR_file=SFR_file_fid)}

    low_hmf = dict(A_tinker = use_Atink_low, a_tinker = use_atink_low, b_tinker = use_btink_low, c_tinker = use_ctink_low)

    low_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, cosmo_input_camb = low_cosmo,
        hmf_pars = low_hmf,
        model_par = low_model_par,
        sigma_scatter = use_astro_sig_scatter_low,
        dndL_Lcut=use_astro_Lcut_low,)}

    up_cosmo = dict(
	f_NL=0, H0=67.67, cosmomc_theta=None,
	ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
	num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
	YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
	TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
	theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
	dark_energy_model='ppf',As=2.105e-09, 
    # !!! Planck 2018 eq. 16/17/18 arXiv:1807.06211 
    ns=ns_up, nrun=nrun_up, nrunrun=nrunrun_up, 
    r=0.0, nt=None, ntrun=0.0, 
	pivot_scalar=0.05, pivot_tensor=0.05,
	parameterization=2,halofit_version='mead')

    up_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_up, \
            beta=use_astro_beta_up, dMF = astro_dMF_fid,\
            sig_SFR=use_astro_sig_SFR_up,SFR_file=SFR_file_fid)}

    up_hmf = dict(A_tinker = use_Atink_up, a_tinker = use_atink_up, b_tinker = use_btink_up, c_tinker = use_ctink_up)

    up_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, cosmo_input_camb = up_cosmo,
        hmf_pars = up_hmf,
        model_par = up_model_par,
        sigma_scatter = use_astro_sig_scatter_up,
        dndL_Lcut=use_astro_Lcut_up,)}

    low = update_VID(obs_pars(z),low_data)[0]
    up = update_VID(obs_pars(z),up_data)[0]

    return low, up


def fisher_VID(model_id, developer, \
        derive_pars = ['ns','nrun','nrunrun','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut','a_tinker'],\
    	save_VID_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True, \
	    ns_prior = 0.0056, \
        a_t_prior = 0.2):


    create_dir(save_fig_dir + 'VID/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_VID/' + developer + '/')

    create_dir(save_fig_dir + 'VID/' + developer + '/' + model_id + '/')
    create_dir(save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + '/')

    survey_name, model_name, obs_pars, model_z, z_vals  = set_survey(model_id,developer,estimator='VID')

    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    for z in z_vals:

        fid_par = get_fiducials(derive_pars, developer, z)[0]

        direct = save_fig_dir + 'VID/' + developer + '/' + model_id + '/z' + str(z)
        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + '/z' + str(z)

        create_dir(direct)
        create_dir(direct_fisher)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')

        model = model_z(z)

        N_z_shells = int(model.Delta_nu.value / (model.dnu.value * 1e-3))

        filename = direct + 'Bi_fid'
        if os.path.exists(filename) and import_allowed:
            fiducial_Ti, fiducial_Bi = import_VID(filename, 'Import Bi fiducial', model, max_Ti, min_Ti, True)
        else:		
            print('Running Bi fiducial')

            Ti_full, Bi_full = model.get_VID(\
            Tbin = False, T = model.T, Tmean = model.Tmean,\
            PT = model.PT,\
            PT_zero=model.PT_zero,\
            Nbin = Npt_i, minBi = min_Bi)

            if save_VID_flag:
                save_VID(filename,Ti_full,Bi_full)

            fiducial_Ti, fiducial_Bi =  cut_T(model, Ti_full, Bi_full, max_Ti, min_Ti)

        # FLUCTUATION NOISE
        sigma2_i = np.zeros(len(fiducial_Bi)) 
        for i in range(len(sigma2_i)):
            if error_type == 'poisson':
                sigma2_i[i] = float(fiducial_Bi[i])
            elif error_type == 'neg_bin':
                sigma2_i[i] = float(fiducial_Bi[i] * (1. + fiducial_Bi[i] / ( model.Nvox)))
            else:
                print('Check error type')
                return
                                        
        der_Bi = []

        ################################################	
        # DERIVATIVES PARAMETERS
        ################################################	
        for i in range(len(derive_pars)):

            filename_low = direct + 'Bi_low_' + str(derive_pars[i])			
            filename_up = direct + 'Bi_up_' + str(derive_pars[i])
            
            model_low, model_up = set_der_par(derive_pars[i],z,developer,obs_pars)

            if os.path.exists(filename_low) and import_allowed:
                Bi_low = import_VID(filename_low, 'Import Bi low', model_low, max_Ti, min_Ti, False)[1]
            else:
                print('Running Bi low')

                Ti_low_full, Bi_low_full = model_low.get_VID(\
                Tbin = fiducial_Ti, T = model_low.T, Tmean = model_low.Tmean,\
                PT = model_low.PT,\
                PT_zero=model_low.PT_zero,\
                Nbin = Npt_i, minBi = min_Bi)

                if save_VID_flag:
                    save_VID(filename_low,Ti_low_full,Bi_low_full)

                Bi_low = Bi_low_full

            if os.path.exists(filename_up) and import_allowed:
                Bi_up = import_VID(filename_up, 'Import Bi up', model_up, max_Ti, min_Ti, False)[1]

            else:
                print('Running Bi up')

                Ti_up_full, Bi_up_full = model_up.get_VID(\
                Tbin = fiducial_Ti, T = model_up.T, Tmean = model_up.Tmean,\
                PT = model_up.PT,\
                PT_zero=model_up.PT_zero,\
                Nbin = Npt_i, minBi = min_Bi)

                if save_VID_flag:
                    save_VID(filename_up,Ti_up_full,Bi_up_full)

                Bi_up = Bi_up_full

            der_Bi.append([])
            for t in range(len(fiducial_Ti)):
                if derive_pars[i] != 'fNL':
                    try:
                        der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                            (2*delta*fid_par[i].value)))
                    except:
                        der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                            (2*delta*fid_par[i])))
                else:
                    try:
                        der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                            (20*delta*fid_par[i].value)))
                    except:
                        der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                            (20*delta*fid_par[i])))

        Fisher_Ti = []
        Fisher_z = np.zeros((len(der_Bi),len(der_Bi)))
        for t in range(len(fiducial_Ti)):
            Fisher_Ti.append(np.zeros((len(der_Bi),len(der_Bi))))

        for i in range(len(der_Bi)):
            for j in range(len(der_Bi)):
                for t in range(len(fiducial_Ti)):
                    Fisher_Ti[t][i,j] = der_Bi[i][t]*der_Bi[j][t] / (sigma2_i[t]) if sigma2_i[t] != 0. else 0.

        Fisher_Ti = np.asarray(Fisher_Ti)
        for i in range(len(der_Bi)):
            for j in range(len(der_Bi)):
                for t in range(len(fiducial_Ti)):
                    Fisher_z[i,j] +=  Fisher_Ti[t][i,j]
        Fisher += N_z_shells * Fisher_z
		 	
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + 'fisher_tot' 
        save_fisher(filename, derive_pars,Fisher)

    prior_matrix = np.zeros((len(derive_pars),len(derive_pars)))
		
    if ns_prior: 
        for i in range(len(derive_pars)):
            if 'ns' in derive_pars[i]:
                prior_matrix[i,i] = ns_prior**-2       
        Fisher += prior_matrix

    if a_t_prior:
        for i in range(len(derive_pars)):
            if 'a_tinker' in derive_pars[i]:
                prior_matrix[i,i] = (a_t_prior*a_tinker_fid(0))**-2
        Fisher += prior_matrix

    print('Fisher with a_tinker prior')

    inv_Fisher = np.linalg.inv(Fisher)

    sigma_tot = np.zeros(len(der_Bi))
    perc_error = np.zeros(len(der_Bi))
    for i in range(len(der_Bi)):
        sigma_tot[i] = np.sqrt(inv_Fisher[i,i])
        perc_error[i] = round(sigma_tot[i] / fid_par[i] * 100,2)

    print('Percentage errors: ' + str(perc_error) + '%')

    if plot_sigma_flag:

        print('Doing contour plots')
        names = get_fiducials(derive_pars, developer, z_vals[1])[1]

        plot_dist = GaussianND(fid_par, inv_Fisher, names=names)			

        settings = gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = 35
        settings.legend_fontsize = 30
        settings.lab_fontsize = 30

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14

        #g.triangle_plot([plot_dist], names, filled=True,\
        #    legend_labels=[r'$B_i\ \rm{constraints,\ COS3\ high}-z$'],legend_loc='upper right',\
        #    contour_colors=['darkblue'], line_args=[{'lw':2, 'color':'darkblue'}],
        #    markers={'x2':0})

        g.triangle_plot([plot_dist], names, filled=True,\
            legend_labels=[r'$B_i\ \rm{constraints} $'],legend_loc='upper right',\
            contour_colors=['darkblue'], line_args=[{'lw':2, 'color':'darkblue'}],
            markers={'x2':0})

        plt.savefig(save_fig_dir + 'ellipse_' + model_id + '_prior.pdf')
        plt.show()

    for i in range(len(fid_par)):
        if type(fid_par[i]) is not np.float64 and type(fid_par[i]) is not float:
            fid_par[i] = fid_par[i].value

    print('Marginal errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fid_par * 100)

    return Fisher, sigma_tot