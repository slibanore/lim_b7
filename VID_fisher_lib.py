# SL: last update 01/22/2023 

from LIM_b7 import *
from LIM_b7.fiducial_pars_pmf import astrocosmo_dict

save_fig_dir = './results/PMF/'


create_dir(save_fig_dir)
create_dir(save_fig_dir + 'VID/')
create_dir(save_fig_dir + 'fisher_VID/')

direct_id = lambda developer, model_id, z: save_fig_dir + 'VID/' + developer + '/' + model_id + '/z' + str(z)

direct_fisher_id = lambda developer, model_id, z: save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + '/z' + str(z)


def set_der_par(derive_par,z,developer,obs_pars):

    print('\nDoing derivative wrt: ' + str(derive_par))
    if derive_par == 'a1': 	
        use_CDS_alpha_1_low = CDS_alpha_1_fid*(1. - delta)
        use_CDS_alpha_1_up = CDS_alpha_1_fid*(1. + delta)
    else:
        use_CDS_alpha_1_low = CDS_alpha_1_fid
        use_CDS_alpha_1_up = CDS_alpha_1_fid

    if derive_par == 'a2': 	
        use_CDS_alpha_2_low = CDS_alpha_2_fid*(1. - delta)
        use_CDS_alpha_2_up = CDS_alpha_2_fid*(1. + delta)
    else:
        use_CDS_alpha_2_low = CDS_alpha_2_fid 
        use_CDS_alpha_2_up = CDS_alpha_2_fid

    if derive_par == 'a3': 	
        use_CDS_alpha_3_low = CDS_alpha_3_fid*(1. - delta)
        use_CDS_alpha_3_up = CDS_alpha_3_fid*(1. + delta)
    else:
        use_CDS_alpha_3_low = CDS_alpha_3_fid 
        use_CDS_alpha_3_up = CDS_alpha_3_fid

    if derive_par == 'a4': 	
        use_CDS_alpha_4_low = CDS_alpha_4_fid*(1. - delta)
        use_CDS_alpha_4_up = CDS_alpha_4_fid*(1. + delta)
    else:
        use_CDS_alpha_4_low = CDS_alpha_4_fid 
        use_CDS_alpha_4_up = CDS_alpha_4_fid

    if derive_par == 'f_FDM': 	
        use_f_axion_low = f_axion_fid*(1. - delta) if f_axion_fid != 0. else -delta
        use_f_axion_up = f_axion_fid*(1. + delta) if f_axion_fid != 0. else delta
    else:
        use_f_axion_low = f_axion_fid 
        use_f_axion_up = f_axion_fid
    
    if derive_par == 'n_B': 	
        use_n_B_low = n_B_fid*(1. - delta) if n_B_fid != 0. else -delta
        use_n_B_up = n_B_fid*(1. + delta) if n_B_fid != 0. else delta
    else:
        use_n_B_low = n_B_fid 
        use_n_B_up = n_B_fid
    
    if derive_par == 'sigma_B_0': 	
        use_sigma_B_0_low = sigma_B_0_fid*(1. - delta) if sigma_B_0_fid != 0. else -delta
        use_sigma_B_0_up = sigma_B_0_fid*(1. + delta) if sigma_B_0_fid != 0. else delta
    else:
        use_sigma_B_0_low = sigma_B_0_fid 
        use_sigma_B_0_up = sigma_B_0_fid
    
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

    if hmf_fid == 'NG_Riotto' and derive_par == 'fNL':
        use_fNL_low = fNL_fid*(1. - 10*delta) if fNL_fid != 0. else -10*delta
        use_fNL_up = fNL_fid*(1. + 10*delta) if fNL_fid != 0. else 10*delta
    else:
        use_fNL_low = fNL_fid
        use_fNL_up = fNL_fid


    low_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_low, \
            beta=use_astro_beta_low, dMF = astro_dMF_fid, \
            sig_SFR=use_astro_sig_SFR_low,SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        low_hmf = dict(A_tinker = use_Atink_low, a_tinker = use_atink_low, b_tinker = use_btink_low, c_tinker = use_ctink_low)
    elif hmf_fid == 'NG_Riotto':
        low_hmf = dict(A_tinker = use_Atink_low, a_tinker = use_atink_low, b_tinker = use_btink_low, c_tinker = use_ctink_low,fNL=use_fNL_low)

    low_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, 
        CDS_alpha_1 = use_CDS_alpha_1_low, \
        CDS_alpha_2 = use_CDS_alpha_2_low, 
        CDS_alpha_3 = use_CDS_alpha_3_low, 
        CDS_alpha_4 = use_CDS_alpha_4_low, 
        hmf_pars = low_hmf,
        f_axion = use_f_axion_low ,
        model_par = low_model_par,
        sigma_scatter = use_astro_sig_scatter_low,
        dndL_Lcut=use_astro_Lcut_low,
        n_B = use_n_B_low,
        sigma_B_0 = use_sigma_B_0_low)}

    up_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_up, \
            beta=use_astro_beta_up, dMF = astro_dMF_fid,\
            sig_SFR=use_astro_sig_SFR_up,SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        up_hmf = dict(A_tinker = use_Atink_up, a_tinker = use_atink_up, b_tinker = use_btink_up, c_tinker = use_ctink_up)
    elif hmf_fid == 'NG_Riotto':
        up_hmf = dict(A_tinker = use_Atink_up, a_tinker = use_atink_up, b_tinker = use_btink_up, c_tinker = use_ctink_up, fNL=use_fNL_up)

    up_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, 
        CDS_alpha_1 = use_CDS_alpha_1_up, \
        CDS_alpha_2 = use_CDS_alpha_2_up, 
        CDS_alpha_3 = use_CDS_alpha_3_up, 
        CDS_alpha_4 = use_CDS_alpha_4_up, 
        hmf_pars = up_hmf,
        f_axion = use_f_axion_up ,
        model_par = up_model_par,
        sigma_scatter = use_astro_sig_scatter_up,
        dndL_Lcut=use_astro_Lcut_up,
        n_B = use_n_B_up,
        sigma_B_0 = use_sigma_B_0_up)}

    low = update_VID(obs_pars(z),low_data)[0]
    up = update_VID(obs_pars(z),up_data)[0]

    return low, up


def fisher_VID(model_id, developer, \
        derive_pars = ['a1','a2','a3','a4','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut','a_tinker'],\
    	save_VID_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True, \
	    a1_prior = 0.03, \
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
        if developer == 'axions':
            direct += '_' + str(f_axion_fid) 
        if developer == 'PMF':
            direct += '_' + str(n_B_fid) + ',' + str(sigma_B_0_fid)
        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_VID/' + developer + '/' + model_id + '/z' + str(z)
        if developer == 'axions':
            direct_fisher += '_' + str(f_axion_fid) 
        if developer == 'PMF':
            direct_fisher += '_' + str(n_B_fid) + ',' + str(sigma_B_0_fid)

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
		
    if a1_prior: 
        for i in range(len(derive_pars)):
            if 'a1' in derive_pars[i]:
                prior_matrix[i,i] = a1_prior**-2       
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
            legend_labels=[r'$B_i\ \rm{constraints,\ wide} $'],legend_loc='upper right',\
            contour_colors=['darkblue'], line_args=[{'lw':2, 'color':'darkblue'}],
            markers={'x2':0})

        plt.savefig(save_fig_dir + 'ellipse_' + model_id + '.pdf')
        plt.show()

    for i in range(len(fid_par)):
        if type(fid_par[i]) is not np.float64 and type(fid_par[i]) is not float:
            fid_par[i] = fid_par[i].value

    print('Marginal errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fid_par * 100)

    return Fisher, sigma_tot