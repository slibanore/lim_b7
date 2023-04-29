from LIM_b7 import *
from LIM_b7.fiducial_pars_pmf import astrocosmo_dict
from datetime import datetime
 
save_fig_dir = './results/PMF/' 

create_dir(save_fig_dir)

direct_id = lambda n_B, sigma_B_0: save_fig_dir + str(round(n_B,3)) + ',' + str(round(sigma_B_0,3)) + '/'

direct_id_model = lambda n_B, sigma_B_0, model_id: direct_id(n_B, sigma_B_0) + model_id + '/' 

direct_id_VID = lambda n_B, sigma_B_0, model_id, z: direct_id_model(n_B, sigma_B_0, model_id) + 'z' + str(z) + '/'


# RUN THE ANALYSIS 
def fisher_VID_pmf(n_B, sigma_B_0, smooth_scale,
        model_id = 'COS3_hz', \
        derive_pars = ['n_B','sigma_B_0','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut'],
        #,'a_tinker'],\
    	save_VID_flag = False, \
        save_fisher_flag = False, \
	    plot_VID_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True,
        prior_matrix = False):
    
    developer = 'PMF'

    create_dir(direct_id(n_B,sigma_B_0))
    create_dir(direct_id_model(n_B,sigma_B_0,model_id))

    model_data = dict(\
        developer = developer,
        n_B = n_B,         
        sigma_B_0 = sigma_B_0,    
        smooth_scale = smooth_scale)

    set_model_data = lambda z: {**astrocosmo_dict(model_data['developer'],z), **model_data}

    if model_id == 'COS3_lz':
        z_vals = [2.6, 3]
        obs_pars = lambda z: deepcopy(obs_params_lowz(z,'COS3'))	
        model_z = lambda z: update_VID(obs_pars(z), set_model_data(z))[0]

    elif model_id == 'COS3_hz':
        z_vals = [7.25] #5.3, 6.25, 7.25]
        obs_pars = lambda z: deepcopy(obs_params_highz(z,'COS3'))	
        model_z = lambda z: update_VID(obs_pars(z), set_model_data(z))[0]

    else:
        print('Check the model input!')
        return
    
    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    for z in z_vals:

        fid_par = get_fiducials_pmf(n_B,sigma_B_0,derive_pars, developer, z)[0]

        direct_VID_z = direct_id_VID(n_B,sigma_B_0,model_id,z)
        create_dir(direct_VID_z)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')

        model = model_z(z)

        N_z_shells = int(model.Delta_nu.value / (model.dnu.value * 1e-3))

        filename = direct_VID_z + 'Bi_fid'
        if os.path.exists(filename) and import_allowed:
            fiducial_Ti, fiducial_Bi = import_VID(filename, 'Import Bi fiducial', model, max_Ti, min_Ti, True)

            # for i in range(len(fiducial_Ti)):
            #    if fiducial_Bi[i] < 0:
            #        fiducial_Bi[i] = 0.

            # plt.figure(z_vals.index(z))
            # plt.loglog(fiducial_Ti, fiducial_Bi,label=r'$n_B=%g$'%n_B)
            # #plt.loglog(fiducial_Ti, fiducial_Bi,label=r'$\sigma=%g$'%sigma_B_0)
            # plt.legend()
            # #plt.title(r'$z = %g,\,$'%z + r'$n_B = %g$'%n_B)
            # plt.title(r'$z = %g,\,$'%z + r'$\sigma_B = %g$'%sigma_B_0)

        else:		
            print('Running Bi fiducial')

            return
            Ti_full, Bi_full = model.get_VID(\
            Tbin = False, T = model.T + model.Tmean, Tmean = model.Tmean,\
            PT = model.PT,\
            PT_zero=model.PT_zero,\
            Nbin = Npt_i, minBi = min_Bi)

            if save_VID_flag:
                save_VID(filename,Ti_full,Bi_full)

            fiducial_Ti, fiducial_Bi =  cut_T(model, Ti_full, Bi_full, max_Ti, min_Ti)

            if plot_VID_flag:

                print('Saving fiducial plot')

                plot_fiducial(model, model.T, model.PT, fiducial_Ti, fiducial_Bi, n_B, sigma_B_0, z , model_id, False)
                if z == z_vals[0]:
                    plot_fiducial(model, model.T, model.PT, fiducial_Ti, fiducial_Bi, n_B, sigma_B_0, z , model_id, True)

        # FLUCTUATION NOISE
        sigma2_i = np.zeros(len(fiducial_Bi)) 
        for i in range(len(sigma2_i)):
            sigma2_i[i] = float(fiducial_Bi[i]) 
                                        
        der_Bi = []

        ################################################	
        # DERIVATIVES PARAMETERS
        ################################################	
        for i in range(len(derive_pars)):
            
            print('Arrived here at time: ' + str(datetime.now()))
            filename_low = direct_VID_z + 'Bi_low_' + str(derive_pars[i])			
            filename_up = direct_VID_z + 'Bi_up_' + str(derive_pars[i])
            
            model_low, model_up = set_der_par(n_B,sigma_B_0,derive_pars[i],z,developer,obs_pars)

            if os.path.exists(filename_low) and import_allowed:
                Bi_low = import_VID(filename_low, 'Import Bi low', model_low, max_Ti, min_Ti, False)[1]

                # for bb in range(len(fiducial_Ti)):
                #     if Bi_low[bb] < 0:
                #         Bi_low[bb] = 0.
                
            else:
                print('Running Bi low')
                return 
                Ti_low_full, Bi_low_full = model_low.get_VID(\
                Tbin = fiducial_Ti, T = model_low.T + model_low.Tmean, Tmean = model_low.Tmean,\
                PT = model_low.PT,\
                PT_zero=model_low.PT_zero,\
                Nbin = Npt_i, minBi = min_Bi)


                if save_VID_flag:
                    save_VID(filename_low,Ti_low_full,Bi_low_full)

                Bi_low = Bi_low_full

            if os.path.exists(filename_up) and import_allowed:
                Bi_up = import_VID(filename_up, 'Import Bi up', model_up, max_Ti, min_Ti, False)[1]

                # for bb in range(len(fiducial_Ti)):
                #     if Bi_up[bb] < 0:
                #         Bi_up[bb] = 0.

            else:
                print('Running Bi up')

                return 
                Ti_up_full, Bi_up_full = model_up.get_VID(\
                Tbin = fiducial_Ti, T = model_up.T + model_up.Tmean, Tmean = model_up.Tmean,\
                PT = model_up.PT,\
                PT_zero=model_up.PT_zero,\
                Nbin = Npt_i, minBi = min_Bi)

                if save_VID_flag:
                    save_VID(filename_up,Ti_up_full,Bi_up_full)

                Bi_up = Bi_up_full

            der_Bi.append([])
            for t in range(len(fiducial_Ti)):
                if derive_pars[i] != 'fNL':
                    if fid_par[i] != 0.:
                        try:
                            der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                                (2*delta*fid_par[i].value)))
                        except:
                            der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                                (2*delta*fid_par[i])))
                    else:
                        try:
                            der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                                (delta)))
                        except:
                            der_Bi[i].append(float((Bi_up[t] - Bi_low[t]) / \
                                (delta)))
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
                    Fisher_Ti[t][i,j] = der_Bi[i][t]*der_Bi[j][t] / (sigma2_i[t]) if sigma2_i[t] == 0. else 0.

        Fisher_Ti = np.asarray(Fisher_Ti)
        for i in range(len(der_Bi)):
            for j in range(len(der_Bi)):
                for t in range(len(fiducial_Ti)):
                    Fisher_z[i,j] +=  Fisher_Ti[t][i,j]
        Fisher += N_z_shells * Fisher_z
		 	
    if save_fisher_flag:
        filename = save_fig_dir + str(round(n_B,3)) + ',' + str(round(sigma_B_0,3)) + '_TEST_fisher_' + model_id 
        save_fisher(filename, derive_pars, Fisher)

    if prior_matrix:
        Fisher += prior_matrix

    inv_Fisher = np.linalg.inv(Fisher)

    sigma_tot = np.zeros(len(der_Bi))
    perc_error = np.zeros(len(der_Bi))
    print('\n---------------------------------------------------\nMarginalized errors:')
    for i in range(len(der_Bi)):
        sigma_tot[i] = np.sqrt(inv_Fisher[i,i])
        perc_error[i] = round(sigma_tot[i] / fid_par[i] * 100,2)
        if round(fid_par[i],3) != 0.:
            print(str(derive_pars[i]) + ': ' + str(round(fid_par[i],3)) + ' +- ' + str(round(sigma_tot[i],3)) + ' (' + str(round(perc_error[i],3)) + '%)')
        else:
            print(str(derive_pars[i]) + ': ' + str(fid_par[i]) + ' +- ' + str(sigma_tot[i]) + ' (' + str(perc_error[i]) + '%)')
    print('----------------------------------------------------\n')

    if save_fisher_flag:
        filename = save_fig_dir + str(round(n_B,3)) + ',' + str(round(sigma_B_0,3)) + '_TEST_sigmapar_' + model_id 
        save_sigmapar(filename, derive_pars, fid_par, sigma_tot)

    if plot_sigma_flag:

        print('Doing contour plots')
        names = get_fiducials_pmf(n_B,sigma_B_0,derive_pars, developer, z_vals[1])[1]

        plot_dist = GaussianND(fid_par, inv_Fisher, names=names)			

        settings = gp.GetDistPlotSettings()

        settings.fig_width_inch = 10
        settings.norm_prob_label = False
        settings.axes_fontsize = 30
        settings.legend_fontsize = 30
        settings.lab_fontsize = 30

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add = 0.4
        g.settings.title_limit_fontsize = 14

        model_label = r'$,\,{\rm COS3\, high\,} z$' if model_id == 'COS3_hz' else r'$,\,{\rm COS3\, low\,} z$' if model_id == 'COS3_lz' else ''

        g.triangle_plot([plot_dist], names, filled=True,\
            legend_labels=[r'$B_i$' + 
            model_label],legend_loc='upper right',\
            contour_colors=[suoko], line_args=[{'lw':2, 'color':suoko}],
            markers={'x2':0})

        plt.tight_layout()
        plt.savefig(save_fig_dir + str(round(n_B,3)) + ',' + str(round(sigma_B_0,3)) + '_ellipse_' + model_id + '.pdf')

    return


############################################
############################################
############################################

###################################
##       FUNCTIONS YOU NEED      ##
###################################

# DERIVATIVES
def set_der_par(n_B,sigma_B_0,derive_par,z,developer,obs_pars):

    print('\nDoing derivative wrt: ' + str(derive_par))
    
    if derive_par == 'n_B': 	
        use_n_B_low = n_B*(1. - delta)
        use_n_B_up = n_B*(1. + delta)
    else:
        use_n_B_low = n_B 
        use_n_B_up = n_B
    
    if derive_par == 'sigma_B_0':
        use_sigma_B_0_low = sigma_B_0*(1. - delta)
        use_sigma_B_0_up = sigma_B_0*(1. + delta)
        # TalAdi: in case we use sigma_B_0 = 0 as fiducial
        if sigma_B_0_fid == 0.:
            use_sigma_B_0_up = delta
            use_sigma_B_0_low = sigma_B_0
    else:
        use_sigma_B_0_low = sigma_B_0 
        use_sigma_B_0_up = sigma_B_0
    
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
        use_fNL_low = fNL_fid*(1. - 10*delta)
        use_fNL_up = fNL_fid*(1. + 10*delta)
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
        hmf_pars = low_hmf,
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
        hmf_pars = up_hmf,
        model_par = up_model_par,
        sigma_scatter = use_astro_sig_scatter_up,
        dndL_Lcut=use_astro_Lcut_up,
        n_B = use_n_B_up,
        sigma_B_0 = use_sigma_B_0_up)}

    low = update_VID(obs_pars(z),low_data)[0]
    up = update_VID(obs_pars(z),up_data)[0]

    return low, up


# FIDUCIAL PLOT
def plot_fiducial(model, T_fid, pT_fid, Ti_fid, Bi_fid, n_B, sigma_B_0, redshift, model_id, save_gen_figure):

	pk_fid = model.Pm[0]
	k_fid = model.k

	n_fid = model.dndM
	M_fid = model.M

	##########################################
	# plot procedure
	##########################################

	plt.figure()

	outer = gridspec.GridSpec(\
		2, 2,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	ax1 = plt.subplot(outer[0])
	ax2 = plt.subplot(outer[1])
	ax3 = plt.subplot(outer[2])
	ax4 = plt.subplot(outer[3])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)
	ax3.get_yaxis().set_label_coords(-0.12,0.5)
	ax4.get_yaxis().set_label_coords(-0.12,0.5)
			

	##########################################
	# 1) plot the matter power spectrum
	##########################################

	label = r'$n_B = %g$'%n_B + r'$,\,\sigma_{B,0} = %g$'%sigma_B_0 

	ax1.loglog(k_fid, pk_fid, linestyle='-',label = label, color = aonibi)
	
	ax1.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	ax1.set_ylabel(r'$P(k,z=%g)\, {\rm [Mpc^{3}]}$'%redshift)
	ax1.legend(loc=3, ncol = 1,frameon=False)
	ax1.set_ylim(3e-5,3e4)
	ax1.set_xlim(k_fid[0].value-0.1,k_fid[-1].value+50)

	##########################################
	# 2) plot the halo mass function
	##########################################

	ax2.loglog(M_fid, n_fid, color = aonibi,linestyle='-')

	ax2.set_xlabel(r'$M_h\, {\rm [}M_{\odot}{\rm ]}$')
	ax2.set_ylabel(r'$\frac{dn}{dM_{h}}\, {\rm [Mpc^{-3}}\,M_{\odot} {\rm ]}$')
	ax2.set_ylim(1.5e-25,1e-8)

	##########################################
	# 3) plot PT
	##########################################

	ax3.loglog(T_fid+model.Tmean, pT_fid,color = aonibi,linestyle='-')
	
	ax3.set_xlabel(r'$T\, {\rm [\mu K]}$')
	ax3.set_ylabel(r'$\mathcal{P}(T)$')
	ax3.set_xlim(model.Tmin_VID,model.Tmax_VID)
	#ax3.set_ylim(8e-7,5e1)

	##########################################
	# 4) plot VID
	##########################################

	ax4.loglog(Ti_fid,Bi_fid,color = aonibi,linestyle='-')
	
	ax4.set_xlabel(r'$T_i\, {\rm [\mu K]}$')
	ax4.set_ylabel(r'$B_i(z=%g)$'%redshift)
	ax4.set_xlim(min_Ti,max(Ti_fid.value)*2.)
	ax4.set_ylim(min(Bi_fid)/2.,max(Bi_fid)*2.)
		
	##########################################

	plt.savefig(direct_id_VID(n_B,sigma_B_0,model_id,redshift) + 'fid_plot.pdf')

	if save_gen_figure:		
		plt.savefig(save_fig_dir + str(round(n_B,3)) + ',' + str(round(sigma_B_0,3)) + '_fidplot_' + model_id + '.pdf')

	return 

# SAVE MARGINALIZED ERRORS
def save_sigmapar(filename, derive_pars, fid_par, sigma):

    f = open(filename,'w')
    f.write('parameter\tfiducial\tsigma_abs\tsigma_perc\n\n')
    for i in range(len(derive_pars)):
        f.write(str(derive_pars[i]) + '\t')
        f.write(str(fid_par[i]) + '\t')
        f.write(str(sigma[i]) + '\t')
        f.write('(' + str(round(sigma[i]/fid_par[i]*100,2)) + '%)'+ '\n')
    f.close()

    return 


def get_fiducials_pmf(n_B, sigma_B_0, derive_pars, developer, z):

	fiducials = []
	names = []
	wrong_pars = False

	for i in range(len(derive_pars)):
						
		if derive_pars[i] == 'astro_alpha':
			fiducials.append(astro_alpha_fid)
			names.append(r'$\alpha$')
			
		if derive_pars[i] == 'astro_beta':
			fiducials.append(astro_beta_fid)
			names.append(r'$\beta$')
			
		if derive_pars[i] == 'astro_sig_SFR':
			fiducials.append(astro_sig_SFR_fid)
			names.append(r'$\sigma_{SFR}$')
			
		if derive_pars[i] == 'astro_sig_scatter':
			fiducials.append(astro_sig_scatter_fid)
			names.append(r'$\sigma_{S}$')
			
		if derive_pars[i] == 'astro_Lcut':
			fiducials.append(astro_Lcut_fid.value)
			names.append(r'$L_{cut}$')
			
		if derive_pars[i] == 'A_tinker':
			fiducials.append(A_tinker_fid(z))
			names.append(r'$A^T$')
			
		if derive_pars[i] == 'a_tinker':
			fiducials.append(a_tinker_fid(z))
			names.append(r'$a^T$')
			
		if derive_pars[i] == 'b_tinker':
			fiducials.append(b_tinker_fid(z))
			names.append(r'$b^T$')
			
		if derive_pars[i] == 'c_tinker':
			fiducials.append(c_tinker_fid(z))
			names.append(r'$c^T$')

		if derive_pars[i] == 'n_B':
			fiducials.append(n_B)
			names.append(r'$n_B$')

		if derive_pars[i] == 'sigma_B_0':
			fiducials.append(sigma_B_0)
			names.append(r'$\sigma_{B,0}$')


		if derive_pars[i] == 'fNL':
			fiducials.append(fNL_fid)
			names.append(r'$f_{NL}$')

	if wrong_pars:
		print('\n-------------------------')
		print('Some parameter in the list are not allowed in the model you are using')
		print('-------------------------\n')
		return 
	
	return fiducials, names
