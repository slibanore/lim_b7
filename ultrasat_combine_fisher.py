# !!!
import ultrasat_VID_fisher_lib as vid
import ultrasat_Pkline_fisher as pk

def sum_fisher(z = 1.135, par_comm = ['fNL','A_lya','B_lya','D_lya'], par_vid = ['astro_sig_scatter','astro_Lcut'], par_pk = ['Pshot','sNL','alpha_par','alpha_perp'],plot_sigma = True):
 
    all_par_vid = par_comm + par_vid
    all_par_pk = par_comm + par_pk
    
    F_vid = vid.fisher_VID(z, 'CDS', \
        derive_pars = all_par_vid,\
    	save_VID_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True)[0]

    F_pk = pk.fisher_Pkline(z, 'CDS', \
        derive_model = par_comm,\
        derive_Pkline = par_pk,	
        multipole = [0,2,4],\
    	save_Pkline_flag = True, \
        save_fisher_flag = True, \
	    plot_sigma_flag = False, \
        import_allowed = True)[0]
        
    all_par = par_comm + par_vid + par_pk

    F_tot = vid.np.zeros((len(all_par),len(all_par)))
    for i in range(len(all_par)):
        for j in range(len(all_par)):

            if all_par[i] in par_comm:
                vid_i = all_par_vid.index(all_par[i])
                pk_i = all_par_pk.index(all_par[i])
            
                if all_par[j] in par_comm:
                    vid_j = all_par_vid.index(all_par[j])
                    pk_j = all_par_pk.index(all_par[j])
                    F_tot[i,j] = F_vid[vid_i,vid_j] + F_pk[pk_i,pk_j]
                elif all_par[j] in par_vid:
                    vid_j = all_par_vid.index(all_par[j])
                    F_tot[i,j] = F_vid[vid_i,vid_j]
                elif all_par[j] in par_pk:
                    pk_j = all_par_pk.index(all_par[j])
                    F_tot[i,j] = F_pk[pk_i,pk_j]
            
            elif all_par[i] in par_vid:
                vid_i = all_par_vid.index(all_par[i])
                if all_par[j] in par_comm or all_par[j] in par_vid:
                    vid_j = all_par_vid.index(all_par[j])
                    F_tot[i,j] = F_vid[vid_i,vid_j] 
            elif all_par[i] in par_pk:
                pk_i = all_par_pk.index(all_par[i])
                if all_par[j] in par_comm or all_par[j] in par_pk:
                    pk_j = all_par_pk.index(all_par[j])
                    F_tot[i,j] = F_pk[pk_i,pk_j]

    phi_norm = vid.np.zeros((len(F_tot),len(F_tot)))
    for i in range(len(F_tot)):
        #for j in range(len(Fisher)):
        #if i == j:
        phi_norm[i,i] = 10**-(vid.np.log10(abs(F_tot[i,i])))/10

    try:
        new_Fisher = vid.np.linalg.multi_dot([phi_norm,F_tot])
        normalized = True
        print('Fisher normalized')
    except:
        new_Fisher = F_tot
        normalized = False

    start_inv_Fisher = vid.np.linalg.inv(new_Fisher)
    if normalized:
        inv_Ftot = vid.np.linalg.multi_dot([phi_norm,start_inv_Fisher])
    else:
        inv_Ftot = start_inv_Fisher

    inv_Ftot = vid.np.linalg.inv(F_tot)
    sigma_tot = vid.np.zeros(len(all_par))
    for i in range(len(all_par)):
        sigma_tot[i] = vid.np.sqrt(inv_Ftot[i,i])

    print('\nMarginalized errors:')
    print(sigma_tot)

    print('\nMarginalized error on fnl:')
    print(sigma_tot[0])

    if plot_sigma:
        print('Doing contour plots')
        fid_comm, names_comm = vid.get_fiducials(par_comm, 'CDS', z)
        fid_vid, names_vid = vid.get_fiducials(par_vid, 'CDS', z)
        fid_pk, names_pk = pk.get_fiducials_Pkline(par_pk, pk.obs_params_ULTRASAT(z), z)

        fid_par = fid_comm + fid_vid + fid_pk
        for i in range(len(fid_par)):
            try:
                fid_par[i] = fid_par[i].value
            except:
                fid_par[i] = fid_par[i]

        names = names_comm + names_vid + names_pk

        plot_dist = vid.GaussianND(fid_par, inv_Ftot, names=names)			

        settings = vid.gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = 35
        settings.legend_fontsize = 30
        settings.lab_fontsize = 30

        g = vid.gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14

        g.triangle_plot([plot_dist], names, filled=True,\
            legend_labels=[r'$B_i\ \rm{constraints,}\ deep$'],legend_loc='upper right',\
            contour_colors=['darkred'], line_args=[{'lw':2, 'color':'darkred'}],
            markers={'x2':0})

        vid.plt.show()

    return F_tot, F_pk, F_vid