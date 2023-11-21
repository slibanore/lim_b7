# SL: last update 01/22/2023

from LIM_b7 import *
from LIM_b7.ultrasat_fiducial_pars import astrocosmo_dict

save_fig_dir = './results/ultrasat/'

if hmf_fid == 'NG_Riotto':
    save_fig_dir += 'fNL' + str(fNL_fid) + '/'
    

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

    '''
    if derive_par == 'z0_Lyalpha':
        use_z0_Lyalpha_low = z0_Lyalpha_fid*(1. - delta)
        use_z0_Lyalpha_up = z0_Lyalpha_fid*(1. + delta)
    else:
        use_z0_Lyalpha_low = z0_Lyalpha_fid
        use_z0_Lyalpha_up = z0_Lyalpha_fid

    if derive_par == 'f0_Lyalpha':
        use_f0_Lyalpha_low = f0_Lyalpha_fid*(1. - delta)
        use_f0_Lyalpha_up = f0_Lyalpha_fid*(1. + delta)
    else:
        use_f0_Lyalpha_low = f0_Lyalpha_fid
        use_f0_Lyalpha_up = f0_Lyalpha_fid

    if derive_par == 't_Lyalpha':
        use_t_Lyalpha_low = t_Lyalpha_fid*(1. - delta)
        use_t_Lyalpha_up = t_Lyalpha_fid*(1. + delta)
    else:
        use_t_Lyalpha_low = t_Lyalpha_fid
        use_t_Lyalpha_up = t_Lyalpha_fid
    '''
    if hmf_fid == 'NG_Riotto' and derive_par == 'fNL':
        use_fNL_low = fNL_fid*(1. - 10*delta) if fNL_fid != 0. else -10*delta
        use_fNL_up = fNL_fid*(1. + 10*delta) if fNL_fid != 0. else 10*delta
    else:
        use_fNL_low = fNL_fid
        use_fNL_up = fNL_fid

    if derive_par == 'A_lya':
        use_A_Lyalpha_low = A_Lyalpha_fid*(1. - delta)
        use_A_Lyalpha_up = A_Lyalpha_fid*(1. + delta)
    else:
        use_A_Lyalpha_low = A_Lyalpha_fid
        use_A_Lyalpha_up = A_Lyalpha_fid

    if derive_par == 'B_lya':
        use_B_Lyalpha_low = B_Lyalpha_fid*(1. - delta)
        use_B_Lyalpha_up = B_Lyalpha_fid*(1. + delta)
    else:
        use_B_Lyalpha_low = B_Lyalpha_fid
        use_B_Lyalpha_up = B_Lyalpha_fid

    if derive_par == 'D_lya':
        use_D_Lyalpha_low = D_Lyalpha_fid*(1. - delta)
        use_D_Lyalpha_up = D_Lyalpha_fid*(1. + delta)
    else:
        use_D_Lyalpha_low = D_Lyalpha_fid
        use_D_Lyalpha_up = D_Lyalpha_fid


    #low_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_low, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_low, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_low,  dndL_Lcut=use_astro_Lcut_low, SFR_file=SFR_file_fid)}
    low_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_low,B_lya=use_B_Lyalpha_low,D_lya=use_D_Lyalpha_low,  dndL_Lcut=use_astro_Lcut_low, SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        low_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z))

    elif hmf_fid == 'NG_Riotto':
        low_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z),
        fNL = use_fNL_low)

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
        n_B = n_B_fid,
        sigma_B_0 = sigma_B_0_fid)}

    #up_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_up, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_up, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_up,  dndL_Lcut=use_astro_Lcut_up, SFR_file=SFR_file_fid)}
    up_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_up,B_lya=use_B_Lyalpha_up,D_lya=use_D_Lyalpha_up,  dndL_Lcut=use_astro_Lcut_up, SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        up_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z))

    elif hmf_fid == 'NG_Riotto':
        up_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z),
        fNL = use_fNL_up)

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
        n_B = n_B_fid,
        sigma_B_0 = sigma_B_0_fid)}

    low = update_VID(obs_pars(z),low_data)[0]
    up = update_VID(obs_pars(z),up_data)[0]

    return low, up


def fisher_VID(z, developer, \
        derive_pars = ['a1','a2','a3','a4','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut','a_tinker'],\
    	save_VID_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):


    create_dir(save_fig_dir + 'VID/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_VID/' + developer + '/')

    create_dir(save_fig_dir + 'VID/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_VID/' + developer + '/')

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))
    #obs_pars = lambda z_val: deepcopy(obs_params_SPHEREx(z_val))

    model = update_VID(obs_pars(z), astrocosmo_dict(developer,z))[0]

    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    fid_par = get_fiducials(derive_pars, developer, z)[0]

    direct = save_fig_dir + 'VID/' + developer + '/z' + str(z)
    if developer == 'axions':
        direct += '_' + str(f_axion_fid) 

    direct += '/'
    direct_fisher = save_fig_dir + 'fisher_VID/' + developer + '/z' + str(z)
    if developer == 'axions':
        direct_fisher += '_' + str(f_axion_fid) 

    create_dir(direct)
    create_dir(direct_fisher)

    print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')

    filename = direct + 'Bi_fid'
    if os.path.exists(filename) and import_allowed:
        fiducial_Ti = []
        fiducial_Bi = []
        f = open(filename,'r')
        lines = f.readlines()
        for x in range(len(lines)-1):
            words = lines[x+1].split('\t')
            fiducial_Ti.append(float(words[0]))
            fiducial_Bi.append(float(words[1]))

        fiducial_Ti *= u.Jy
        fiducial_Ti /= u.sr

    else:		
        print('Running Bi fiducial')

        fiducial_Ti, fiducial_Bi = model.get_VID(\
        Tbin = False, T = model.T, Tmean = model.Tmean,\
        PT = model.PT,\
        PT_zero=model.PT_zero,\
        Nbin = Npt_i, minBi = min_Bi)

        if save_VID_flag:
            save_VID(filename,fiducial_Ti.value,fiducial_Bi)

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

            Bi_low = []
            f = open(filename_low,'r')
            lines = f.readlines()
            for x in range(len(lines)-1):
                words = lines[x+1].split('\t')
                Bi_low.append(float(words[1]))

        else:
            print('Running Bi low')

            Ti_low_full, Bi_low = model_low.get_VID(\
            Tbin = fiducial_Ti, T = model_low.T, Tmean = model_low.Tmean,\
            PT = model_low.PT,\
            PT_zero=model_low.PT_zero,\
            Nbin = Npt_i, minBi = min_Bi)

            if save_VID_flag:
                save_VID(filename_low,Ti_low_full,Bi_low)

        if os.path.exists(filename_up) and import_allowed:
            Bi_up = []
            f = open(filename_up,'r')
            lines = f.readlines()
            for x in range(len(lines)-1):
                words = lines[x+1].split('\t')
                Bi_up.append(float(words[1]))

        else:
            print('Running Bi up')

            Ti_up_full, Bi_up = model_up.get_VID(\
            Tbin = fiducial_Ti, T = model_up.T, Tmean = model_up.Tmean,\
            PT = model_up.PT,\
            PT_zero=model_up.PT_zero,\
            Nbin = Npt_i, minBi = min_Bi)

            if save_VID_flag:
                save_VID(filename_up,Ti_up_full,Bi_up)

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
    Fisher = np.zeros((len(der_Bi),len(der_Bi)))
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
                Fisher[i,j] +=  Fisher_Ti[t][i,j]

    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_VID/' + developer + '/' + 'fisher_tot' 
        save_fisher(filename, derive_pars,Fisher)

    phi_norm = np.zeros((len(Fisher),len(Fisher)))
    for i in range(len(Fisher)):
        #for j in range(len(Fisher)):
        #if i == j:
        phi_norm[i,i] = 10**-(np.log10(abs(Fisher[i,i])))/10
    print(Fisher)
                
    try:
        new_Fisher = np.linalg.multi_dot([phi_norm,Fisher])
        normalized = True
        print('Fisher normalized')
    except:
        new_Fisher = Fisher
        normalized = False

    inv_Fisher = np.linalg.inv(Fisher)

    sigma_tot = np.zeros(len(der_Bi))
    perc_error = np.zeros(len(der_Bi))
    for i in range(len(der_Bi)):
        sigma_tot[i] = np.sqrt(inv_Fisher[i,i])
        perc_error[i] = round(sigma_tot[i] / fid_par[i] * 100,2)

    print('Percentage errors: ' + str(perc_error) + '%')

    if plot_sigma_flag:

        print('Doing contour plots')
        names = get_fiducials(derive_pars, developer, z)[1]

        plot_dist = GaussianND(fid_par, inv_Fisher, names=names)			

        settings = gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = 40
        settings.legend_fontsize = 40
        settings.lab_fontsize = 40

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = True
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14
        g.settings.axis_tick_y_rotation = 65

        g.triangle_plot([plot_dist], names, filled=True,\
            legend_labels=[r'${\rm Ly\alpha}\, B_i\, {\rm forecasts}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'lw':2, 'color':kitsune}], markers={'x2':0})


        plt.savefig(save_fig_dir + 'ellipse.pdf')
        plt.show()

    for i in range(len(fid_par)):
        if type(fid_par[i]) is not np.float64 and type(fid_par[i]) is not float:
            fid_par[i] = fid_par[i].value

    print('Marginal errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fid_par * 100)

    return Fisher, sigma_tot


def error_dndL(z, developer):

    derive_pars = ['z0_Lyalpha','f0_Lyalpha','t_Lyalpha', 'astro_sig_scatter','astro_Lcut']

    Fisher = fisher_VID(z, developer, \
        derive_pars,\
    	save_VID_flag = True, \
        save_fisher_flag = True, \
	    plot_sigma_flag = False, \
        import_allowed = True)[0]

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))

    dz0_model_up = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid*(1+delta), zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    z0_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = dz0_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    dz0_model_low = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid*(1-delta), zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    z0_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = dz0_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    z0_low = update_VID(obs_pars(z),z0_low_data)[0]
    z0_up = update_VID(obs_pars(z),z0_up_data)[0]

    dn_dz0 = (z0_up.dndL.value-z0_low.dndL.value)/(2*delta*z0_Lyalpha_fid)

    df0_model_up = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid*(1+delta), SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    f0_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = df0_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    df0_model_low = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid*(1-delta), SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    f0_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = df0_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    f0_low = update_VID(obs_pars(z),f0_low_data)[0]
    f0_up = update_VID(obs_pars(z),f0_up_data)[0]

    dn_df0 = (f0_up.dndL.value-f0_low.dndL.value)/(2*delta*f0_Lyalpha_fid)

    t_model_up = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid*(1+delta), dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    t_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = t_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    t_model_low = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid*(1-delta),  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    t_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = t_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    t_low = update_VID(obs_pars(z),t_low_data)[0]
    t_up = update_VID(obs_pars(z),t_up_data)[0]

    dn_dt = (t_up.dndL.value-t_low.dndL.value)/(2*delta*t_Lyalpha_fid)

    sigma_model = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    sigma_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = sigma_model,
        sigma_scatter = astro_sig_scatter_fid*(1-delta),
        dndL_Lcut=astro_Lcut_fid)}

    sigma_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = sigma_model,
        sigma_scatter = astro_sig_scatter_fid*(1+delta),
        dndL_Lcut=astro_Lcut_fid)}

    sigma_low = update_VID(obs_pars(z),sigma_low_data)[0]
    sigma_up = update_VID(obs_pars(z),sigma_up_data)[0]

    dn_dsigma = (sigma_up.dndL.value-sigma_low.dndL.value)/(2*delta*astro_sig_scatter_fid)

    Lcut_model_up = dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid*(1+delta), SFR_file=SFR_file_fid)

    Lcut_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = Lcut_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid*(1+delta))}

    Lcut_model_low =dict(csi=csi_Lyalpha_fid,z0=z0_Lyalpha_fid, zeta=zeta_Lyalpha_fid, f0 = f0_Lyalpha_fid, SFR0 = SFR0_Lyalpha_fid, t = t_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid*(1-delta), SFR_file=SFR_file_fid)

    Lcut_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = Lcut_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid*(1-delta))}

    Lcut_low = update_VID(obs_pars(z),Lcut_low_data)[0]
    Lcut_up = update_VID(obs_pars(z),Lcut_up_data)[0]

    dn_dLcut = (Lcut_up.dndL.value-Lcut_low.dndL.value)/(2*delta*astro_Lcut_fid.value)


    sigma_dndl = np.zeros(len(dn_dz0))

    for i in range(len(sigma_dndl)):
        dz0_dn = dn_dz0[i]**-1 if dn_dz0[i] != 0. else 0.
        df0_dn = dn_df0[i]**-1 if dn_df0[i] != 0. else 0.
        dt_dn = dn_dt[i]**-1 if dn_dt[i] != 0. else 0.
        dLcut_dn = dn_dLcut[i]**-1 if dn_dLcut[i] != 0. else 0.
        dsigma_dn = dn_dsigma[i]**-1 if dn_dsigma[i] != 0. else 0.
        ders = np.array((dn_dz0[i], dn_df0[i], dn_dt[i], dn_dsigma[i], dn_dLcut[i]))
        new_cov = np.linalg.multi_dot([ders,np.linalg.inv(Fisher),ders.T])
        sigma_dndl[i] = np.sqrt(new_cov) if new_cov != 0. else 0.

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))

    model = update_VID(obs_pars(z), astrocosmo_dict(developer,z))[0]

    L = model.L.value
    fid_dn = model.dndL.value
    plt.figure(figsize=(15,9))
    plt.fill_between(L,fid_dn-sigma_dndl, fid_dn+sigma_dndl,color=kitsune,alpha=0.3)
    #plt.fill_between(L,fid_dn-2*sigma_dndl, fid_dn+2*sigma_dndl,color=kitsune,alpha=0.2)
    plt.loglog(L,fid_dn,color=kitsune)

    plt.xlabel(r'$L_{\rm Ly\alpha}\, [L_\odot]$',fontsize=40)
    plt.ylabel(r'$dn/dL_{\rm Ly\alpha}\, [{\rm Mpc}^{-3}L_\odot^{-1}]$',fontsize=40)
    plt.xlim(3e7,1e10)
    plt.ylim(1e-16,1e-9)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)

    plt.tight_layout()
    plt.savefig(save_fig_dir + 'error_dndl.pdf')
    plt.show()

    return sigma_dndl



def error_dndL_new(z, developer):

    derive_pars = ['A_lya', 'B_lya', 'D_lya', 'astro_sig_scatter','astro_Lcut']

    Fisher = fisher_VID(z, developer, \
        derive_pars,\
    	save_VID_flag = True, \
        save_fisher_flag = True, \
	    plot_sigma_flag = False, \
        import_allowed = True)[0]

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))

    dA_model_up = dict(A_lya=A_Lyalpha_fid*(1+delta), B_lya=B_Lyalpha_fid,D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    A_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = dA_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    dA_model_low = dict(A_lya=A_Lyalpha_fid*(1-delta), B_lya=B_Lyalpha_fid,D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    A_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = dA_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    A_low = update_VID(obs_pars(z),A_low_data)[0]
    A_up = update_VID(obs_pars(z),A_up_data)[0]

    dn_dA = (A_up.dndL.value-A_low.dndL.value)/(2*delta*A_Lyalpha_fid)

    dB_model_up = dict(A_lya=A_Lyalpha_fid, B_lya=B_Lyalpha_fid*(1+delta),D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    B_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = dB_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    dB_model_low = dict(A_lya=A_Lyalpha_fid, B_lya=B_Lyalpha_fid*(1-delta),D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    B_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = dB_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    B_low = update_VID(obs_pars(z),B_low_data)[0]
    B_up = update_VID(obs_pars(z),B_up_data)[0]

    dn_dB = (A_up.dndL.value-B_low.dndL.value)/(2*delta*B_Lyalpha_fid)

    dD_model_up = dict(A_lya=A_Lyalpha_fid, B_lya=B_Lyalpha_fid,D_lya=D_Lyalpha_fid*(1+delta),  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    D_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = dD_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    dD_model_low = dict(A_lya=A_Lyalpha_fid, B_lya=B_Lyalpha_fid,D_lya=D_Lyalpha_fid*(1-delta),  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    D_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = dD_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

    D_low = update_VID(obs_pars(z),D_low_data)[0]
    D_up = update_VID(obs_pars(z),D_up_data)[0]

    dn_dD = (D_up.dndL.value-D_low.dndL.value)/(2*delta*D_Lyalpha_fid)

    sigma_model = dict(A_lya=A_Lyalpha_fid,B_lya=B_Lyalpha_fid, D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

    sigma_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = sigma_model,
        sigma_scatter = astro_sig_scatter_fid*(1-delta),
        dndL_Lcut=astro_Lcut_fid)}

    sigma_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = sigma_model,
        sigma_scatter = astro_sig_scatter_fid*(1+delta),
        dndL_Lcut=astro_Lcut_fid)}

    sigma_low = update_VID(obs_pars(z),sigma_low_data)[0]
    sigma_up = update_VID(obs_pars(z),sigma_up_data)[0]

    dn_dsigma = (sigma_up.dndL.value-sigma_low.dndL.value)/(2*delta*astro_sig_scatter_fid)

    Lcut_model_up =dict(A_lya=A_Lyalpha_fid,B_lya=B_Lyalpha_fid, D_lya=D_Lyalpha_fid,
    dndL_Lcut=astro_Lcut_fid*(1+delta), SFR_file=SFR_file_fid)

    Lcut_up_data = {**astrocosmo_dict(developer,z), **dict(model_par = Lcut_model_up,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid*(1+delta))}

    Lcut_model_low =dict(A_lya=A_Lyalpha_fid,B_lya=B_Lyalpha_fid, D_lya=D_Lyalpha_fid,  dndL_Lcut=astro_Lcut_fid*(1-delta), SFR_file=SFR_file_fid)

    Lcut_low_data = {**astrocosmo_dict(developer,z), **dict(model_par = Lcut_model_low,
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid*(1-delta))}

    Lcut_low = update_VID(obs_pars(z),Lcut_low_data)[0]
    Lcut_up = update_VID(obs_pars(z),Lcut_up_data)[0]

    dn_dLcut = (Lcut_up.dndL.value-Lcut_low.dndL.value)/(2*delta*astro_Lcut_fid.value)


    sigma_dndl = np.zeros(len(dn_dA))

    for i in range(len(sigma_dndl)):

        ders = np.array((dn_dA[i], dn_dB[i], dn_dD[i], dn_dsigma[i], dn_dLcut[i]))
        #new_cov = np.linalg.multi_dot([ders,np.linalg.inv(Fisher),ders.T])
        new_cov = np.linalg.multi_dot([ders,np.linalg.inv(Fisher),ders.T])
        sigma_dndl[i] = np.sqrt(new_cov) if new_cov != 0. else 0.

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))

    model = update_VID(obs_pars(z), astrocosmo_dict(developer,z))[0]

    L = model.L.value
    fid_dn = model.dndL.value
    plt.figure(figsize=(15,9))
    plt.fill_between(L,fid_dn-sigma_dndl, fid_dn+sigma_dndl,color=kitsune,alpha=0.3)
    plt.fill_between(L,fid_dn-2*sigma_dndl, fid_dn+2*sigma_dndl,color=kitsune,alpha=0.2)
    plt.loglog(L,fid_dn,color=kitsune)

    plt.xlabel(r'$L_{\rm Ly\alpha}\, [L_\odot]$',fontsize=40)
    plt.ylabel(r'$dn/dL_{\rm Ly\alpha}\, [{\rm Mpc}^{-3}L_\odot^{-1}]$',fontsize=40)
    plt.xlim(3e7,1e10)
    plt.ylim(1e-16,1e-9)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)

    plt.tight_layout()
    plt.savefig(save_fig_dir + 'error_dndl.pdf')
    plt.show()

    return sigma_dndl








def dndl_ext():

    ders = ['z0_Lyalpha','f0_Lyalpha','t_Lyalpha','astro_sig_scatter','astro_Lcut']

    Fisher_us = [[3286119.1682657185,-51043502.01681407,2777673.5541573227,-53138637.822572134,-0.0022699042629170055],
[-51043502.01681407,801763566.3865196,-43643334.68054139,848177036.2728131,0.03480224170042015],
[2777673.5541573227,-43643334.68054139,2375736.1508579375,-46186961.64556567,-0.001877698650681227],
[-53138637.822572134,848177036.2728131,-46186961.64556567,917772889.925487,0.03695665512524026],
[-0.0022699042629170055,0.03480224170042015,-0.001877698650681227,0.03695665512524026,1.0186612596869506e-11]
]


    Fisher_sp = [[0.08325385965355703,-12.646840021251373,0.4667457453050122,-0.11047524495238079,7.396742728753806e-11],
[-12.646840021251373,1934.643579658818,-72.63571358264359,25.072562998702963,-3.5456406979343895e-08],
[0.4667457453050122,-72.63571358264359,2.842049757992649,-1.673488415447301,4.018971696073125e-09],
[-0.11047524495238079,25.072562998702963,-1.673488415447301,5.2809702394821,-1.2968746470492966e-08],
[7.396742728753806e-11,-3.5456406979343895e-08,4.018971696073125e-09,-1.2968746470492966e-08,1.390349430014801e-16]
]

    Fisher = [Fisher_us,Fisher_sp]
    developer = 'CDS'

    z = [1.135,5.74]
    z0 = [z0_Lyalpha_fid,3,125]
    f0 = [f0_Lyalpha_fid,0.18]
    t = [t_Lyalpha_fid,0.875]
    color=[kitsune,'k']
    style=['-','--']
 
    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))
    plt.figure(figsize=(15,9))

    for j in range(len(z)):
        dz0_model_up = dict(csi=csi_Lyalpha_fid,z0=z0[j]*(1+delta), zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        z0_up_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = dz0_model_up,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        dz0_model_low = dict(csi=csi_Lyalpha_fid,z0=z0[j]*(1-delta), zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        z0_low_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = dz0_model_low,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        z0_low = update_VID(obs_pars(z[j]),z0_low_data)[0]
        z0_up = update_VID(obs_pars(z[j]),z0_up_data)[0]

        dn_dz0 = (z0_up.dndL.value-z0_low.dndL.value)/(2*delta*z0[j])

        df0_model_up = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j]*(1+delta), SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        f0_up_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = df0_model_up,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        df0_model_low = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j]*(1-delta), SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        f0_low_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = df0_model_low,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        f0_low = update_VID(obs_pars(z[j]),f0_low_data)[0]
        f0_up = update_VID(obs_pars(z[j]),f0_up_data)[0]

        dn_df0 = (f0_up.dndL.value-f0_low.dndL.value)/(2*delta*f0[j])

        t_model_up = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j]*(1+delta), dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        t_up_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = t_model_up,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        t_model_low = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j]*(1-delta),  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        t_low_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = t_model_low,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid)}

        t_low = update_VID(obs_pars(z[j]),t_low_data)[0]
        t_up = update_VID(obs_pars(z[j]),t_up_data)[0]

        dn_dt = (t_up.dndL.value-t_low.dndL.value)/(2*delta*t[j])

        sigma_model = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)

        sigma_low_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = sigma_model,
            sigma_scatter = astro_sig_scatter_fid*(1-delta),
            dndL_Lcut=astro_Lcut_fid)}

        sigma_up_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = sigma_model,
            sigma_scatter = astro_sig_scatter_fid*(1+delta),
            dndL_Lcut=astro_Lcut_fid)}

        sigma_low = update_VID(obs_pars(z[j]),sigma_low_data)[0]
        sigma_up = update_VID(obs_pars(z[j]),sigma_up_data)[0]

        dn_dsigma = (sigma_up.dndL.value-sigma_low.dndL.value)/(2*delta*astro_sig_scatter_fid)

        Lcut_model_up = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid*(1+delta), SFR_file=SFR_file_fid)

        Lcut_up_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = Lcut_model_up,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid*(1+delta))}

        Lcut_model_low =dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid*(1-delta), SFR_file=SFR_file_fid)

        Lcut_low_data = {**astrocosmo_dict(developer,z[j]), **dict(model_par = Lcut_model_low,
            sigma_scatter = astro_sig_scatter_fid,
            dndL_Lcut=astro_Lcut_fid*(1-delta))}

        Lcut_low = update_VID(obs_pars(z[j]),Lcut_low_data)[0]
        Lcut_up = update_VID(obs_pars(z[j]),Lcut_up_data)[0]

        dn_dLcut = (Lcut_up.dndL.value-Lcut_low.dndL.value)/(2*delta*astro_Lcut_fid.value)

        sigma_dndl = np.zeros(len(dn_dz0))

        for i in range(len(sigma_dndl)):
            ders = np.array((dn_dz0[i], dn_df0[i], dn_dt[i], dn_dsigma[i], dn_dLcut[i]))
            new_cov = np.linalg.multi_dot([ders,np.linalg.inv(Fisher[j]),ders.T])
            sigma_dndl[i] = np.sqrt(new_cov) if new_cov != 0. else 0.

        fid_mod = {**astrocosmo_dict('CDS',z[j]), **dict(model_par = dict(csi=csi_Lyalpha_fid,z0=z0[j], zeta=zeta_Lyalpha_fid, f0 = f0[j], SFR0 = SFR0_Lyalpha_fid, t = t[j],  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid),
        sigma_scatter = astro_sig_scatter_fid,
        dndL_Lcut=astro_Lcut_fid)}

        fid = update_VID(obs_pars(z[j]),fid_mod)[0]

        L = fid.L.value
        fid_dn = fid.dndL.value
        plt.fill_between(L,fid_dn-sigma_dndl, fid_dn+sigma_dndl,color=color[j],alpha=0.3)
        plt.loglog(fid.L,fid.dndL,color=color[j],linestyle=style[j])

    plt.xlabel(r'$L_{\rm Ly\alpha}\, [L_\odot]$',fontsize=40)
    plt.ylabel(r'$dn/dL_{\rm Ly\alpha}\, [{\rm Mpc}^{-3}L_\odot^{-1}]$',fontsize=40)
    plt.xlim(3e7,1e10)
    plt.ylim(1e-16,1e-9)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)

    plt.tight_layout()
    #plt.savefig(save_fig_dir + 'error_dndl.pdf')
    plt.show()
    return

#    plt.errorbar(z[0],[f_esc(z[0],z0_Lyalpha_fid,f0_Lyalpha_fid,t_Lyalpha_fid)],yerr=[sigma_dndl_us],color=kitsune)

#    plt.errorbar(z[1],[f_esc(z[1],3.125,0.18,0.875)],yerr=[sigma_dndl_spherex],color=kitsune)

    plt.xlabel(r'$z$')
    plt.ylabel(r'$f_{\rm esc}$')

    plt.show()

    return sigma_dndl_us, sigma_dndl_spherex