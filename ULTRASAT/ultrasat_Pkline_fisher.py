# SL: last update 01/23/2023

from LIM_b7 import *
from LIM_b7.ultrasat_fiducial_pars import astrocosmo_dict

save_fig_dir = './results_ultrasat/new/'#spherex/' 

if hmf_fid == 'NG_Riotto':
    save_fig_dir += 'fNL' + str(fNL_fid) + '/'

create_dir(save_fig_dir)
create_dir(save_fig_dir + 'Pkline/')
create_dir(save_fig_dir + 'fisher_Pkline/')

direct_id = lambda developer, model_id, z: save_fig_dir + 'Pkline/' + developer + '/' + model_id + '/z' + str(z)

direct_fisher_id = lambda developer, model_id, z: save_fig_dir + 'fisher_Pkline/' + developer + '/' + model_id + '/z' + str(z)


# compute power spectrum with redshift space distortions
def Pk_AP(model, z, alpha_par = 1., alpha_perp = 1.):

    Tfs8, Tbs8, s8 = get_fiducials_Pkline(['Tfs8','Tbs8','s8'], model, z)[0]

    Pshot = model.L2mean*model.CLT**2. 
    sNL = model.sigma_NL

    F = alpha_par/alpha_perp
    prefac = 1./alpha_perp**2./alpha_par # volume modification

    kprime = np.zeros((len(model.mu),len(model.k)))
    mu_prime = model.mui_grid/F/np.sqrt(1.+model.mui_grid**2.*(F**-2-1.))
    for imu in range(model.nmu):
        kprime[imu,:] = (model.k /alpha_perp)*np.sqrt(1.+model.mu[imu]**2*(F**-2-1.)) 

    #Obtain the corresponding P_m for the "real" k
    Pmprime = log_interp1d(model.k,model.Pm[0,:])(kprime)*model.Pm.unit

    #Apply RSD for the "real" mu
    kaiser = (Tbs8 + Tfs8*mu_prime**2)**2
    loren = 1./(1.+0.5*(kprime*mu_prime*sNL.value)**2)

    #Units change depend on the survey. Autoconsistent definition in terms of Mpc/h
    Pkclus = (kaiser*Pmprime*loren**2/s8**2) 

    Pkclus = prefac*((Pkclus + Pshot)*(model.Wkmin*model.Wkmax)) 

    return Pkclus


def analytical_derivatives(derive_pars, model,z,Pkclus,get_multipole):

    ders = []

    Tfs8, Tbs8, fs8, s8 = get_fiducials_Pkline(['Tfs8','Tbs8','fs8','s8'], model, z)[0]

    for i in range(len(derive_pars)):
        if derive_pars[i] == 'Tfs8':
            analitic_der =  (2.*Pkclus*(model.mui_grid**2/ (Tbs8+Tfs8*model.mui_grid**2)))	

            der_marg_rsd = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
            (2*get_multipole+1)/2*np.trapz(analitic_der*legendre(get_multipole)(model.mui_grid),model.mu,axis=0)

            ders.append(der_marg_rsd.value)

        if derive_pars[i] == 'Tbs8':
 
            analitic_der = (2.*Pkclus/(Tbs8+Tfs8*model.mui_grid**2))

            der_marg_rsd = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
            (2*get_multipole+1)/2*np.trapz(analitic_der*legendre(get_multipole)(model.mui_grid),model.mu,axis=0)
 
            ders.append(der_marg_rsd.value)
 
        if derive_pars[i] == 'Pshot':
 
            ders.append(np.ones(model.k.shape))

        if derive_pars[i] == 'sNL':
 
            sNL = model.sigma_NL

            analitic_der =  (-2.*Pkclus*sNL*(model.ki_grid*model.mui_grid)**2./\
                (1.+0.5*(sNL*model.ki_grid*model.mui_grid)**2.))
                
            der_marg_rsd =  0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
            (2*get_multipole+1)/2*np.trapz(analitic_der*legendre(get_multipole)(model.mui_grid),model.mu,axis=0)
            
            ders.append(der_marg_rsd.value) 


        if derive_pars[i] == 'alpha_par':

            Pkclus_up = Pk_AP(model,z, 1.+delta, 1.)
            Pkclus_low = Pk_AP(model,z, 1.-delta, 1.)

            analitic_der = (Pkclus_up - Pkclus_low).value / (2*delta)
            derivative = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                        (2*get_multipole+1)/2*np.trapz(analitic_der*legendre(get_multipole)(model.mui_grid),model.mu,axis=0) 

            ders.append(derivative) 

        if derive_pars[i] == 'alpha_perp':

            Pkclus_up = Pk_AP(model,z, 1., 1.+delta)
            Pkclus_low = Pk_AP(model,z, 1., 1.-delta)

            analitic_der = (Pkclus_up - Pkclus_low).value / (2*delta)
            derivative =  0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                        (2*get_multipole+1)/2*np.trapz(analitic_der*legendre(get_multipole)(model.mui_grid),model.mu,axis=0)

            ders.append(derivative) 


    return ders 


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
        use_f_axion_low = f_axion_fid*(1. - delta)
        use_f_axion_up = f_axion_fid*(1. + delta)
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

    if hmf_fid == 'NG_Riotto' and derive_par == 'fNL':
        use_fNL_low = fNL_fid*(1. - 10*delta)
        use_fNL_up = fNL_fid*(1. + 10*delta)
    else:
        use_fNL_low = fNL_fid
        use_fNL_up = fNL_fid


#    low_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_low, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_low, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_low,  dndL_Lcut=use_astro_Lcut_low, SFR_file=SFR_file_fid)}
    low_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_low,B_lya=use_B_Lyalpha_low, D_lya=use_D_Lyalpha_low, dndL_Lcut=use_astro_Lcut_low, SFR_file=SFR_file_fid)}
                     
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

#    up_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_up, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_up, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_up,  dndL_Lcut=use_astro_Lcut_up, SFR_file=SFR_file_fid)}
    up_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_up,B_lya=use_B_Lyalpha_up, D_lya=use_D_Lyalpha_up, dndL_Lcut=use_astro_Lcut_up, SFR_file=SFR_file_fid)}

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

    low = update_Pkline(obs_pars(z),low_data)[0]
    up = update_Pkline(obs_pars(z),up_data)[0]

    return low, up


def fisher_Pkline(z, developer, \
        derive_model = ['a1','a2','a3','a4','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','a_tinker'],\
        derive_Pkline = ['Tfs8','Tbs8','Pshot','sNL','alpha_par','alpha_perp'],	multipole = [0,2,4],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):

    if 'fNL' in derive_Pkline:
        print('!!! Derivative wrt fNL is computed numerically !!!')

    create_dir(save_fig_dir + 'Pkline/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_Pkline/' + developer + '/')

    create_dir(save_fig_dir + 'Pkline/' + developer + '/' )
    create_dir(save_fig_dir + 'fisher_Pkline/' + developer + '/')

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT(z_val))
    #obs_pars = lambda z_val: deepcopy(obs_params_SPHEREx(z_val))
    
    model = update_Pkline(obs_pars(z), astrocosmo_dict(developer,z))[0]

    derive_pars = derive_model + derive_Pkline
    Fisher = np.zeros((len(derive_pars),len(derive_pars)))


    direct = save_fig_dir + 'Pkline/' + developer  + '/z' + str(z)
    if developer == 'axions':
        direct += '_' + str(f_axion_fid) 

    direct += '/'
    direct_fisher = save_fig_dir + 'fisher_Pkline/' + developer  + '/z' + str(z)
    if developer == 'axions':
        direct_fisher += '_' + str(f_axion_fid) 

    create_dir(direct)
    create_dir(direct_fisher)

    print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')


    fid_par, name_fid = get_fiducials(derive_model, developer, z)

    vecpar, name_vec = get_fiducials_Pkline(derive_Pkline, model, z)

    fiducials = fid_par + vecpar
    names = name_fid + name_vec

    fiducial_k = model.k.value 

    fiducial_k = np.asarray(fiducial_k) 
    Pkclus_fid = Pk_AP(model,z,1.,1.)

    #integrand_cov = (model.Pk+model.Pnoise)/(model.Nmodes**0.5) # Pnoise includes Vvox --> N_z_shells
    integrand_cov = (Pkclus_fid+model.Pnoise)/(model.Nmodes**0.5) # Pnoise includes Vvox --> N_z_shells
    cov_02 = (5./2.*np.trapz(integrand_cov**2*legendre(2)(model.mui_grid)**2,model.mu,axis=0)).value
    if 4 in multipole:
        cov_04 = (9./2.*np.trapz(integrand_cov**2*legendre(4)(model.mui_grid)**2,model.mu,axis=0)).value
        cov_24 = (45./2.*np.trapz(integrand_cov**2*legendre(2)(model.mui_grid)*legendre(4)(model.mui_grid),model.mu,axis=0)).value

        cov = [[[],cov_02,cov_04],
            [cov_02,[],cov_24],
            [cov_04,cov_24,[]]]
        der_Pk_multi = [[],[],[]]
    else:
        cov = [[[],cov_02],
            [cov_02,[]]]
        der_Pk_multi = [[],[]]

    for l in range(len(multipole)):

        print('--------- Doing L = ' + str(multipole[l]) + ' ---------')
        print('Running Pk fiducial')
        fiducial_Pk = 0.5*np.trapz(Pkclus_fid, model.mu, axis=0) if multipole[l] == 0 else (2*multipole[l]+1)/2.*np.trapz(Pkclus_fid*legendre(multipole[l])(model.mui_grid),model.mu,axis=0) 

        cov[l][l] = np.zeros(len(fiducial_Pk))
                
        sigma2 = 0.5*np.trapz(integrand_cov**2,model.mu,axis=0) if multipole[l] == 0 else\
                (2*multipole[l]+1)*(2*multipole[l]+1)/2.*np.trapz(integrand_cov**2*legendre(multipole[l])(model.mui_grid)*\
                legendre(multipole[l])(model.mui_grid),model.mu,axis=0) 

        for kval in range(len(sigma2)):
            cov[l][l][kval] = sigma2[kval].value


        all_k_low = []
        all_k_up = []
        all_Pk_low = []
        all_Pk_up = []
        for i in range(len(derive_model)):
            filename_low = direct + 'Pkline_' + str(multipole[l])+ '_low_' + str(derive_model[i])			
            filename_up = direct + 'Pkline_' + str(multipole[l]) + 'up_' + str(derive_model[i])

            if os.path.exists(filename_low) and import_allowed:
                k_low, Pk_low = import_Pkline(filename_low, '\nDoing derivative wrt: ' + str(derive_model[i]) + '\nImport Pk low')
            else:
                model_low, model_up = set_der_par(derive_model[i],z,developer,obs_pars)

                Pkclus_low = Pk_AP(model_low,z,1.,1.)
                Pkclus_up = Pk_AP(model_up,z,1.,1.)

                Pk_low = 0.5*np.trapz(Pkclus_low , model_low.mu, axis=0) if multipole[l] == 0 else\
                (2*multipole[l]+1)/2.*np.trapz(Pkclus_low*legendre(multipole[l])(model_low.mui_grid),model_low.mu,axis=0)

                k_low = model_low.k

                if save_Pkline_flag:
                    np.savetxt(filename_low,(k_low,Pk_low))

            all_k_low.append(np.asarray(k_low)) 
            all_Pk_low.append(np.asarray(Pk_low))

            if os.path.exists(filename_up) and import_allowed:
                k_up, Pk_up = import_Pkline(filename_up, 'Import Pk up')
            else:
                Pk_up = 0.5*np.trapz(Pkclus_up, model_up.mu, axis=0) if multipole[l] == 0 else\
                        (2*multipole[l]+1)/2.*np.trapz(Pkclus_up*legendre(multipole[l])(model_up.mui_grid),model_up.mu,axis=0) 
                k_up = model_up.k

                if save_Pkline_flag:
                    np.savetxt(filename_up, (k_up,Pk_up))

            all_k_up.append(np.asarray(k_up)) 
            all_Pk_up.append(np.asarray(Pk_up))

        pkprop_ders = analytical_derivatives(derive_Pkline, model,z,Pkclus_fid,multipole[l])

        der_k = []
        der_Pk = []
        for i in range(len(derive_pars)):

            if derive_pars[i] in derive_model:
                id_der = derive_model.index(derive_pars[i])
                der_k.append([])
                der_Pk.append([])
                for t in range(len(fiducial_k)):
                    der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / (2*delta*fid_par[id_der]))

            elif derive_pars[i] in derive_Pkline:
                id_der = derive_Pkline.index(derive_pars[i])
                all_k_up.append(model.k)
                der_Pk.append(pkprop_ders[id_der])

            der_Pk_multi[l] = der_Pk

    inv_cov = []
    multi_der_Pk_ordered = []
    if 4 in multipole:
        for kval in range(len(cov[0][0])):
            use_cov = np.asarray([[cov[0][0][kval],cov[0][1][kval],cov[0][2][kval]],[cov[1][0][kval],cov[1][1][kval],cov[1][2][kval]],[cov[2][0][kval],cov[2][1][kval],cov[2][2][kval]]])
            inv_cov.append(np.linalg.inv(use_cov))
        for a in range(len(derive_pars)):
            multi_der_Pk_ordered.append([])	
            for kval in range(len(der_Pk_multi[0][a])):
                multi_der_Pk_ordered[a].append([der_Pk_multi[0][a][kval],der_Pk_multi[1][a][kval],der_Pk_multi[2][a][kval]])
    else:
        for kval in range(len(cov[0][0])):
            use_cov = np.asarray([[cov[0][0][kval],cov[0][1][kval]],[cov[1][0][kval],cov[1][1][kval]]])
            inv_cov.append(np.linalg.inv(use_cov))
        for a in range(len(derive_pars)):
            multi_der_Pk_ordered.append([])	
            for kval in range(len(der_Pk_multi[0][a])):
                multi_der_Pk_ordered[a].append([der_Pk_multi[0][a][kval],der_Pk_multi[1][a][kval]])

    inv_cov = np.asarray(inv_cov)
    multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)

    #return inv_cov, multi_der_Pk_ordered
    Fisher_k = []
    Fisher = np.zeros((len(der_Pk),len(der_Pk)))
    for t in range(len(fiducial_k)):
        Fisher_k.append(np.zeros((len(der_Pk),len(der_Pk))))		
        for i in range(len(der_Pk)):
            for j in range(len(der_Pk)):
                Fisher_k[t][i,j] = np.linalg.multi_dot(\
                        [multi_der_Pk_ordered[i][t],inv_cov[t],\
                        multi_der_Pk_ordered[j][t]])
                
    Fisher_k = np.asarray(Fisher_k)
    for i in range(len(der_Pk)):
        for j in range(len(der_Pk)):
            for t in range(len(fiducial_k)):
                Fisher[i,j] +=  Fisher_k[t][i,j]
        
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkline/fisher_tot' 
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

    print(new_Fisher)
    start_inv_Fisher = np.linalg.inv(new_Fisher)
    if normalized:
        inv_Fisher = np.linalg.multi_dot([phi_norm,start_inv_Fisher])
    else:
        inv_Fisher = start_inv_Fisher

    sigma_tot = np.zeros(len(derive_pars))
    for i in range(len(derive_pars)):
        try:
            sigma_tot[i] = np.sqrt(inv_Fisher[i,i]).value
        except:
            sigma_tot[i] = np.sqrt(inv_Fisher[i,i])

    print(sigma_tot)
    for i in range(len(fiducials)):
        if type(fiducials[i]) is not np.float64 and type(fiducials[i]) is not float:
            fiducials[i] = fiducials[i].value


    if plot_sigma_flag:

        print('Doing contour plots')

        plot_dist = GaussianND(fiducials, inv_Fisher, names=names)			

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
            legend_labels=[r'${\rm Ly\alpha}\, Pk\, {\rm forecasts}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'lw':2, 'color':kitsune}], markers={'x2':0})


        plt.savefig(save_fig_dir + 'ellipse_pk.pdf')
        plt.show()
        '''
        plot_dist = GaussianND(fiducials, inv_Fisher, names=names)			

        settings = gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = 20
        settings.legend_fontsize = 25
        settings.lab_fontsize = 25

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = False
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14

        g.triangle_plot([plot_dist], names, filled=True,\
            legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers={'x2':0})

        plt.savefig(save_fig_dir + 'ellipse_pk.pdf')
        plt.show()
        '''
        new_names = [r'$10^5D_A(z)(1+z)^2$',r'$H(z)/(1+z)^{3/2}$']

        inv_Fisher_DAH = [[inv_Fisher[-1][-1],inv_Fisher[-1][-2]],[inv_Fisher[-2][-1],inv_Fisher[-2][-2]]]

        plot_dist = GaussianND([7.9e3/1e5,41.8], inv_Fisher_DAH, names=new_names)			
        # plot_dist = GaussianND([5.6e4,37.9], inv_Fisher_DAH, names=new_names)			

        print('\nError on DA/rs = ' + str(np.sqrt(inv_Fisher[-1][-1])*11.61 * 100))
        print('Error on Hrs = ' + str(np.sqrt(inv_Fisher[-2][-2]) * 100) + '\n')
        # print('\nError on DA/rs = ' + str(np.sqrt(inv_Fisher[-1][-1])*8.22*100))
        # print('Error on Hrs = ' + str(np.sqrt(inv_Fisher[-2][-2]) *100)+ '\n')

        settings = gp.GetDistPlotSettings()

        settings.fig_width_inch = 7
        settings.norm_prob_label = False
        settings.axes_fontsize = 20
        settings.legend_fontsize = 25
        settings.lab_fontsize = 25

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = False
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14

        g.plot_2d([plot_dist], new_names, filled=True,\
            #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers=[{'x2':0}])

        plt.tight_layout()
        plt.savefig(save_fig_dir + 'ellipse_pk_DAH.pdf')
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot


def plot_combo():


    print('Doing contour plots')
    new_names = [r'$D_A(z)(1+z)^2$',r'$H(z)/(1+z)^{3/2}$']


    fisher_us = [[107908034.31026058,-1579893712.0981967,99170374.71607204,-42.0550626566643,-33672307.140302755,-2.850759309522287e+16,26628.715457146416,170495899.45658666,327036276.5957075],
[-1579893712.0981967,23140686519.10109,-1452693624.95489,582.0309567247424,493000077.0916768,4.204429005104437e+17,-391708.5258562882,-2496249489.3734527,-4795707424.797494],
[99170374.71607207,-1452693624.9548898,91197467.9261914,-36.00215476901637,-30945754.5648389,-2.643950341789131e+16,24616.541321071563,156690282.2605156,301146516.9510186],
[-42.0550626566643,582.0309567247423,-36.00215476901637,0.2124470259550733,12.95309062837356,8895691854.577436,-0.00854598982008228,-69.83909592479358,-122.01130391895383],
[-33672307.140302755,493000077.09167683,-30945754.564838894,12.953090628373559,10507320.331474854,8895691854191547.0,-8309.394658108082,-53202618.65960785,-102050477.44797397],
[-2.8507593095222856e+16,4.2044290051044365e+17,-2.6439503417891308e+16,8895691854.577436,8895691854191545.0,1.3105571238723815e+28,-12113002924889.908,-4.504243439572226e+16,-8.910753197486542e+16],
[26628.715457146413,-391708.5258562882,24616.541321071563,-0.008545989820082282,-8309.394658108082,-12113002924889.912,7.892304010563655,42073.70114779437,82296.69459526648],
[170495899.45658666,-2496249489.373453,156690282.26051566,-69.8390959247936,-53202618.65960785,-4.5042434395722264e+16,42073.70114779437,269385480.5509899,516721131.91791344],
[327036276.5957075,-4795707424.797495,301146516.95101863,-122.01130391895383,-102050477.44797395,-8.910753197486544e+16,82296.69459526648,516721131.9179133,997349826.8810705]
]
    
    fisher_spherex= [[37101.42200735512,-15370552.272697376,1921946.0892772682,-95.87903074312617,-3172734.2456136155,-9800781.19500189,3297.804672466154,1544590.6843691724,3036537.3455087463],
[-15370552.272697376,6371860945.105638,-797006943.2173152,27730.081230509044,1314418247.4761074,4078727310.87035,-1364692.7776885305,-639917919.9420037,-1260023839.049314],
[1921946.0892772682,-797006943.2173152,99708589.10611138,-2686.004846949351,-164356094.3904687,-511201497.8679478,170542.5652966904,80016987.64793757,157686222.1328999],
[-95.8790307431262,27730.081230509044,-2686.0048469493504,804.2468031068785,8002.986975357364,3083.53421497587,-11.371702697855742,-3912.4450396647644,-5114.76899367081],
[-3172734.245613616,1314418247.4761076,-164356094.3904687,8002.986975357362,271316942.14551014,838120852.6677353,-282011.61749735166,-132085947.40978815,-259670660.36670962],
[-9800781.195001893,4078727310.8703513,-511201497.8679479,3083.53421497587,838120852.6677351,4831257378.327812,-798621.4213982024,-408100783.0354253,-812067846.6243827],
[3297.804672466154,-1364692.7776885307,170542.56529669036,-11.37170269785574,-282011.6174973516,-798621.4213982022,302.91391593723273,137286.30005696716,269065.47471036],
[1544590.684369173,-639917919.9420037,80016987.64793757,-3912.4450396647644,-132085947.40978816,-408100783.0354253,137286.3000569672,64303829.62711403,126424610.00065374],
[3036537.3455087463,-1260023839.049314,157686222.13289988,-5114.76899367081,-259670660.36670962,-812067846.6243826,269065.47471036005,126424610.00065374,249562786.40800065]
]

    z = np.linspace(0,7)
    detector_params = lambda z: obs_params_ULTRASAT(1.135)
    model_data = dict(\
        developer = 'CDS',               
        CDS_alpha_1 = CDS_alpha_1_fid, 
            CDS_alpha_2 = CDS_alpha_2_fid, 
        CDS_alpha_3 = CDS_alpha_3_fid, 
            CDS_alpha_4 = CDS_alpha_4_fid, 
        k_12 = k_12_fid, 
            k_23 = k_23_fid, 
            k_34 = k_34_fid )

    data = lambda developer, z:{**astrocosmo_dict(developer,z), **model_data}
    mod = update_Pkline(\
	        detector_params(1.135),\
		    data(model_data['developer'],1.135))[0]

    cosmo = mod.cosmo

    Da = np.zeros(len(z))
    H = np.zeros(len(z))
    for i in range(len(z)):
        Da[i] = cosmo.angular_diameter_distance(z[i])*(1+z[i])**2/1e5
        H[i] = cosmo.hubble_parameter(z[i])/(1+z[i])**(3/2)
    Da_z1 =  (cosmo.angular_diameter_distance(1.135)*(1+1.135)**2)/1e5
    Da_z5 =  (cosmo.angular_diameter_distance(5.74)*(1+5.74)**2)/1e5
    H_z1 = cosmo.hubble_parameter(1.135)/(1+1.135)**(3/2)
    H_z5 = cosmo.hubble_parameter(5.74)/(1+5.74)**(3/2)

    inv_Fisher_us = np.linalg.inv(fisher_us)
    inv_Fisher_spherex = np.linalg.inv(fisher_spherex)

    inv_Fisher_DAH_ultrasat = [[inv_Fisher_us[-1][-1],inv_Fisher_us[-1][-2]],[inv_Fisher_us[-2][-1],inv_Fisher_us[-2][-2]]]

    inv_Fisher_DAH_spherex = [[inv_Fisher_spherex[-1][-1],inv_Fisher_spherex[-1][-2]],[inv_Fisher_spherex[-2][-1],inv_Fisher_spherex[-2][-2]]]

    plot_dist_ultrasat = GaussianND([Da_z1-Da_z1,H_z1-H_z1], inv_Fisher_DAH_ultrasat, names=new_names)			
    plot_dist = GaussianND([Da_z5-Da_z5,H_z5-H_z5], inv_Fisher_DAH_spherex, names=new_names)			

    settings = gp.GetDistPlotSettings()

    settings.fig_width_inch = 7
    settings.norm_prob_label = False
    settings.axes_fontsize = 20
    settings.legend_fontsize = 25
    settings.lab_fontsize = 25

    g = gp.get_subplot_plotter(settings=settings)

    g.settings.figure_legend_frame = False
    g.settings.alpha_filled_add=0.4
    g.settings.title_limit_fontsize = 14

   # gp.plt.plot(Da,H,'k')

    g.plot_2d([plot_dist_ultrasat,plot_dist,], new_names, filled=True,\
        #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
        contour_colors=[kitsune,'k',], line_args=[{'ls':'-', 'color':kitsune},{'ls':'-', 'color':'k'},],
        markers=[{'x2':0},{'x2':0}])
    #gp.plt.plot(Da,H,'k:')
    #gp.plt.xlim(0.05,0.6)
    #gp.plt.ylim(36,45)

    plt.tight_layout()
    plt.savefig(save_fig_dir + 'ellipse_pk_DAH_combo.pdf')
    plt.show()

    return 


import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def plot_constraints():

    z = np.linspace(0,7)
    detector_params = lambda z: obs_params_ULTRASAT(1.135)
    model_data = dict(\
        developer = 'CDS',               
        CDS_alpha_1 = CDS_alpha_1_fid, 
            CDS_alpha_2 = CDS_alpha_2_fid, 
        CDS_alpha_3 = CDS_alpha_3_fid, 
            CDS_alpha_4 = CDS_alpha_4_fid, 
        k_12 = k_12_fid, 
            k_23 = k_23_fid, 
            k_34 = k_34_fid )

    data = lambda developer, z:{**astrocosmo_dict(developer,z), **model_data}
    mod = update_Pkline(\
	        detector_params(1.135),\
		    data(model_data['developer'],1.135))[0]

    cosmo = mod.cosmo

    Da = np.zeros(len(z))
    H = np.zeros(len(z))
    for i in range(len(z)):
        Da[i] = cosmo.angular_diameter_distance(z[i])
    Da_z1 =  (cosmo.angular_diameter_distance(1.135))
    Da_z5 =  (cosmo.angular_diameter_distance(5.74))

    fig, ax = plt.subplots(figsize = (25,14))

    plt.plot(z,Da,'k--',linewidth=.9)

    # !!!
    plt.plot([0.38],[1512/(1+0.38)],'D',color='purple',label=r'$\rm BOSS/eBOSS$',alpha=0.4)
    plt.errorbar([0.38],[1512/(1+0.38)],yerr=22,color='purple',capsize=5,alpha=0.4,)
    plt.plot([0.51],[1975/(1+0.51)],'D',color='purple',alpha=0.4)
    plt.errorbar([0.51],[1975/(1+0.51)],yerr=27,color='purple',capsize=5,alpha=0.4,)
    plt.plot([0.61],[2307/(1+0.61)],'D',color='purple',alpha=0.4)
    plt.errorbar([0.61],[2307/(1+0.61)],yerr=33,color='purple',capsize=5,alpha=0.4,)    
    

    plt.plot([0.65],[cosmo.angular_diameter_distance(0.65)],'D',color='r',label=r'$\rm DESI$',alpha=0.4)
    plt.errorbar([0.65],[cosmo.angular_diameter_distance(0.65)],yerr=0.82*cosmo.angular_diameter_distance(0.65)/100,color='r',capsize=5,alpha=0.4)

    desi = [0.75,0.85,0.95,1.15,1.25,1.35,1.45,1.65,1.75,1.85,1.96,2.12,2.28,2.43,2.59,2.75,2.91,3.07,3.23,3.39,3.55]
    err_desi = [0.69,0.69,0.73,0.89,0.94,0.96,1.5,1.59,1.90,2.88,4.64,4.71,2.69,1.95,2.18,2.46,2.86,3.40,4.21,5.29,7.10,10.46,15.91]
    for i in range(len(desi)):
        plt.plot([desi[i]],[cosmo.angular_diameter_distance(desi[i])],'D',color='r',alpha=0.4)
        plt.errorbar([desi[i]],[cosmo.angular_diameter_distance(desi[i])],yerr=err_desi[i]*cosmo.angular_diameter_distance(desi[i])/100,color='r',capsize=5,alpha=0.4)

    plt.plot([0.55],[cosmo.angular_diameter_distance(0.55)],'o',color='green',label=r'$\rm SPHEREx\,H\alpha$',alpha=0.4)
    plt.errorbar([0.55],[cosmo.angular_diameter_distance(0.55)],yerr=4.7*cosmo.angular_diameter_distance(0.55)/100,color='green',capsize=5,alpha=0.4,)

    patches = [mpatches.FancyBboxPatch(
        [1.07,1650],1.2-1.07, 1830-1650,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02))]
    
    collection = PatchCollection(patches, alpha=0.9, edgecolor=kitsune,facecolor='none')
    ax.add_collection(collection)

    plt.plot([1.135],[Da_z1],'o',color=kitsune,label=r'$\rm ULTRASAT$')
    plt.errorbar(1.135,Da_z1,yerr=4.1*Da_z1/100,color=kitsune,capsize=5,)

    plt.plot([1.90],[cosmo.angular_diameter_distance(1.90)],'o',color='green',alpha=0.4)
    plt.errorbar([1.90],[cosmo.angular_diameter_distance(1.90)],yerr=3.1*cosmo.angular_diameter_distance(1.90)/100,color='green',capsize=5,alpha=0.4,)

    # plt.plot([2.73],[cosmo.angular_diameter_distance(2.73)],'o',color='c',label=r'$\rm IMS3\, CO$')
    # plt.errorbar([2.73],[cosmo.angular_diameter_distance(2.73)],yerr=1.1*cosmo.angular_diameter_distance(2.73)/100,color='c',capsize=5,)


    plt.plot([2.84],[cosmo.angular_diameter_distance(2.84)],'o',color='b',label=r'$\rm COMAP2$',alpha=0.4)
    plt.errorbar([2.84],[cosmo.angular_diameter_distance(2.84)],yerr=5.2*cosmo.angular_diameter_distance(2.84)/100,color='b',capsize=5,alpha=0.4,)

    plt.plot([3.20],[cosmo.angular_diameter_distance(3.20)],'o',color='green',alpha=0.4)
    plt.errorbar([3.20],[cosmo.angular_diameter_distance(3.20)],yerr=4*cosmo.angular_diameter_distance(3.20)/100,color='green',capsize=5,alpha=0.4,)

    # plt.plot([4.01],[cosmo.angular_diameter_distance(4.01)],'o',color='c')
    # plt.errorbar([4.01],[cosmo.angular_diameter_distance(4.01)],yerr=1.1*cosmo.angular_diameter_distance(4.01)/100,color='c',capsize=5,)

    plt.plot([4.52],[cosmo.angular_diameter_distance(4.52)],'o',color='green',alpha=0.4)
    plt.errorbar([4.52],[cosmo.angular_diameter_distance(4.52)],yerr=7*cosmo.angular_diameter_distance(4.52)/100,color='green',capsize=5,alpha=0.4,)

    #plt.plot([5.30],[cosmo.angular_diameter_distance(5.30)],'o',color='c')
    #plt.errorbar([5.30],[cosmo.angular_diameter_distance(5.30)],yerr=1.5*cosmo.angular_diameter_distance(4.01)/100,color='c',capsize=5,)

    patches = [mpatches.FancyBboxPatch(
        [5.65,1150],5.8-5.65, 1320-1150,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02))]
    
    collection = PatchCollection(patches, alpha=0.5, edgecolor='k',facecolor='none')
    ax.add_collection(collection)

    plt.plot([5.7],[cosmo.angular_diameter_distance(5.7)],'o',color='grey',label=r'$\rm SPHEREx\, + ULTRASAT$')
    plt.errorbar(5.7,cosmo.angular_diameter_distance(5.7),yerr=2.8*Da_z5/100,color='grey',capsize=5,linestyle='--')

    plt.plot([5.74],[Da_z5],'o',color='k',label=r'$\rm SPHEREx\, Ly\alpha$')
    plt.errorbar(5.74,Da_z5,yerr=5.8*Da_z5/100,color='k',capsize=5,)

    plt.legend(fontsize=40,loc=1,markerscale=2)
    plt.xlabel(r'$z$',fontsize=47)
    plt.ylabel(r'$D_A(z)\,{\rm [Mpc]}$',fontsize=47)
    plt.xticks(fontsize=45)
    plt.yticks(fontsize=45)
    plt.ylim(1050,1890)
    #plt.ylim(0,2500)
    plt.xlim(0.3,6)
    #plt.xlim(0.,9)
    plt.tight_layout()
    plt.savefig(save_fig_dir + 'DA_all.pdf')
    plt.show()

    return 

def prior_on_spherex():

    astro = ['z0_Lyalpha','f0_Lyalpha','t_Lyalpha']

    fisher_us = [[107908034.31026058,-1579893712.0981967,99170374.71607204,-42.0550626566643,-33672307.140302755,-2.850759309522287e+16,26628.715457146416,170495899.45658666,327036276.5957075],
[-1579893712.0981967,23140686519.10109,-1452693624.95489,582.0309567247424,493000077.0916768,4.204429005104437e+17,-391708.5258562882,-2496249489.3734527,-4795707424.797494],
[99170374.71607207,-1452693624.9548898,91197467.9261914,-36.00215476901637,-30945754.5648389,-2.643950341789131e+16,24616.541321071563,156690282.2605156,301146516.9510186],
[-42.0550626566643,582.0309567247423,-36.00215476901637,0.2124470259550733,12.95309062837356,8895691854.577436,-0.00854598982008228,-69.83909592479358,-122.01130391895383],
[-33672307.140302755,493000077.09167683,-30945754.564838894,12.953090628373559,10507320.331474854,8895691854191547.0,-8309.394658108082,-53202618.65960785,-102050477.44797397],
[-2.8507593095222856e+16,4.2044290051044365e+17,-2.6439503417891308e+16,8895691854.577436,8895691854191545.0,1.3105571238723815e+28,-12113002924889.908,-4.504243439572226e+16,-8.910753197486542e+16],
[26628.715457146413,-391708.5258562882,24616.541321071563,-0.008545989820082282,-8309.394658108082,-12113002924889.912,7.892304010563655,42073.70114779437,82296.69459526648],
[170495899.45658666,-2496249489.373453,156690282.26051566,-69.8390959247936,-53202618.65960785,-4.5042434395722264e+16,42073.70114779437,269385480.5509899,516721131.91791344],
[327036276.5957075,-4795707424.797495,301146516.95101863,-122.01130391895383,-102050477.44797395,-8.910753197486544e+16,82296.69459526648,516721131.9179133,997349826.8810705]]
    
    fisher_spherex= [[37101.42200735512,-15370552.272697376,1921946.0892772682,-95.87903074312617,-3172734.2456136155,-9800781.19500189,3297.804672466154,1544590.6843691724,3036537.3455087463],
[-15370552.272697376,6371860945.105638,-797006943.2173152,27730.081230509044,1314418247.4761074,4078727310.87035,-1364692.7776885305,-639917919.9420037,-1260023839.049314],
[1921946.0892772682,-797006943.2173152,99708589.10611138,-2686.004846949351,-164356094.3904687,-511201497.8679478,170542.5652966904,80016987.64793757,157686222.1328999],
[-95.8790307431262,27730.081230509044,-2686.0048469493504,804.2468031068785,8002.986975357364,3083.53421497587,-11.371702697855742,-3912.4450396647644,-5114.76899367081],
[-3172734.245613616,1314418247.4761076,-164356094.3904687,8002.986975357362,271316942.14551014,838120852.6677353,-282011.61749735166,-132085947.40978815,-259670660.36670962],
[-9800781.195001893,4078727310.8703513,-511201497.8679479,3083.53421497587,838120852.6677351,4831257378.327812,-798621.4213982024,-408100783.0354253,-812067846.6243827],
[3297.804672466154,-1364692.7776885307,170542.56529669036,-11.37170269785574,-282011.6174973516,-798621.4213982022,302.91391593723273,137286.30005696716,269065.47471036],
[1544590.684369173,-639917919.9420037,80016987.64793757,-3912.4450396647644,-132085947.40978816,-408100783.0354253,137286.3000569672,64303829.62711403,126424610.00065374],
[3036537.3455087463,-1260023839.049314,157686222.13289988,-5114.76899367081,-259670660.36670962,-812067846.6243826,269065.47471036005,126424610.00065374,249562786.40800065]
]

    new_fisher = deepcopy(fisher_spherex)
    for i in range(len(astro)):
        for j in range(len(astro)):
            if i == j:
                new_fisher[i][j] = fisher_spherex[i][j] + fisher_us[i][j]

    inv = np.linalg.inv(new_fisher)

    print('\nError on DA/rs = ' + str(np.sqrt(inv[-1][-1])*8.22*100))    

    return 