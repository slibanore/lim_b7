# SL: last update 01/23/2023

from LIM_b7 import *
from LIM_b7.fiducial_pars import astrocosmo_dict

save_fig_dir = './results/nu_mass/' 


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
                        2.5*np.trapz(analitic_der*legendre(2)(model.mui_grid),model.mu,axis=0)\
                        if get_multipole == 2 else -1 

            ders.append(der_marg_rsd.value)

        if derive_pars[i] == 'Tbs8':
 
            analitic_der = (2.*Pkclus/(Tbs8+Tfs8*model.mui_grid**2))

            der_marg_rsd = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                        2.5*np.trapz(analitic_der*legendre(2)(model.mui_grid),model.mu,axis=0)\
                        if get_multipole == 2 else -1 
 
            ders.append(der_marg_rsd.value)
 
        if derive_pars[i] == 'Pshot':
 
            ders.append(np.ones(model.k.shape))

 
        if derive_pars[i] == 'fNL':
 
            if model.camb_pars.num_nu_massive != 0:
                Tk = model.transfer_cb(model.k.value)
            else:
                Tk = model.transfer_m(model.k.value)
            k2Tk = model.k**2*Tk

            fNL = 0.
            factor = 3.*model.camb_pars.omegam*(100.*model.hubble*(u.km/u.s/u.Mpc))**2./   \
                        (cu.c.to(u.km/u.s)**2.*k2Tk)
                        
            analitic_der = (2.*Pkclus*(model.bavg-1.)*factor/(model.bavg+(model.bavg-1)*fNL*factor+
                    fs8(model,z)/s8(model,z)))

            der_marg_rsd = 0.5*np.trapz(analitic_der ,model.mu,axis=0) if get_multipole == 0 else\
                        2.5*np.trapz(analitic_der*legendre(2)(model.mui_grid),model.mu,axis=0)\
                        if get_multipole == 2 else -1 
 
            ders.append(der_marg_rsd.value)

 
        if derive_pars[i] == 'sNL':
 
            sNL = model.sigma_NL

            analitic_der =  (-2.*Pkclus*sNL*(model.ki_grid*model.mui_grid)**2./\
                (1.+0.5*(sNL*model.ki_grid*model.mui_grid)**2.))
                
            der_marg_rsd = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                        2.5*np.trapz(analitic_der*legendre(2)(model.mui_grid),model.mu,axis=0)\
                        if get_multipole == 2 else -1 
 
            ders.append(der_marg_rsd.value) 


        if derive_pars[i] == 'alpha_par':

            Pkclus_up = Pk_AP(model,z, 1.+delta, 1.)
            Pkclus_low = Pk_AP(model,z, 1.-delta, 1.)

            analitic_der = (Pkclus_up - Pkclus_low).value / (2*delta)
            derivative = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                    2.5*np.trapz(analitic_der*legendre(2\
                    )(model.mui_grid),model.mu,axis=0)\
                    if get_multipole == 2 else -1  

            ders.append(derivative) 

        if derive_pars[i] == 'alpha_perp':

            Pkclus_up = Pk_AP(model,z, 1., 1.+delta)
            Pkclus_low = Pk_AP(model,z, 1., 1.-delta)

            analitic_der = (Pkclus_up - Pkclus_low).value / (2*delta)
            derivative = 0.5*np.trapz(analitic_der,model.mu,axis=0) if get_multipole == 0 else\
                    2.5*np.trapz(analitic_der*legendre(2\
                    )(model.mui_grid),model.mu,axis=0)\
                if get_multipole == 2 else -1 

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

    if derive_par == 'n_B': 	
        use_n_B_low = n_B_fid*(1. - delta)
        use_n_B_up = n_B_fid*(1. + delta)
    else:
        use_n_B_low = n_B_fid 
        use_n_B_up = n_B_fid
    
    if derive_par == 'sigma_B_0': 	
        use_sigma_B_0_low = sigma_B_0_fid*(1. - delta)
        use_sigma_B_0_up = sigma_B_0_fid*(1. + delta)
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

    if derive_par == 'fNL':
        use_fNL_low = fNL_fid*(1. - delta)
        use_fNL_up = fNL_fid*(1. + delta)
    else:
        use_fNL_low = fNL_fid
        use_fNL_up = fNL_fid


    low_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_low, \
            beta=use_astro_beta_low, dMF = astro_dMF_fid, \
            sig_SFR=use_astro_sig_SFR_low,SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        low_hmf = dict(A_tinker = use_Atink_low, a_tinker = use_atink_low, b_tinker = use_btink_low, c_tinker = use_ctink_low)

    elif hmf_fid == 'NG_Riotto':
        low_hmf = dict(A_tinker = use_Atink_low, a_tinker = use_atink_low, b_tinker = use_btink_low, c_tinker = use_ctink_low,
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
        n_B = use_n_B_low,
        sigma_B_0 = use_sigma_B_0_low)}

    up_model_par = {**astrocosmo_dict(developer,z), **dict(alpha=use_astro_alpha_up, \
            beta=use_astro_beta_up, dMF = astro_dMF_fid,\
            sig_SFR=use_astro_sig_SFR_up,SFR_file=SFR_file_fid)}

    if hmf_fid == 'Tinker':
        up_hmf = dict(A_tinker = use_Atink_up, a_tinker = use_atink_up, b_tinker = use_btink_up, c_tinker = use_ctink_up)

    elif hmf_fid == 'NG_Riotto':
        up_hmf = dict(A_tinker = use_Atink_up, a_tinker = use_atink_up, b_tinker = use_btink_up, c_tinker = use_ctink_up,
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
        n_B = use_n_B_up,
        sigma_B_0 = use_sigma_B_0_up)}

    low = update_Pkline(obs_pars(z),low_data)[0]
    up = update_Pkline(obs_pars(z),up_data)[0]

    return low, up


def fisher_Pkline(model_id, developer, \
        derive_model = ['a1','a2','a3','a4','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','a_tinker'],\
        derive_Pkline = ['Tfs8','Tbs8','Pshot','sNL','alpha_par','alpha_perp'],	multipole = [0,2],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True, \
	    a1_prior = 0.03, \
        a_t_prior = 0.2):

    create_dir(save_fig_dir + 'Pkline/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_Pkline/' + developer + '/')

    create_dir(save_fig_dir + 'Pkline/' + developer + '/' + model_id + '/')
    create_dir(save_fig_dir + 'fisher_Pkline/' + developer + '/' + model_id + '/')

    survey_name, model_name, obs_pars, model_z, z_vals  = set_survey(model_id,developer,estimator='Pkline')

    derive_pars = derive_model + derive_Pkline
    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    for z in z_vals:

        direct = save_fig_dir + 'Pkline/' + developer + '/' + model_id + '/z' + str(z)
        if developer == 'axions':
            direct += '_' + str(f_axion_fid) 
        if developer == 'PMF':
            direct += '_' + str(n_B_fid) + ',' + str(sigma_B_0_fid)
        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_Pkline/' + developer  + '/' + model_id + '/z' + str(z)
        if developer == 'axions':
            direct_fisher += '_' + str(f_axion_fid) 
        if developer == 'PMF':
            direct_fisher += '_' + str(n_B_fid) + ',' + str(sigma_B_0_fid)

        create_dir(direct)
        create_dir(direct_fisher)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')

        model = model_z(z)

        fid_par, name_fid = get_fiducials(derive_model, developer, z)

        vecpar, name_vec = get_fiducials_Pkline(derive_Pkline, model, z)

        fiducials = fid_par + vecpar
        names = name_fid + name_vec

        fiducial_k = model.k.value 

        fiducial_k = np.asarray(fiducial_k) 
        Pkclus_fid = Pk_AP(model,z,1.,1.)

        integrand_cov = (model.Pk+model.Pnoise)/model.Nmodes**0.5 # Pnoise includes Vvox --> N_z_shells
        cov_02 = (5./2.*np.trapz(integrand_cov**2*legendre(2)(model.mui_grid)**2,model.mu,axis=0)).value

        cov = [[[],cov_02],[cov_02,[]]]

        der_Pk_multi = [[],[]]

        for l in range(len(multipole)):

            print('--------- Doing L = ' + str(multipole[l]) + ' ---------')
            print('Running Pk fiducial')
            fiducial_Pk = 0.5*np.trapz(Pkclus_fid, model.mu, axis=0) if multipole[l] == 0 else 2.5*np.trapz(Pkclus_fid*legendre(2)(model.mui_grid),model.mu,axis=0) if multipole[l] == 2 else -1 

            cov[l][l] = np.zeros(len(fiducial_Pk))
                
            sigma2 = 0.5*np.trapz(integrand_cov**2,model.mu,axis=0) if multipole[l] == 0 else\
                    25/2.*np.trapz(integrand_cov**2*legendre(2)(model.mui_grid)*\
                    legendre(2)(model.mui_grid),model.mu,axis=0) if multipole[l] == 2 else -1  

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
                    2.5*np.trapz(Pkclus_low*legendre(2)(model_low.mui_grid),model_low.mu,axis=0)\
                    if multipole[l] == 2 else -1 

                    k_low = model_low.k

                    if save_Pkline_flag:
                        np.savetxt(filename_low,(k_low,Pk_low))

                all_k_low.append(np.asarray(k_low)) 
                all_Pk_low.append(np.asarray(Pk_low))

                if os.path.exists(filename_up) and import_allowed:
                    k_up, Pk_up = import_Pkline(filename_up, 'Import Pk up')
                else:
                    Pk_up = 0.5*np.trapz(Pkclus_up, model_up.mu, axis=0) if multipole[l] == 0 else\
                            2.5*np.trapz(Pkclus_up*legendre(2)(model_up.mui_grid),model_up.mu,axis=0)\
                            if multipole[l] == 2 else -1 
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
        for kval in range(len(cov[0][0])):
            use_cov = np.asarray([[cov[0][0][kval],cov[0][1][kval]],\
                        [cov[1][0][kval],cov[1][1][kval]]])
            inv_cov.append(np.linalg.inv(use_cov))
            
        for a in range(len(derive_pars)):
            multi_der_Pk_ordered.append([])	
            for kval in range(len(der_Pk_multi[0][a])):
                multi_der_Pk_ordered[a].append([der_Pk_multi[0][a][kval],der_Pk_multi[1][a][kval]])

        inv_cov = np.asarray(inv_cov)
        multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)

        #return inv_cov, multi_der_Pk_ordered
        Fisher_k = []
        Fisher_z = np.zeros((len(der_Pk),len(der_Pk)))
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
                    Fisher_z[i,j] +=  Fisher_k[t][i,j]

        Fisher += Fisher_z 
         
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkline/' + model_id + 'fisher_tot' 
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

    inv_Fisher = np.linalg.inv(Fisher)
    sigma_tot = np.zeros(len(derive_pars))
    for i in range(len(derive_pars)):
        try:
            sigma_tot[i] = np.sqrt(inv_Fisher[i,i]).value
        except:
            sigma_tot[i] = np.sqrt(inv_Fisher[i,i])

    for i in range(len(fiducials)):
        if type(fiducials[i]) is not np.float64 and type(fiducials[i]) is not float:
            fiducials[i] = fiducials[i].value

    if plot_sigma_flag:

        print('Doing contour plots')

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
            legend_labels=[r'$PS\, wide$'],legend_loc='upper right',\
            contour_colors=['green'], line_args=[{'ls':'--', 'color':'green'}],
            #title_limit=1,
            markers={'x2':0})

        plt.savefig(save_fig_dir + 'ellipse_' + model_id + '_pk.pdf')

        plt.suptitle(survey_name,fontsize=25)
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot
