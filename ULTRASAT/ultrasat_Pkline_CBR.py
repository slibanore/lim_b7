# SL: last update 01/23/2023

from LIM_b7 import *
from LIM_b7.ultrasat_fiducial_pars import astrocosmo_dict

save_fig_dir = './results_ultrasat/new/CBR/' 

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



# compute power spectrum with redshift space distortions
def Pk_cross(model, z, bLy, bg, alpha_par = 1., alpha_perp = 1.):

    sNL = model.sigma_NL

    F = alpha_par/alpha_perp
    prefac = 1./alpha_perp**2./alpha_par # volume modification

    kprime = np.zeros((len(model.mu),len(model.k)))
    mu_prime = model.mui_grid/F/np.sqrt(1.+model.mui_grid**2.*(F**-2-1.))
    for imu in range(model.nmu):
        kprime[imu,:] = (model.k /alpha_perp)*np.sqrt(1.+model.mu[imu]**2*(F**-2-1.)) 

    #Obtain the corresponding P_m for the "real" k
    Pmprime = log_interp1d(model.k,model.Pm[0,:])(kprime)*model.Pm.unit

    sigma = model.cosmo.get_sigma8()[::-1]
    sigma_interp = interp1d(model.zcosmo,sigma)
    s8 = sigma_interp(z)	

    fsigma = model.cosmo.get_fsigma8()[::-1]
    fsigma_interp = interp1d(model.zcosmo,fsigma)

    #Apply RSD for the "real" mu
    kaiser = (model.Tmean*bLy*s8 + model.Tmean*fsigma_interp(z)*mu_prime**2)*(bg*s8 + fsigma_interp(z)*mu_prime**2)
    loren = 1./(1.+0.5*(kprime*mu_prime*sNL.value)**2)

    #Units change depend on the survey. Autoconsistent definition in terms of Mpc/h
    Pkclus = (kaiser*Pmprime*loren**2/s8**2) 

    Pkclus = prefac*(Pkclus*(model.Wkmin*model.Wkmax)) 

    return Pkclus

# compute power spectrum with redshift space distortions
def Pk_gal(model, z,  bg, alpha_par = 1., alpha_perp = 1.):

    sNL = model.sigma_NL

    F = alpha_par/alpha_perp
    prefac = 1./alpha_perp**2./alpha_par # volume modification

    kprime = np.zeros((len(model.mu),len(model.k)))
    mu_prime = model.mui_grid/F/np.sqrt(1.+model.mui_grid**2.*(F**-2-1.))
    for imu in range(model.nmu):
        kprime[imu,:] = (model.k /alpha_perp)*np.sqrt(1.+model.mu[imu]**2*(F**-2-1.)) 

    #Obtain the corresponding P_m for the "real" k
    Pmprime = log_interp1d(model.k,model.Pm[0,:])(kprime)*model.Pm.unit

    sigma = model.cosmo.get_sigma8()[::-1]
    sigma_interp = interp1d(model.zcosmo,sigma)
    s8 = sigma_interp(z)	

    fsigma = model.cosmo.get_fsigma8()[::-1]
    fsigma_interp = interp1d(model.zcosmo,fsigma)

    #Apply RSD for the "real" mu
    kaiser = (bg*s8 + fsigma_interp(z)*mu_prime**2)**2
    loren = 1./(1.+0.5*(kprime*mu_prime*sNL.value)**2)

    #Units change depend on the survey. Autoconsistent definition in terms of Mpc/h
    Pkclus = (kaiser*Pmprime*loren**2/s8**2) 

    Pkclus = prefac*(Pkclus*(model.Wkmin*model.Wkmax)) 

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

    if derive_par == 'w0':
        w0_low = -(1-delta)
        w0_up = -(1+delta)
    else:
        w0_low = -1
        w0_up = -1

    if derive_par == 'wa':
        wa_low = 0.
        wa_up = delta
    else:
        wa_low = 0.
        wa_up = 0.

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


    low_cosmo = dict(
	f_NL=0, H0=67.67, cosmomc_theta=None,
	ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
	num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
	YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
	TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
	theta_H0_range=[10, 100], w=w0_low, wa=wa_low, cs2=1.0, 
	dark_energy_model='ppf',As=2.105e-09, 
    ns=0.96,r=0.0, nt=None, ntrun=0.0, 
	pivot_scalar=0.05, pivot_tensor=0.05,
	parameterization=2,halofit_version='mead')

#    low_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_low, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_low, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_low,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)}
    low_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_low,B_lya=use_B_Lyalpha_low, D_lya=use_D_Lyalpha_low, dndL_Lcut=use_astro_Lcut_low, SFR_file=SFR_file_fid)}
                     
    low_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z))

    low_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, 
        cosmo_input_camb = low_cosmo,
        hmf_pars = low_hmf,
        model_par = low_model_par,
        sigma_scatter = use_astro_sig_scatter_low,
        dndL_Lcut=use_astro_Lcut_low)}

    up_cosmo = dict(
	f_NL=0, H0=67.67, cosmomc_theta=None,
	ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
	num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
	YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
	TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
	theta_H0_range=[10, 100], w=w0_up, wa=wa_up, cs2=1.0, 
	dark_energy_model='ppf',As=2.105e-09, 
    ns=0.96,r=0.0, nt=None, ntrun=0.0, 
	pivot_scalar=0.05, pivot_tensor=0.05,
	parameterization=2,halofit_version='mead')

#    up_model_par = {**astrocosmo_dict(developer,z), **dict(csi=csi_Lyalpha_fid,z0=use_z0_Lyalpha_up, zeta=zeta_Lyalpha_fid, f0 = use_f0_Lyalpha_up, SFR0 = SFR0_Lyalpha_fid, t = use_t_Lyalpha_up,  dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)}
    up_model_par = {**astrocosmo_dict(developer,z), **dict(A_lya=use_A_Lyalpha_up,B_lya=use_B_Lyalpha_up, D_lya=use_D_Lyalpha_up, dndL_Lcut=use_astro_Lcut_up, SFR_file=SFR_file_fid)}
                     
    up_hmf = dict(A_tinker = A_tinker_fid(z), a_tinker = a_tinker_fid(z), b_tinker = b_tinker_fid(z), c_tinker = c_tinker_fid(z))

    up_data = {**astrocosmo_dict(developer,z), **dict(developer = developer, 
        cosmo_input_camb = up_cosmo,
        hmf_pars = up_hmf,
        model_par = up_model_par,
        sigma_scatter = use_astro_sig_scatter_up,
        dndL_Lcut=use_astro_Lcut_up)}

    low = update_Pkline(obs_pars(z),low_data)[0]
    up = update_Pkline(obs_pars(z),up_data)[0]

    return low, up


def fisher_Pkline(z, developer, \
        derive_model = ['w0','wa'],\
        derive_Pkline = ['Tfs8','Tbs8'],multipole = [0,2,4],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):

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
    if 2 in multipole:
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
    else:
        cov = [[[]]]
        der_Pk_multi = [[]]

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
                    denom = (2*delta*fid_par[id_der]) if fid_par[id_der] != 0. else delta
                    der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)

            elif derive_pars[i] in derive_Pkline:
                id_der = derive_Pkline.index(derive_pars[i])
                all_k_up.append(model.k)
                der_Pk.append(pkprop_ders[id_der])

            der_Pk_multi[l] = der_Pk

    inv_cov = []
    multi_der_Pk_ordered = []
    if 2 in multipole:
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
    else:
        for kval in range(len(cov[0][0])):
            use_cov = np.asarray(cov[0][0][kval])
            inv_cov.append(use_cov**-1)
        for a in range(len(derive_pars)):
            multi_der_Pk_ordered.append([])
            for kval in range(len(der_Pk_multi[0][a])):
                multi_der_Pk_ordered[a].append(der_Pk_multi[0][a][kval])        

    inv_cov = np.asarray(inv_cov)
    multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)
    
    Fisher_k = []
    Fisher = np.zeros((len(der_Pk),len(der_Pk)))
    for t in range(len(fiducial_k)):
        Fisher_k.append(np.zeros((len(der_Pk),len(der_Pk))))		
        for i in range(len(der_Pk)):
            for j in range(len(der_Pk)):
                if len(multipole) > 1:
                    Fisher_k[t][i,j] = np.linalg.multi_dot(\
                        [multi_der_Pk_ordered[i][t],inv_cov[t],\
                        multi_der_Pk_ordered[j][t]])
                else:
                    Fisher_k[t][i,j] = multi_der_Pk_ordered[i][t]*\
                        inv_cov[t]*multi_der_Pk_ordered[j][t]

    Fisher_k = np.asarray(Fisher_k)
    for i in range(len(der_Pk)):
        for j in range(len(der_Pk)):
            for t in range(len(fiducial_k)):
                Fisher[i,j] +=  Fisher_k[t][i,j]
        
    
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkline/fisher_Lya' 
        save_fisher(filename, derive_pars,Fisher)

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

        print('\nError on w0 = ' + str(np.sqrt(inv_Fisher[0][0])))
        print('Error on wa = ' + str(np.sqrt(inv_Fisher[1][1]) ) + '\n')

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

        g.plot_2d([plot_dist], names, filled=True,\
            #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers=[{'x2':0}])

        plt.tight_layout()
        plt.savefig(save_fig_dir + 'ellipse_DE.pdf')
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot


def fisher_Pkcross(z_arr, Nbin, developer, bg =2.,\
        derive_model = ['w0','wa'],\
        derive_Pkline = ['Tfs8','Tbs8'],multipole = [0,2,4],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):

    fsky = 1/3.
    create_dir(save_fig_dir + 'Pkcross/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_Pkcross/' + developer + '/')

    create_dir(save_fig_dir + 'Pkcross/' + developer + '/' )
    create_dir(save_fig_dir + 'fisher_Pkcross/' + developer + '/')

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT_CBR(z_val,Nbin))
    obs_pars_line = lambda z_val: deepcopy(obs_params_ULTRASAT((z_arr[0]+z_arr[-1])/2.))
    #obs_pars = lambda z_val: deepcopy(obs_params_SPHEREx(z_val))

    derive_pars = derive_model + derive_Pkline

    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    model_line = update_Pkline(obs_pars_line((z_arr[0]+z_arr[-1])/2.), astrocosmo_dict(developer,(z_arr[0]+z_arr[-1])/2.))[0]

    # tilde_Pkline_fid = Pk_AP(model_line,(z_arr[0]+z_arr[-1])/2., 1.,1.) + model_line.Pnoise

    for z in z_arr:

        model = update_Pkline(obs_pars(z), astrocosmo_dict(developer,z))[0]

        Fisher_z = np.zeros((len(derive_pars),len(derive_pars)))


        direct = save_fig_dir + 'Pkcross/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct += '_' + str(f_axion_fid) 

        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_Pkcross/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct_fisher += '_' + str(f_axion_fid) 

        create_dir(direct)
        create_dir(direct_fisher)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')


        fid_par, name_fid = get_fiducials(derive_model, developer, z)

        vecpar, name_vec = get_fiducials_Pkline(derive_Pkline, model, z)

        if 'bLy' and 'bg' in derive_model:
            fiducials = [model.bavg[0],bg] + fid_par + vecpar
            names = [r'$b_{\rm Ly\alpha}$',r'$b_g$'] + name_fid + name_vec
        else:    
            fiducials = fid_par + vecpar
            names = name_fid + name_vec

        fiducial_k = model.k.value 

        fiducial_k = np.asarray(fiducial_k) 
        Pkline_fid = Pk_AP(model_line,z, 1.,1.) + model_line.Pnoise    
        Pkcross_fid = Pk_cross(model,z,model.bavg[0],bg, 1.,1.)
        Pkgal_fid = Pk_gal(model,z,bg, 1.,1.)

        integrand_cov = (Pkcross_fid**2 + (Pkline_fid+model.Pnoise)*(Pkgal_fid+model.nbar**-1))/(model.Nmodes**0.5) 
        #integrand_cov = ((Pkcross_fid**2 + tilde_Pkline_fid*(Pkgal_fid+model.nbar**-1))/(2*model.Nmodes))**0.5 

        if 2 in multipole:
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
        else:
            cov = [[[]]]
            der_Pk_multi = [[]]

        for l in range(len(multipole)):

            print('--------- Doing L = ' + str(multipole[l]) + ' ---------')
            print('Running Pk fiducial')
            fiducial_Pk = 0.5*np.trapz(Pkcross_fid, model.mu, axis=0) if multipole[l] == 0 else (2*multipole[l]+1)/2.*np.trapz(Pkcross_fid*legendre(multipole[l])(model.mui_grid),model.mu,axis=0) 

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

                filename_low = direct + 'Pkcross_' + str(multipole[l])+ '_low_' + str(derive_model[i])			
                filename_up = direct + 'Pkcross_' + str(multipole[l]) + 'up_' + str(derive_model[i])

                model_low, model_up = set_der_par(derive_model[i],z,developer,obs_pars)

                if os.path.exists(filename_low) and import_allowed:
                    k_low, Pkcross_low = import_Pkline(filename_low, '\nDoing derivative wrt: ' + str(derive_model[i]) + '\nImport Pk low')
                else:
                    k_low = model_low.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_low = model.bavg[0]*(1-delta)
                        else:
                            bLy_low = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_low = bg*(1-delta)
                        else:
                            bg_low = bg
                        Pkclus_low = Pk_cross(model,z,bLy_low,bg_low, 1.,1.)

                    else:

                        Pkclus_low = Pk_cross(model_low,z,model.bavg[0],bg,1.,1.)

                    Pkcross_low = 0.5*np.trapz(Pkclus_low , model_low.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_low*legendre(multipole[l])(model_low.mui_grid),model_low.mu,axis=0)

                    if save_Pkline_flag:
                        np.savetxt(filename_low,(k_low,Pkcross_low))

                all_k_low.append(np.asarray(k_low)) 
                all_Pk_low.append(np.asarray(Pkcross_low))

                if os.path.exists(filename_up) and import_allowed:
                    k_up, Pkcross_up = import_Pkline(filename_up, 'Import Pk up')

                else:
                    k_up = model_up.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_up = model.bavg[0]*(1+delta)
                        else:
                            bLy_up = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_up = bg*(1+delta)
                        else:
                            bg_up = bg

                        Pkclus_up = Pk_cross(model,z,bLy_up,bg_up, 1.,1.)

                    else:

                        Pkclus_up = Pk_cross(model_up,z,model.bavg[0],bg,1.,1.)

                    Pkcross_up = 0.5*np.trapz(Pkclus_up, model_up.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_up*legendre(multipole[l])(model_up.mui_grid),model_up.mu,axis=0) 

                    if save_Pkline_flag:
                        np.savetxt(filename_up, (k_up,Pkcross_up))

                all_k_up.append(np.asarray(k_up)) 
                all_Pk_up.append(np.asarray(Pkcross_up))


                pkprop_ders = analytical_derivatives(derive_Pkline, model,z,Pkcross_fid,multipole[l])

            der_k = []
            der_Pk = []
            for i in range(len(derive_pars)):

                if derive_pars[i] in derive_model:
                    der_k.append([])
                    der_Pk.append([])
                    id_der = derive_model.index(derive_pars[i])

                    if derive_pars[i] == 'bLy':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*model.bavg[0])
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    elif derive_pars[i] == 'bg':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*bg)
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    else:      
                        if 'bLy' and 'bg' in derive_model:
                            denom = (2*delta*fid_par[id_der-2]) if fid_par[id_der-2] != 0. else delta
                        else:
                            denom = (2*delta*fid_par[id_der]) if fid_par[id_der] != 0. else delta
                        for t in range(len(fiducial_k)):
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)

                elif derive_pars[i] in derive_Pkline:
                    id_der = derive_Pkline.index(derive_pars[i])
                    all_k_up.append(model.k)
                    der_Pk.append(pkprop_ders[id_der])

                der_Pk_multi[l] = der_Pk

        inv_cov = []
        multi_der_Pk_ordered = []
        if 2 in multipole:
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
        else:
            for kval in range(len(cov[0][0])):
                use_cov = np.asarray(cov[0][0][kval])
                inv_cov.append(use_cov**-1)
            for a in range(len(derive_pars)):
                multi_der_Pk_ordered.append([])
                for kval in range(len(der_Pk_multi[0][a])):
                    multi_der_Pk_ordered[a].append(der_Pk_multi[0][a][kval])        

        inv_cov = np.asarray(inv_cov)
        multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)
        
        Fisher_k = []
        Fisher = np.zeros((len(der_Pk),len(der_Pk)))
        for t in range(len(fiducial_k)):
            Fisher_k.append(np.zeros((len(der_Pk),len(der_Pk))))		
            for i in range(len(der_Pk)):
                for j in range(len(der_Pk)):
                    if len(multipole) > 1:
                        Fisher_k[t][i,j] = np.linalg.multi_dot(\
                            [multi_der_Pk_ordered[i][t],inv_cov[t],\
                            multi_der_Pk_ordered[j][t]])
                    else:
                        Fisher_k[t][i,j] = multi_der_Pk_ordered[i][t]*\
                            inv_cov[t]*multi_der_Pk_ordered[j][t]

        Fisher_k = np.asarray(Fisher_k)
        for i in range(len(der_Pk)):
            for j in range(len(der_Pk)):
                for t in range(len(fiducial_k)):
                    Fisher_z[i,j] +=  Fisher_k[t][i,j]
        
        Fisher += Fisher_z * fsky
    
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkcross/fisher_Lya' 
        save_fisher(filename, derive_pars,Fisher)

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

        print('\nError on w0 = ' + str(np.sqrt(inv_Fisher[0][0])))
        print('Error on wa = ' + str(np.sqrt(inv_Fisher[1][1]) ) + '\n')

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

        g.plot_2d([plot_dist], names, filled=True,\
            #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers=[{'x2':0}])

        plt.tight_layout()
        plt.savefig(save_fig_dir + 'ellipse_DE.pdf')
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot


def fisher_Pkgal(z_arr, Nbin, developer, bg =2.,\
        derive_model = ['w0','wa'],\
        derive_Pkline = ['Tfs8','Tbs8'],multipole = [0,2,4],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):

    fsky = 1/3.
    create_dir(save_fig_dir + 'Pkgal/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_Pkgal/' + developer + '/')

    create_dir(save_fig_dir + 'Pkgal/' + developer + '/' )
    create_dir(save_fig_dir + 'fisher_Pkgal/' + developer + '/')

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT_CBR(z_val,Nbin))

    derive_pars = derive_model + derive_Pkline

    Fisher = np.zeros((len(derive_pars),len(derive_pars)))


    for z in z_arr:

        model = update_Pkline(obs_pars(z), astrocosmo_dict(developer,z))[0]

        Fisher_z = np.zeros((len(derive_pars),len(derive_pars)))

        direct = save_fig_dir + 'Pkgal/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct += '_' + str(f_axion_fid) 

        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_Pkgal/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct_fisher += '_' + str(f_axion_fid) 

        create_dir(direct)
        create_dir(direct_fisher)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')


        fid_par, name_fid = get_fiducials(derive_model, developer, z)

        vecpar, name_vec = get_fiducials_Pkline(derive_Pkline, model, z)

        if 'bLy' and 'bg' in derive_model:
            fiducials = [model.bavg[0],bg] + fid_par + vecpar
            names = [r'$b_{\rm Ly\alpha}$',r'$b_g$'] + name_fid + name_vec
        else:    
            fiducials = fid_par + vecpar
            names = name_fid + name_vec

        fiducial_k = model.k.value 

        fiducial_k = np.asarray(fiducial_k) 
        Pkgal_fid = Pk_gal(model,z,bg, 1.,1.)

        integrand_cov = (Pkgal_fid+model.nbar**-1)/(model.Nmodes**0.5) 

        if 2 in multipole:
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
        else:
            cov = [[[]]]
            der_Pk_multi = [[]]

        for l in range(len(multipole)):

            print('--------- Doing L = ' + str(multipole[l]) + ' ---------')
            print('Running Pk fiducial')
            fiducial_Pk = 0.5*np.trapz(Pkgal_fid, model.mu, axis=0) if multipole[l] == 0 else (2*multipole[l]+1)/2.*np.trapz(Pkgal_fid*legendre(multipole[l])(model.mui_grid),model.mu,axis=0) 

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

                filename_low = direct + 'Pkgal_' + str(multipole[l])+ '_low_' + str(derive_model[i])			
                filename_up = direct + 'Pkgal_' + str(multipole[l]) + 'up_' + str(derive_model[i])

                model_low, model_up = set_der_par(derive_model[i],z,developer,obs_pars)

                if os.path.exists(filename_low) and import_allowed:
                    k_low, Pkgal_low = import_Pkline(filename_low, '\nDoing derivative wrt: ' + str(derive_model[i]) + '\nImport Pk low')
                else:
                    k_low = model_low.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_low = model.bavg[0]*(1-delta)
                        else:
                            bLy_low = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_low = bg*(1-delta)
                        else:
                            bg_low = bg
                        Pkclus_low = Pk_gal(model,z,bg_low, 1.,1.)

                    else:

                        Pkclus_low = Pk_gal(model_low,z,bg,1.,1.)

                    Pkgal_low = 0.5*np.trapz(Pkclus_low , model_low.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_low*legendre(multipole[l])(model_low.mui_grid),model_low.mu,axis=0)

                    if save_Pkline_flag:
                        np.savetxt(filename_low,(k_low,Pkgal_low))

                all_k_low.append(np.asarray(k_low)) 
                all_Pk_low.append(np.asarray(Pkgal_low))

                if os.path.exists(filename_up) and import_allowed:
                    k_up, Pkgal_up = import_Pkline(filename_up, 'Import Pk up')

                else:
                    k_up = model_up.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_up = model.bavg[0]*(1+delta)
                        else:
                            bLy_up = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_up = bg*(1+delta)
                        else:
                            bg_up = bg

                        Pkclus_up = Pk_gal(model,z,bg_up, 1.,1.)

                    else:

                        Pkclus_up = Pk_gal(model_up,z,bg,1.,1.)

                    Pkgal_up = 0.5*np.trapz(Pkclus_up, model_up.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_up*legendre(multipole[l])(model_up.mui_grid),model_up.mu,axis=0) 

                    if save_Pkline_flag:
                        np.savetxt(filename_up, (k_up,Pkgal_up))

                all_k_up.append(np.asarray(k_up)) 
                all_Pk_up.append(np.asarray(Pkgal_up))


                pkprop_ders = analytical_derivatives(derive_Pkline, model,z,Pkgal_fid,multipole[l])

            der_k = []
            der_Pk = []
            for i in range(len(derive_pars)):

                if derive_pars[i] in derive_model:
                    der_k.append([])
                    der_Pk.append([])
                    id_der = derive_model.index(derive_pars[i])

                    if derive_pars[i] == 'bLy':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*model.bavg[0])
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    elif derive_pars[i] == 'bg':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*bg)
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    else:      
                        if 'bLy' and 'bg' in derive_model:
                            denom = (2*delta*fid_par[id_der-2]) if fid_par[id_der-2] != 0. else delta
                        else:
                            denom = (2*delta*fid_par[id_der]) if fid_par[id_der] != 0. else delta
                        for t in range(len(fiducial_k)):
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)

                elif derive_pars[i] in derive_Pkline:
                    id_der = derive_Pkline.index(derive_pars[i])
                    all_k_up.append(model.k)
                    der_Pk.append(pkprop_ders[id_der])

                der_Pk_multi[l] = der_Pk

        inv_cov = []
        multi_der_Pk_ordered = []
        if 2 in multipole:
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
        else:
            for kval in range(len(cov[0][0])):
                use_cov = np.asarray(cov[0][0][kval])
                inv_cov.append(use_cov**-1)
            for a in range(len(derive_pars)):
                multi_der_Pk_ordered.append([])
                for kval in range(len(der_Pk_multi[0][a])):
                    multi_der_Pk_ordered[a].append(der_Pk_multi[0][a][kval])        

        inv_cov = np.asarray(inv_cov)
        multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)
        
        Fisher_k = []
        Fisher = np.zeros((len(der_Pk),len(der_Pk)))
        for t in range(len(fiducial_k)):
            Fisher_k.append(np.zeros((len(der_Pk),len(der_Pk))))		
            for i in range(len(der_Pk)):
                for j in range(len(der_Pk)):
                    if len(multipole) > 1:
                        Fisher_k[t][i,j] = np.linalg.multi_dot(\
                            [multi_der_Pk_ordered[i][t],inv_cov[t],\
                            multi_der_Pk_ordered[j][t]])
                    else:
                        Fisher_k[t][i,j] = multi_der_Pk_ordered[i][t]*\
                            inv_cov[t]*multi_der_Pk_ordered[j][t]

        Fisher_k = np.asarray(Fisher_k)
        for i in range(len(der_Pk)):
            for j in range(len(der_Pk)):
                for t in range(len(fiducial_k)):
                    Fisher_z[i,j] +=  Fisher_k[t][i,j]
        
        Fisher += Fisher_z *fsky
    
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkgal/fisher_gal' 
        save_fisher(filename, derive_pars,Fisher)

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

        print('\nError on w0 = ' + str(np.sqrt(inv_Fisher[0][0])))
        print('Error on wa = ' + str(np.sqrt(inv_Fisher[1][1]) ) + '\n')

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

        g.plot_2d([plot_dist], names, filled=True,\
            #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers=[{'x2':0}])

        plt.tight_layout()
        plt.savefig(save_fig_dir + 'ellipse_DE.pdf')
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot


def fisher_Pkline_CBR(z_arr, Nbin, developer, bg =2.,\
        derive_model = ['w0','wa'],\
        derive_Pkline = ['Tfs8','Tbs8'],multipole = [0,2,4],\
    	save_Pkline_flag = False, \
        save_fisher_flag = False, \
	    plot_sigma_flag = False, \
        import_allowed = True):

    create_dir(save_fig_dir + 'Pkcross/' + developer + '/')
    create_dir(save_fig_dir + 'fisher_Pkcross/' + developer + '/')

    create_dir(save_fig_dir + 'Pkcross/' + developer + '/' )
    create_dir(save_fig_dir + 'fisher_Pkcross/' + developer + '/')

    obs_pars = lambda z_val: deepcopy(obs_params_ULTRASAT_CBR(z_val,Nbin))
    obs_pars_line = lambda z_val: deepcopy(obs_params_ULTRASAT((z_arr[0]+z_arr[-1])/2.))
    #obs_pars = lambda z_val: deepcopy(obs_params_SPHEREx(z_val))

    derive_pars = derive_model + derive_Pkline

    Fisher = np.zeros((len(derive_pars),len(derive_pars)))

    model_line = update_Pkline(obs_pars_line((z_arr[0]+z_arr[-1])/2.), astrocosmo_dict(developer,(z_arr[0]+z_arr[-1])/2.))[0]

    # tilde_Pkline_fid = Pk_AP(model_line,(z_arr[0]+z_arr[-1])/2., 1.,1.) + model_line.Pnoise

    for z in z_arr:

        model = update_Pkline(obs_pars(z), astrocosmo_dict(developer,z))[0]

        Fisher_z = np.zeros((len(derive_pars),len(derive_pars)))


        direct = save_fig_dir + 'Pkcross/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct += '_' + str(f_axion_fid) 

        direct += '/'
        direct_fisher = save_fig_dir + 'fisher_Pkcross/' + developer  + '/z' + str(z)
        if developer == 'axions':
            direct_fisher += '_' + str(f_axion_fid) 

        create_dir(direct)
        create_dir(direct_fisher)

        print('\n----------------------------------------------------\n---- Doing z = ' + str(z) + ' ----\n----------------------------------------------------\n')


        fid_par, name_fid = get_fiducials(derive_model, developer, z)

        vecpar, name_vec = get_fiducials_Pkline(derive_Pkline, model, z)

        if 'bLy' and 'bg' in derive_model:
            fiducials = [model.bavg[0],bg] + fid_par + vecpar
            names = [r'$b_{\rm Ly\alpha}$',r'$b_g$'] + name_fid + name_vec
        else:    
            fiducials = fid_par + vecpar
            names = name_fid + name_vec

        fiducial_k = model.k.value 

        fiducial_k = np.asarray(fiducial_k) 
        Pkline_fid = Pk_AP(model_line,z, 1.,1.) + model_line.Pnoise    

        integrand_cov = (Pkline_fid+model.Pnoise)/(model.Nmodes**0.5) 
        #integrand_cov = ((Pkcross_fid**2 + tilde_Pkline_fid*(Pkgal_fid+model.nbar**-1))/(2*model.Nmodes))**0.5 

        if 2 in multipole:
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
        else:
            cov = [[[]]]
            der_Pk_multi = [[]]

        for l in range(len(multipole)):

            print('--------- Doing L = ' + str(multipole[l]) + ' ---------')
            print('Running Pk fiducial')
            fiducial_Pk = 0.5*np.trapz(Pkline_fid, model.mu, axis=0) if multipole[l] == 0 else (2*multipole[l]+1)/2.*np.trapz(Pkline_fid*legendre(multipole[l])(model.mui_grid),model.mu,axis=0) 

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

                filename_low = direct + 'Pkcross_' + str(multipole[l])+ '_low_' + str(derive_model[i])			
                filename_up = direct + 'Pkcross_' + str(multipole[l]) + 'up_' + str(derive_model[i])

                model_low, model_up = set_der_par(derive_model[i],z,developer,obs_pars)

                if os.path.exists(filename_low) and import_allowed:
                    k_low, Pkcross_low = import_Pkline(filename_low, '\nDoing derivative wrt: ' + str(derive_model[i]) + '\nImport Pk low')
                else:
                    k_low = model_low.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_low = model.bavg[0]*(1-delta)
                        else:
                            bLy_low = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_low = bg*(1-delta)
                        else:
                            bg_low = bg
                        Pkclus_low = Pk_cross(model,z,bLy_low,bg_low, 1.,1.)

                    else:

                        Pkclus_low = Pk_cross(model_low,z,model.bavg[0],bg,1.,1.)

                    Pkcross_low = 0.5*np.trapz(Pkclus_low , model_low.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_low*legendre(multipole[l])(model_low.mui_grid),model_low.mu,axis=0)

                    if save_Pkline_flag:
                        np.savetxt(filename_low,(k_low,Pkcross_low))

                all_k_low.append(np.asarray(k_low)) 
                all_Pk_low.append(np.asarray(Pkcross_low))

                if os.path.exists(filename_up) and import_allowed:
                    k_up, Pkcross_up = import_Pkline(filename_up, 'Import Pk up')

                else:
                    k_up = model_up.k

                    if derive_model[i] == 'bLy' or derive_model[i] == 'bg':
                        if derive_model[i] == 'bLy':
                            bLy_up = model.bavg[0]*(1+delta)
                        else:
                            bLy_up = model.bavg[0]
                        if derive_model[i] == 'bg':
                            bg_up = bg*(1+delta)
                        else:
                            bg_up = bg

                        Pkclus_up = Pk_cross(model,z,bLy_up,bg_up, 1.,1.)

                    else:

                        Pkclus_up = Pk_cross(model_up,z,model.bavg[0],bg,1.,1.)

                    Pkcross_up = 0.5*np.trapz(Pkclus_up, model_up.mu, axis=0) if multipole[l] == 0 else\
                    (2*multipole[l]+1)/2.*np.trapz(Pkclus_up*legendre(multipole[l])(model_up.mui_grid),model_up.mu,axis=0) 

                    if save_Pkline_flag:
                        np.savetxt(filename_up, (k_up,Pkcross_up))

                all_k_up.append(np.asarray(k_up)) 
                all_Pk_up.append(np.asarray(Pkcross_up))


                pkprop_ders = analytical_derivatives(derive_Pkline, model,z,Pkline_fid,multipole[l])

            der_k = []
            der_Pk = []
            for i in range(len(derive_pars)):

                if derive_pars[i] in derive_model:
                    der_k.append([])
                    der_Pk.append([])
                    id_der = derive_model.index(derive_pars[i])

                    if derive_pars[i] == 'bLy':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*model.bavg[0])
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    elif derive_pars[i] == 'bg':
                        for t in range(len(fiducial_k)):
                            denom = (2*delta*bg)
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)
                    else:      
                        if 'bLy' and 'bg' in derive_model:
                            denom = (2*delta*fid_par[id_der-2]) if fid_par[id_der-2] != 0. else delta
                        else:
                            denom = (2*delta*fid_par[id_der]) if fid_par[id_der] != 0. else delta
                        for t in range(len(fiducial_k)):
                            der_Pk[i].append((all_Pk_up[id_der][t] - all_Pk_low[id_der][t]) / denom)

                elif derive_pars[i] in derive_Pkline:
                    id_der = derive_Pkline.index(derive_pars[i])
                    all_k_up.append(model.k)
                    der_Pk.append(pkprop_ders[id_der])

                der_Pk_multi[l] = der_Pk

        inv_cov = []
        multi_der_Pk_ordered = []
        if 2 in multipole:
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
        else:
            for kval in range(len(cov[0][0])):
                use_cov = np.asarray(cov[0][0][kval])
                inv_cov.append(use_cov**-1)
            for a in range(len(derive_pars)):
                multi_der_Pk_ordered.append([])
                for kval in range(len(der_Pk_multi[0][a])):
                    multi_der_Pk_ordered[a].append(der_Pk_multi[0][a][kval])        

        inv_cov = np.asarray(inv_cov)
        multi_der_Pk_ordered = np.asarray(multi_der_Pk_ordered)
        
        Fisher_k = []
        Fisher = np.zeros((len(der_Pk),len(der_Pk)))
        for t in range(len(fiducial_k)):
            Fisher_k.append(np.zeros((len(der_Pk),len(der_Pk))))		
            for i in range(len(der_Pk)):
                for j in range(len(der_Pk)):
                    if len(multipole) > 1:
                        Fisher_k[t][i,j] = np.linalg.multi_dot(\
                            [multi_der_Pk_ordered[i][t],inv_cov[t],\
                            multi_der_Pk_ordered[j][t]])
                    else:
                        Fisher_k[t][i,j] = multi_der_Pk_ordered[i][t]*\
                            inv_cov[t]*multi_der_Pk_ordered[j][t]

        Fisher_k = np.asarray(Fisher_k)
        for i in range(len(der_Pk)):
            for j in range(len(der_Pk)):
                for t in range(len(fiducial_k)):
                    Fisher_z[i,j] +=  Fisher_k[t][i,j]
        
        Fisher += Fisher_z
    
    if save_fisher_flag:
        filename = save_fig_dir + 'fisher_Pkcross/fisher_Lya' 
        save_fisher(filename, derive_pars,Fisher)

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

        print('\nError on w0 = ' + str(np.sqrt(inv_Fisher[0][0])))
        print('Error on wa = ' + str(np.sqrt(inv_Fisher[1][1]) ) + '\n')

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

        g.plot_2d([plot_dist], names, filled=True,\
            #legend_labels=[r'$\tilde{P}_{\rm Ly\alpha}(k)\,{\rm constraints}$'],legend_loc='upper right',\
            contour_colors=[kitsune], line_args=[{'ls':'--', 'color':kitsune}],
            markers=[{'x2':0}])

        plt.tight_layout()
        plt.savefig(save_fig_dir + 'ellipse_DE.pdf')
        plt.show()

    print('Forecasted errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot)
    print('Percentage errors on the parameters ' + str(derive_pars) + ': ')
    print(sigma_tot / fiducials * 100)

    return Fisher, sigma_tot