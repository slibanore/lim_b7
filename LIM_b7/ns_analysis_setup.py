# SL: last update 01/18/2023

from .ns_fiducial_pars import *
from .modelling.limlam_setup import lim
from .VID_setup import lim_VID
from .Pkline_setup import lim_Pkline

# write fiducial pars in the default file
def update_fiducials(file_pars,par_list,par_value):
    
    f = open(file_pars, 'r')
    list_of_lines = f.readlines()
    for i in range(len(par_list)):
        par_found = False
        for j in range(len(list_of_lines)):
            if par_list[i] in list_of_lines[j]:
                par_id = j
                par_found = True
                break
        if not par_found:
            print(par_list[i] + ' parameter is missing in the default file!') 
        if type(par_value) == int:
            par_value = float(par_value)   
        list_of_lines[par_id] = par_list[i] + ' = ' + str(par_value[i]) + '\n'

    f.close()
    f = open(file_pars, 'w')
    f.writelines(list_of_lines)
    f.close()
    
    return 


# update the parameter list 
def update_pars(start_params, what_to_update):

	obs_params = deepcopy(start_params)
	obs_params.update(what_to_update)
	
	model = lim(obs_params)

	return model, obs_params

# update the VID parameter list 
def update_VID(start_params, what_to_update):

	obs_params = update_pars(start_params, what_to_update)[1]

	model = lim_VID(obs_params)

	return model, obs_params

# update the Pkline parameter list 
def update_Pkline(start_params, what_to_update):

	obs_params = update_pars(start_params, what_to_update)[1]
	
	model = lim_Pkline(obs_params)

	return model, obs_params



# define array of fiducial values depending 
# on the input list of parameters
def get_fiducials(derive_pars, developer, z):

	fiducials = []
	names = []
	wrong_pars = False

	for i in range(len(derive_pars)):
		
		if derive_pars[i] == 'a1':
			if developer == 'CDS':
				fiducials.append(CDS_alpha_1_fid)
				names.append(r'$a_1$')
			else:
				wrong_pars = True
				
		if derive_pars[i] == 'a2':
			if developer == 'CDS':
				fiducials.append(CDS_alpha_2_fid)
				names.append(r'$a_2$')
			else:
				wrong_pars = True
				
		if derive_pars[i] == 'a3':
			if developer == 'CDS':
				fiducials.append(CDS_alpha_3_fid)
				names.append(r'$a_3$')
			else:
				wrong_pars = True
				
		if derive_pars[i] == 'a4':
			if developer == 'CDS':
				fiducials.append(CDS_alpha_4_fid)
				names.append(r'$a_4$')
			else:
				wrong_pars = True
				
		if derive_pars[i] == 'f_FDM':
			if developer == 'axions':
				fiducials.append(f_axion_fid)
				names.append(r'$f_{FDM}$')
			else:
				wrong_pars = True
				
		if derive_pars[i] == 'astro_alpha':
			fiducials.append(astro_alpha_fid)
			names.append(r'$\alpha_{astro}$')
			
		if derive_pars[i] == 'astro_beta':
			fiducials.append(astro_beta_fid)
			names.append(r'$\beta_{astro}$')
			
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

		if derive_pars[i] == 'fNL':
			fiducials.append(fNL_fid)
			names.append(r'$f_{NL}$')

		if derive_pars[i] == 'ns':
			fiducials.append(ns_fid)
			names.append(r'$n_s$')

		if derive_pars[i] == 'nrun':
			fiducials.append(nrun_fid)
			names.append(r'$\alpha$')

		if derive_pars[i] == 'nrunrun':
			fiducials.append(nrunrun_fid)
			names.append(r'$\beta$')

	if wrong_pars:
		print('\n-------------------------')
		print('Some parameter in the list are not allowed in the model you are using')
		print('-------------------------\n')
		return 
	
	return fiducials, names


# setup default survey based on input name
def set_survey(model_id, developer, estimator = 'VID'):

	if estimator == 'VID':
		update_and_run = lambda start_params, what_to_update: update_VID(start_params, what_to_update)
	elif estimator == 'Pkline':
		update_and_run = lambda start_params, what_to_update: update_Pkline(start_params, what_to_update)
	else:
		print('Check estimator!')
		return
	
	if model_id == 'COMAP1':
		survey_name = r'$COMAP$'
		z_vals = [2.9]
		obs_pars = lambda z: deepcopy(obs_params_COMAP1(z))

	elif model_id == 'IMS3':
		survey_name = r'$IMS3$'
		z_vals = [5.3]
		obs_pars = lambda z: deepcopy(obs_params_IMS3(z))

	elif model_id == 'EOR_lz':
		survey_name = r'${\rm COMAP-EoR}\, low\, z$'
		z_vals = [2.6, 3.]
		obs_pars = lambda z: deepcopy(obs_params_lowz(z,'EOR'))	

	elif model_id == 'EOR_hz':
		survey_name = r'${\rm COMAP-EoR}\, high,\ z$'
		z_vals = [5.3, 6.25, 7.25]
		obs_pars = lambda z: deepcopy(obs_params_highz(z,'EOR'))	

	elif model_id == 'EOR_lz_deep':
		survey_name = r'${\rm COMAP-EoR}\, low,\ z\, deep$'
		z_vals = [2.6, 3.]
		obs_pars = lambda z: deepcopy(obs_params_lowz(z,'deep'))	

	elif model_id == 'EOR_hz_deep':
		survey_name = r'${\rm COMAP-EoR}\, high\, z\, deep$'
		z_vals = [5.3, 6.25, 7.25]
		obs_pars = lambda z: deepcopy(obs_params_highz(z,'deep'))	

	elif model_id == 'COS3_lz':
		survey_name = r'${\rm COS3}\, low\, z$'
		z_vals = [2.6, 3]
		obs_pars = lambda z: deepcopy(obs_params_lowz(z,'COS3'))	

	elif model_id == 'COS3_lz_wide':
		survey_name = r'${\rm COS3}\, low\, z,\, wide$'
		z_vals = [2.6, 3]
		obs_pars = lambda z: deepcopy(obs_params_lowz(z,'wide'))	

	elif model_id == 'COS3_hz':
		survey_name = r'${\rm COS3}\, high\, z$'
		z_vals = [5.3, 6.25, 7.25]
		obs_pars = lambda z: deepcopy(obs_params_highz(z,'COS3'))	
	
	elif model_id == 'COS3_hz_wide':
		survey_name = r'${\rm COS3}\, high\, z,\, wide$'
		z_vals = [5.3, 6.25, 7.25]
		obs_pars = lambda z: deepcopy(obs_params_highz(z,'wide'))	
	
	elif model_id == 'COS3_hz_deep':
		survey_name = r'${\rm COS3}\, high\, z,\, deep$'
		z_vals = [5.3, 6.25, 7.25]
		obs_pars = lambda z: deepcopy(obs_params_highz(z,'deep'))	
	
	else:
		print('---- Using the model given as an input ----')
		model_id = 'INPUT'
		survey_name = ''

	if survey_name == '':
		model = model_id
	else:
		model = lambda z: update_and_run(obs_pars(z), astrocosmo_dict(developer,z))[0]

	if developer == 'axion':
		model_name = deepcopy(model_id) + '_fDM' + str(f_axion_fid) 
	else:
		model_name = deepcopy(model_id)

	return survey_name, model_name, obs_pars, model, z_vals


# set up derivatives for line power spectrum computation
def get_fiducials_Pkline(derive_pars, model, z):

	fiducials = []
	names = []

	for i in range(len(derive_pars)):
		
		if derive_pars[i] == 'Tfs8':

			fsigma = model.cosmo.get_fsigma8()[::-1]
			fsigma_interp = interp1d(model.zcosmo,fsigma)
			Tfs8 = model.Tmean*fsigma_interp(z)	

			fiducials.append(Tfs8)
			names.append(r'$Tf\sigma_8$')

		if derive_pars[i] == 'Tbs8':

			sigma = model.cosmo.get_sigma8()[::-1]
			sigma_interp = interp1d(model.zcosmo,sigma)
			Tbs8 = model.Tmean*model.bavg[0]*sigma_interp(z)	

			fiducials.append(Tbs8)
			names.append(r'$Tb\sigma_8$')

		if derive_pars[i] == 's8':

			sigma = model.cosmo.get_sigma8()[::-1]
			sigma_interp = interp1d(model.zcosmo,sigma)
			s8 = sigma_interp(z)	

			fiducials.append(s8)
			names.append(r'$\sigma_8$')

		if derive_pars[i] == 'fs8':

			sigma = model.cosmo.get_fsigma8()[::-1]
			fsigma_interp = interp1d(model.zcosmo,sigma)
			fs8 = fsigma_interp(z)	

			fiducials.append(fs8)
			names.append(r'$f\sigma_8$')

		if derive_pars[i] == 'Pshot':

			Pshot = model.L2mean.value*model.CLT.value**2.

			fiducials.append(Pshot)
			names.append(r'$P_{shot}$')

		if derive_pars[i] == 'fNL':

			fiducials.append(0.)
			names.append(r'$f_{NL}$')

		if derive_pars[i] == 'sNL':

			fiducials.append(model.sigma_NL.value)
			names.append(r'$s_{NL}$')

		if derive_pars[i] == 'alpha_par':

			fiducials.append(1.)
			names.append(r'$\alpha_{\parallel}$')

		if derive_pars[i] == 'alpha_par':

			fiducials.append(1.)
			names.append(r'$\alpha_{\perp}$')


	return fiducials, names

