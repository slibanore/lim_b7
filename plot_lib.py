# SL: last update 01/20/2023

from LIM_b7 import *
from LIM_b7.fiducial_pars import astrocosmo_dict


save_fig_dir = './results/standard_code/' 

###########################################
# MODEL PARAMETERS
###########################################

# dictionaries for minimal setup 
# use to check models under analysis
# see default_pars for the full parameter list
model_data = dict(\
       developer = 'CDS',               
       CDS_alpha_1 = CDS_alpha_1_fid, 
	    CDS_alpha_2 = CDS_alpha_2_fid, 
       CDS_alpha_3 = CDS_alpha_3_fid, 
	    CDS_alpha_4 = CDS_alpha_4_fid, 
       k_12 = k_12_fid, 
	    k_23 = k_23_fid, 
	    k_34 = k_34_fid )

# model_data = dict(\
#       developer = 'axions',
#       m_axion = m_axion_fid,
#       f_axion = f_axion_fid)

# model_data = dict(\
#         developer = 'PMF',
#         n_B = n_B_fid,         
#         sigma_B_0 = sigma_B_0_fid,    
#         smooth_scale = smooth_scale_fid)

data = lambda developer, z:{**astrocosmo_dict(developer,z), **model_data}

###########################################
# SURVEY PARAMETERS
###########################################

survey = 'COS3' # EOR, deep, wide
redshift = 6.25 # 2.6, 3, 5.3, 6.25, 7.25
detector_params = lambda z, s: obs_params_lowz(z, s) if z < 4 else obs_params_highz(z, s)

fid_model = update_VID(\
	        detector_params(redshift, survey),\
		    data(model_data['developer'],redshift))[0]
fid_model_pk = update_Pkline(\
	        detector_params(redshift, survey),\
		    data(model_data['developer'],redshift))[0]


###########################################
# COMPARE DIFFERENT MODELS
###########################################

# Input pars: parameter to change, value to use
# Output: plot Pm(k), dndM, PT, Bi
def compare_models(mod_par = 'CDS_alpha_3', 
                    mod_val = 0.5, save_figure = False):

	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

	if type(mod_val) == types.LambdaType:
		mod_val = mod_val(redshift)
		
	# this is the modified model 
	test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

	test_data[mod_par] = mod_val
	
	test_model = update_VID(\
		detector_params(redshift,survey), \
	    test_data)[0]

	##########################################
	# 1) compute the matter power spectrum
	##########################################
	
	print('\n-------------------------------')
	print('Computing matter power spectrum...')
	
	pk_fid = fid_model.Pm[0]
	k_fid = fid_model.k

	pk_test = test_model.Pm[0]
	k_test = test_model.k

	##########################################
	# 2) compute the halo mass function
	##########################################

	print('\n...computing halo mass function...')
	
	n_fid = fid_model.dndM
	M_fid = fid_model.M

	n_test = test_model.dndM
	M_test = test_model.M

	##########################################
	# 3) compute PT
	##########################################

	print('\n...computing fiducial PT...')
	pT_fid = fid_model.PT	
	T_fid = fid_model.T + fid_model.Tmean

	print('\n...computing modified PT...')
	pT_test = test_model.PT	
	T_test = test_model.T + test_model.Tmean
	
	##########################################
	# 4) compute VID per shell
	##########################################
	
	print('\n...computing fiducial VID...')

	Ti_fid, Bi_fid = fid_model.get_VID(\
		Tbin = False, T = fid_model.T, Tmean = fid_model.Tmean,\
	    PT = pT_fid,\
		PT_zero=fid_model.PT_zero,\
	    Nbin = Npt_i, minBi = min_Bi)

	print('\n...computing modified VID')
	Ti_test, Bi_test = test_model.get_VID(\
		Tbin = Ti_fid, T = test_model.T,\
		Tmean = test_model.Tmean,\
	    PT = pT_test, \
		PT_zero = test_model.PT_zero,\
	    Nbin = Npt_i, minBi = min_Bi)

	print('\n-------------------------------')
	print('\nThe plot is on its way!')

	##########################################
	# plot procedure
	##########################################

	plt.figure()

	outer = gridspec.GridSpec(\
		2, 2,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[0],
		wspace=0., hspace=0.)
	ax1 = plt.subplot(inner[0:2,0])
	ax1A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[1],
		wspace=0., hspace=0.)
	ax2 = plt.subplot(inner[0:2,0])
	ax2A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[2],
		wspace=0., hspace=0.)
	ax3 = plt.subplot(inner[0:2,0])
	ax3A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[3],
		wspace=0., hspace=0.)
	ax4 = plt.subplot(inner[0:2,0])
	ax4A = plt.subplot(inner[2,0])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax1A.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)
	ax2A.get_yaxis().set_label_coords(-0.12,0.5)
	ax3.get_yaxis().set_label_coords(-0.12,0.5)
	ax3A.get_yaxis().set_label_coords(-0.12,0.5)
	ax4.get_yaxis().set_label_coords(-0.12,0.5)
	ax4A.get_yaxis().set_label_coords(-0.12,0.5)
			

	##########################################
	# 1) plot the matter power spectrum
	##########################################

	fid_label = r'$\Lambda \rm CDM$'
	try:
		test_label = r'%s'%mod_par + r'$\, = \, %g$'%mod_val
	except:
		test_label = r'$Modified\, model$'
	ax1.loglog(k_fid, pk_fid, linestyle='--',label = fid_label, color = aonibi)
	ax1.loglog(k_test, pk_test, color = chojizome,linestyle='-', label=test_label)
	
	ax1.set_ylabel(r'$P(k,%g)\, {\rm [Mpc^{3}]}$'%redshift)
	ax1.legend(loc=3, ncol = 1,frameon=False)
	ax1.set_ylim(3e-5,3e4)
	ax1.set_xlim(k_fid[0].value-0.1,k_fid[-1].value+50)
	ax1.set_xticks([])

	ax1A.plot(k_fid, (pk_fid-pk_fid)/pk_fid,color = aonibi,linestyle='--')
	ax1A.plot(k_test, (pk_test-pk_fid)/pk_fid,linestyle='-',color=chojizome)

	ax1A.set_xscale('log')
	ax1A.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	ax1A.set_ylabel(r'$\frac{\Delta P}{P}$')
	ax1A.set_xlim(k_fid[0].value-0.1,k_fid[-1].value+50)
	ax1A.set_ylim(max(-1.5,min(2*(pk_test-pk_fid)/pk_fid)),min(1.5,max(2*(pk_test-pk_fid)/pk_fid)))

	##########################################
	# 2) plot the halo mass function
	##########################################

	ax2.loglog(M_fid, n_fid, color = aonibi,linestyle='--')
	ax2.loglog(M_test, n_test, linestyle='-',color=chojizome)

	ax2.set_ylim(1.5e-25,1e-8)
	ax2.set_ylabel(r'$\frac{dn}{dM_{h}}\, {\rm [Mpc^{-3}}\,M_{\odot} {\rm ]}$')
	ax2.set_xticks([])

	ax2A.plot(M_fid, (n_fid-n_fid)/n_fid, color = aonibi, linestyle='--')
	ax2A.plot(M_test, (n_test-n_fid)/n_fid, linestyle='-', color = chojizome)

	ax2A.set_xscale('log')	
	ax2A.set_xlabel(r'$M_h\, {\rm [}M_{\odot}{\rm ]}$')
	ax2A.set_ylabel(r'$\frac{\Delta {\rm hmf}}{{\rm hmf}}$')
	ax2A.set_ylim(max(-1.5,min(2*(n_test-n_fid)/n_fid)),min(1.5,max(2*(n_test-n_fid)/n_fid)))

	##########################################
	# 3) plot PT
	##########################################

	ax3.loglog(T_fid, pT_fid,color = aonibi,linestyle='--')
	ax3.loglog(T_test, pT_test, linestyle = '-',color=chojizome)
	
	ax3.set_ylabel(r'$\mathcal{P}(T)$')
	ax3.set_xticks([])
	ax3.set_xlim(5e-5,100)
	ax3.set_ylim(8e-7,5e1)

	ax3A.plot(T_fid,(pT_fid-pT_fid)/pT_fid,linestyle='--',color=aonibi)
	ax3A.plot(T_test,(pT_test-pT_fid)/pT_fid,linestyle='-',color=chojizome)
	
	ax3A.set_ylabel(r'$\frac{\Delta \mathcal{P}(T)}{\mathcal{P}(T)}$')
	ax3A.set_xlabel(r'$T\, {\rm [\mu K]}$')
	ax3A.set_xlim(5e-5,100)
	ax3A.set_ylim(max(-1.5,min(2*(pT_test-pT_fid)/pT_fid)),min(1.5,max(2*(pT_test-pT_fid)/pT_fid)))
	ax3A.set_xscale('log')

	##########################################
	# 4) plot VID
	##########################################

	ax4.loglog(Ti_fid,Bi_fid,color = aonibi,linestyle='--')
	ax4.loglog(Ti_test,Bi_test,linestyle='-',color=chojizome)

	ax4.set_ylabel(r'$B_i(z\sim %g)$'%round(redshift,0))
	ax4.set_xticks([])
	ax4.set_xlim(1e0,9e1)
	ax4.set_ylim(5e-4,3e2)

	ax4A.plot(Ti_fid,(Bi_fid-Bi_fid)/Bi_fid,linestyle='--',color=aonibi)
	ax4A.plot(Ti_test,(Bi_test-Bi_fid)/Bi_fid,linestyle='-',color=chojizome)
	
	ax4A.set_xscale('log')	
	ax4A.set_xlim(1e0,9e1)
	ax4A.set_ylim(max(-1.5,min(2*(Bi_test-Bi_fid)/Bi_fid)),min(1.5,max(2*(Bi_test-Bi_fid)/Bi_fid)))
	ax4A.set_ylabel(r'$\frac{\Delta B_i}{B_i}$')	
	ax4A.set_xlabel(r'$T_i\, {\rm [\mu K]}$')
		
	##########################################

	if save_figure:		
		plt.savefig(save_fig_dir + 'compare_' + mod_par + '.pdf')

	plt.show()

	return Ti_fid, Bi_fid, Ti_test, Bi_test



def compare_multi_vals(mod_par, mod_list, save_figure, get_SNR = False):
	
	##########################################
	# plot procedure
	##########################################

	plt.figure()

	outer = gridspec.GridSpec(\
		2, 2,left=0.08,bottom=0.1,right= 0.98, 
		#1, 2,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[0],
		wspace=0., hspace=0.)
	ax1 = plt.subplot(inner[0:2,0])
	ax1A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[1],
		wspace=0., hspace=0.)
	ax2 = plt.subplot(inner[0:2,0])
	ax2A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[2],
		wspace=0., hspace=0.)
	ax3 = plt.subplot(inner[0:2,0])
	ax3A = plt.subplot(inner[2,0])

	inner = gridspec.GridSpecFromSubplotSpec(\
		3, 1, subplot_spec=outer[3],
		wspace=0., hspace=0.)
	ax4 = plt.subplot(inner[0:2,0])
	ax4A = plt.subplot(inner[2,0])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax1A.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)
	ax2A.get_yaxis().set_label_coords(-0.12,0.5)
	ax3.get_yaxis().set_label_coords(-0.12,0.5)
	ax3A.get_yaxis().set_label_coords(-0.12,0.5)
	ax4.get_yaxis().set_label_coords(-0.12,0.5)
	ax4A.get_yaxis().set_label_coords(-0.12,0.5)

	fid_label = r'$\Lambda \rm CDM$'
	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

	print('\nComputing fiducial matter power spectrum...')

	pk_fid = fid_model.Pm[0]
	k_fid = fid_model.k

	print('\n...computing fiducial halo mass function...')

	n_fid = fid_model.dndM
	M_fid = fid_model.M

	print('\n...computing fiducial PT...')
	pT_fid = fid_model.PT	
	T_fid = fid_model.T + fid_model.Tmean

	print('\n...computing fiducial VID...')

	Ti_fid, Bi_fid = fid_model.get_VID(\
		Tbin = False, T = fid_model.T, Tmean = fid_model.Tmean,\
		PT = pT_fid,\
		PT_zero=fid_model.PT_zero,\
		Nbin = Npt_i, minBi = min_Bi)

	ax1.loglog(k_fid, pk_fid, linestyle='--',label = fid_label, color = aonibi)

	ax1A.plot(k_fid, (pk_fid-pk_fid)/pk_fid,color = aonibi,linestyle='--')

	ax2.loglog(M_fid, n_fid, color = aonibi,linestyle='--')

	ax2A.plot(M_fid, (n_fid-n_fid)/n_fid, color = aonibi, linestyle='--')
#  
	ax3.loglog(T_fid, pT_fid,color = aonibi,linestyle='--')
	
	ax3A.plot(T_fid,(pT_fid-pT_fid)/pT_fid,linestyle='--',color=aonibi)

	ax4.loglog(Ti_fid,Bi_fid,color = aonibi,linestyle='--')

	ax4A.plot(Ti_fid,(Bi_fid-Bi_fid)/Bi_fid,linestyle='--',color=aonibi)
	
	use_id = 0
	for mod_val in mod_list:

		use_color = colors[use_id]

		if type(mod_val) == types.LambdaType:
			mod_val = mod_val(redshift)

		if mod_par == 'ns':
			cosmo = dict(
				f_NL=0, H0=67.67, cosmomc_theta=None,
				ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
				num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
				YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
				TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
				theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
				dark_energy_model='ppf',As=2.105e-09, 
				# !!! Planck 2018 eq. 16/17/18 arXiv:1807.06211 
				#ns=mod_val, nrun=nrun_fid, nrunrun=nrunrun_fid, 
				r=0.0, nt=None, ntrun=0.0, 
				pivot_scalar=0.05, pivot_tensor=0.05,
				parameterization=2,halofit_version='mead')


			test_data = {**astrocosmo_dict(model_data['developer'],redshift), **dict(developer = model_data['developer'], cosmo_input_camb = cosmo)}

		elif mod_par == 'nrun':
			cosmo = dict(
				f_NL=0, H0=67.67, cosmomc_theta=None,
				ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
				num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
				YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
				TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
				theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
				dark_energy_model='ppf',As=2.105e-09, 
				# !!! Planck 2018 eq. 16/17/18 arXiv:1807.06211 
				#ns=ns_fid, nrun=mod_val, nrunrun=nrunrun_fid, 
				r=0.0, nt=None, ntrun=0.0, 
				pivot_scalar=0.05, pivot_tensor=0.05,
				parameterization=2,halofit_version='mead')

			test_data = {**astrocosmo_dict(model_data['developer'],redshift), **dict(developer = model_data['developer'], cosmo_input_camb = cosmo,)}
			
		elif mod_par == 'nrunrun':

			cosmo = dict(
				f_NL=0, H0=67.67, cosmomc_theta=None,
				ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
				num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
				YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
				TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
				theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
				dark_energy_model='ppf',As=2.105e-09, 
				# !!! Planck 2018 eq. 16/17/18 arXiv:1807.06211 
				ns=ns_fid, nrun=nrun_fid, nrunrun=mod_val, 
				r=0.0, nt=None, ntrun=0.0, 
				pivot_scalar=0.05, pivot_tensor=0.05,
				parameterization=2,halofit_version='mead')

			test_data = {**astrocosmo_dict(model_data['developer'],redshift), **dict(developer = model_data['developer'], cosmo_input_camb = cosmo,)}
			
		# this is the modified model 
		else:
			test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

			test_data[mod_par] = mod_val

		test_model = update_VID(\
			detector_params(redshift,survey), \
			test_data)[0]

		##########################################
		# 1) compute the matter power spectrum
		##########################################

		print('\n-------------------------------')
		print('\n...computing modified matter power spectrum (id = ' + str(use_id) + ')...')

		pk_test = test_model.Pm[0]
		k_test = test_model.k

		##########################################
		# 2) compute the halo mass function
		##########################################

		print('\n...computing modified halo mass function (id = ' + str(use_id) + ')...')

		n_test = test_model.dndM
		M_test = test_model.M

		##########################################
		# 3) compute PT
		##########################################

		print('\n...computing modified PT (id = ' + str(use_id) + ')...')
		pT_test = test_model.PT	
		T_test = test_model.T + test_model.Tmean

		##########################################
		# 4) compute VID per shell
		##########################################

		print('\n...computing modified VID (id = ' + str(use_id) + ')...')
		Ti_test, Bi_test = test_model.get_VID(\
			Tbin = Ti_fid, T = test_model.T,\
			Tmean = test_model.Tmean,\
			PT = pT_test, \
			PT_zero = test_model.PT_zero,\
			Nbin = Npt_i, minBi = min_Bi)

		if get_SNR:
			N_z_shells = int(fid_model.Delta_nu.value / (fid_model.dnu.value * 1e-3))
			SNR = np.sqrt(N_z_shells*sum((Bi_test-Bi_fid)**2/Bi_fid))        

			print('SNR = ' + str(SNR))

		##########################################
		# 1) plot the matter power spectrum
		##########################################
		try:
			test_label = r'%s'%mod_par + r'$\, = \, %g$'%mod_val
		except:
			test_label = r'$Model\, %g$'%use_id

		ax1.loglog(k_test, pk_test, color = use_color,linestyle='-', label=test_label)
		
		ax1A.plot(k_test, (pk_test-pk_fid)/pk_fid,linestyle='-',color=use_color)


		##########################################
		# 2) plot the halo mass function
		##########################################

		ax2.loglog(M_test, n_test, linestyle='-',label=r'$%s\, = $'%mod_par + r'$%g$'%mod_val)#,color=use_color)

	
		ax2A.plot(M_test, (n_test-n_fid)/n_fid, linestyle='-')#, color = use_color)


		##########################################
		# 3) plot PT
		##########################################

		ax3.loglog(T_test, pT_test, linestyle = '-',color=use_color)
		
		ax3A.plot(T_test,(pT_test-pT_fid)/pT_fid,linestyle='-',color=use_color)
		# 
		##########################################
		# 4) plot VID
		##########################################

		ax4.loglog(Ti_test,Bi_test,linestyle='-')#,color=use_color)

		ax4.set_ylabel(r'$B_i(z\sim %g)$'%round(redshift,0))
		ax4.set_xticks([])
		ax4.set_xlim(1e0,9e1)
		ax4.set_ylim(5e-4,3e2)

		ax4A.plot(Ti_test,(Bi_test-Bi_fid)/Bi_fid,linestyle='-')#,color=use_color)
		
		use_id += 1

		##########################################

	print('\n-------------------------------')
	print('\nThe plot is on its way!')

	ax1.set_ylabel(r'$P(k,%g)\, {\rm [Mpc^{3}]}$'%redshift)
	ax1.legend(loc=3, ncol = int(len(mod_list)/2),frameon=False)
	ax1.set_ylim(3e-5,3e4)
	ax1.set_xlim(k_fid[0].value-0.1,k_fid[-1].value+50)
	ax1.set_xticks([])

	ax1A.set_xscale('log')
	ax1A.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	ax1A.set_ylabel(r'$\frac{\Delta P}{P}$')
	ax1A.set_xlim(k_fid[0].value-0.1,k_fid[-1].value+50)
	ax1A.set_ylim(-1.5,1.5)

	ax2.set_ylim(1.5e-25,1e-8)
	ax2.set_ylabel(r'$\frac{dn}{dM_{h}}\, {\rm [Mpc^{-3}}\,M_{\odot} {\rm ]}$')
	ax2.set_xticks([])

	ax2A.set_xscale('log')	
	ax2A.set_xlabel(r'$M_h\, {\rm [}M_{\odot}{\rm ]}$')
	ax2A.set_ylabel(r'$\frac{\Delta {\rm hmf}}{{\rm hmf}}$')
	ax2A.set_ylim(-.5,2)#-1.5,1.5)

	ax3.set_ylabel(r'$\mathcal{P}(T)$')
	ax3.set_xticks([])
	ax3.set_xlim(5e-5,100)
	ax3.set_ylim(8e-7,5e1)

	ax3A.set_ylabel(r'$\frac{\Delta \mathcal{P}(T)}{\mathcal{P}(T)}$')
	ax3A.set_xlabel(r'$T\, {\rm [\mu K]}$')
	ax3A.set_xlim(5e-5,100)
	#ax3A.set_ylim(-1.5,1.5)
	ax3A.set_yscale('log')#set_ylim(-0.1,0.1)
	ax3A.set_xscale('log')

	ax4A.set_xscale('log')	
	ax4A.set_xlim(1e0,9e1)
	#a44A.set_ylim(-1.5,1.5)
	ax4A.set_yscale('log')#set_ylim(-0.1,0.1)
	ax4A.set_ylabel(r'$\frac{\Delta B_i}{B_i}$')	
	ax4A.set_xlabel(r'$T_i\, {\rm [\mu K]}$')
		

	if save_figure:		 
		plt.savefig(save_fig_dir + 'multi_compare_' + mod_par + '.pdf')

	plt.show()

	return 


def compare_scaledependence(mod_par, mod_list, save_figure = False):
		
	##########################################
	# 1) compute PT
	##########################################

	print('\nComputing monopole and bias...')
	k_fid = fid_model_pk.k	
	pk0_fid = fid_model_pk.Pk_0
	bias_fid = fid_model_pk.bavg


	##########################################
	# plot procedure
	##########################################

	plt.figure(figsize=(15,8))

	outer = gridspec.GridSpec(\
		1, 2,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	ax1 = plt.subplot(outer[0])
	ax2 = plt.subplot(outer[1])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)

	##########################################
	# 1) plot P0
	##########################################

	ax1.loglog(k_fid, pk0_fid,color = aonibi,linestyle='--',label=r'$\Lambda {\rm CDM}$')

	ax2.loglog(k_fid, bias_fid,color = aonibi,linestyle='--',label=r'$\Lambda {\rm CDM}$')

	for mod_val in mod_list:

		if type(mod_val) == types.LambdaType:
			mod_val = mod_val(redshift)
			
		# this is the modified model 
		test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

		test_data[mod_par] = mod_val

		test_model = update_Pkline(\
			detector_params(redshift,survey), \
			test_data)[0]

		print('\n...computing modified monopole and bias...')

		k_test = test_model.k	
		pk0_test = test_model.Pk_0
		b_test = test_model.bavg

		ax1.loglog(k_test, pk0_test)
		ax2.loglog(k_test, b_test)
	
	print('\n-------------------------------')
	print('\nThe plot is on its way!')

	
	ax1.set_ylabel(r'$\tilde{P}_{\rm CO}(k,z=%g))$'%round(redshift,2))
	ax1.legend(loc=1)

	ax1.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	# ax1.set_xlim(1e-3,1e1)
	ax1.set_ylim(1e-12,1e4)

	ax2.legend(loc=1)

	ax2.set_ylabel(r'$b_{\rm CO}(k,z=%g))$'%round(redshift,2))
	ax2.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
		
	# ax2.set_xlim(1e-3,1e1)
	# ax2.set_ylim(1e-1,1e5)
	
	##########################################

	if save_figure:		
		plt.savefig(save_fig_dir + 'pk_fnl.pdf')

	plt.show()

	return 
