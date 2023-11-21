# SL: last update 01/20/2023

from LIM_b7 import *
from LIM_b7.ultrasat_fiducial_pars import astrocosmo_dict

save_fig_dir = './results/ultrasat/decayDM/' 

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

data = lambda developer, z:{**astrocosmo_dict(developer,z), **model_data}

###########################################
# SURVEY PARAMETERS
###########################################

survey = 'ULTRASAT' 
redshift = (0.89+1.38)/2. 
detector_params = lambda z: obs_params_ULTRASAT(z)

fid_model = update_VID(\
	        detector_params(redshift),\
		    data(model_data['developer'],redshift))[0]

fid_model_pk = update_Pkline(\
	        detector_params(redshift),\
		    data(model_data['developer'],redshift))[0]



def plot_vid(save_figure = False):

	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

	print('\nComputing PT...')
	I_fid = fid_model.T	+ fid_model.Tmean
	PT_fid = fid_model.PT
		
	z = np.linspace(0.3, 1.1, 5)
	mass = lambda nu,z: ((4*np.pi*cu.hbar*(1+z)*nu.to(u.s**-1)).to(u.eV)).value
	color = ['#d0d1e6', '#a6bddb', '#74a9cf', '#2b8cbe', '#045a8d']

	##########################################
	# plot procedure
	##########################################

	plt.figure(figsize=(12,7))

	for i in range(len(z)):
		# this is the modified model 
		test_data = {**astrocosmo_dict('decayDM',z[i]), **model_data}

		test_data['developer'] = 'decayDM'
		test_data['Theta_DM'] = Theta_DM_fid
		
		test_model = update_VID(\
			detector_params(z[i]), test_data)[0]

		print('\nComputing DM PT...')
		PT_DM = test_model.PT
		print(len(PT_fid),len(PT_DM))

		f_Pfid = fft(PT_fid)#*fid_model.dT
		#f_Pfid = ((f_Pfid*f_Pfid.unit).decompose()).value 

		f_PDM = fft(PT_DM)#*fid_model.dT
		#f_PDM = ((f_PDM*f_PDM.unit).decompose()).value 

		fPT_tot = f_Pfid*f_PDM
					
		PT_tot_un = (ifft(fPT_tot)).real#/fid_model.dT).real
		PT_tot = PT_tot_un / np.trapz(PT_tot_un,fid_model.T)


		print('\n-------------------------------')

		##########################################
		# 2) plot PT
		##########################################

		plt.loglog(I_fid, PT_tot,color = color[i],linestyle='-',label=r'$m_{\rm DM}c^2 = %g$'%round(mass(fid_model.nuObs,z[i]),0) + r'${\rm eV},\, z_{\rm DM} = %g$'%z[i])

	##########################################
	# 2) plot PT
	##########################################

	plt.loglog(I_fid, PT_fid,color = kitsune,linestyle='-',label=r'$z_{\rm Ly\alpha} = %g$'%redshift)
		
	##########################################

	print('\nThe plot is on its way!')

	plt.legend(bbox_to_anchor=(1, 0.7),ncol=1)

	plt.ylabel(r'$\mathcal{P}(I)$')
	plt.xlabel(r'$I\, {\rm [Jy/sr]}$')
	plt.xlim(5e-1,4e2)
	plt.ylim(1e-6,1e-1)
	
	##########################################

	plt.tight_layout()
	if save_figure:		
		plt.savefig(save_fig_dir + 'vid_DM.pdf')

	plt.show()

	return 


def plot_pk(save_figure = False):

	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

	print('\nComputing monopole...')
	k_fid = fid_model_pk.k	
	pk0_fid = fid_model_pk.Pk_0
	
	print('\n...computing quadrupole...')
	pk2_fid = fid_model_pk.Pk_2

	print('\n...computing execapole...')
	pk4_fid = fid_model_pk.Pk_4

	z = np.linspace(0.3, 1.1, 5)
	mass = lambda nu,z: ((4*np.pi*cu.hbar*(1+z)*nu.to(u.s**-1)).to(u.eV)).value
	color = ['#d0d1e6', '#a6bddb', '#74a9cf', '#2b8cbe', '#045a8d']

	##########################################
	# plot procedure
	##########################################

	plt.figure(figsize=(20,6))

	outer = gridspec.GridSpec(\
		1, 3,left=0.08,bottom=0.1,right= 0.8, 
		top=0.98,wspace=.25,hspace=.25)

	ax1 = plt.subplot(outer[0])
	ax2 = plt.subplot(outer[1])
	ax3 = plt.subplot(outer[2])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)
	ax3.get_yaxis().set_label_coords(-0.12,0.5)

	for i in range(len(z)):
		# this is the modified model 
		test_data = {**astrocosmo_dict('decayDM',z[i]), **model_data}

		test_data['developer'] = 'decayDM'
		test_data['Theta_DM'] = Theta_DM_fid
		
		test_model = update_Pkline(\
			detector_params(z[i]), test_data)[0]

		print('\nComputing monopole...')
		k_test = test_model.k	
		pk0_test = test_model.Pk_0
		
		print('\n...computing quadrupole...')
		pk2_test = test_model.Pk_2

		print('\n...computing execapole...')
		pk4_test = test_model.Pk_4

		ax1.loglog(k_test, pk0_test,color = color[i],linestyle='-',)#label=r'$m_{\rm DM}c^2 = %g\,$'%round(mass(fid_model.nuObs,z[i]),0) + r'${\rm eV},\, z_{\rm DM} = %g$'%z[i])
		ax1.loglog(k_test, -pk0_test,color = color[i],linestyle='--')
		
		ax2.loglog(k_test, pk2_test,color = color[i],linestyle='-',)#label=r'$m_{\rm DM}c^2 = %g\,$'%round(mass(fid_model.nuObs,z[i]),0) + r'${\rm eV},\, z_{\rm DM} = %g$'%z[i])
		ax2.loglog(k_test, -pk2_test,color = color[i],linestyle='--')

		ax3.loglog(k_test, pk4_test,color = color[i],linestyle='-',label=r'$m_{\rm DM}c^2 = %g\,$'%round(mass(fid_model.nuObs,z[i]),0) + r'${\rm eV},\, z_{\rm DM} = %g$'%z[i])
		ax3.loglog(k_test, -pk4_test,color = color[i],linestyle='--')


	##########################################
	# 1) plot P0
	##########################################

	ax1.loglog(k_fid, pk0_fid,color = kitsune,linestyle='-',label=r'$z_{\rm Ly\alpha} = %g$'%redshift)
	
	ax1.set_ylabel(r'$\tilde{P}_{\rm Ly\alpha}^0(k,z=%g)$'%round(redshift,2))
	#ax1.legend(loc=1)

	ax1.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	ax1.set_xlim(1e-3,1e1)
	ax1.set_ylim(1e-1,1e6)

	##########################################
	# 2) plot P2
	##########################################

	ax2.loglog(k_fid,pk2_fid,color = kitsune,label=r'$z_{\rm Ly\alpha} = %g$'%redshift)
	ax2.loglog(k_fid,-pk2_fid,linestyle='--',color = kitsune)

	ax2.set_ylabel(r'$|\tilde{P}_{\rm Ly\alpha}^2(k,z=%g)|$'%round(redshift,2))
	
	#ax2.legend(loc=1)

	ax2.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
		
	ax2.set_xlim(1e-3,1e1)
	ax2.set_ylim(1e-1,1e6)
	
	##########################################
	# 3) plot P4
	##########################################

	ax3.loglog(k_fid,pk4_fid,color = kitsune,label=r'$z_{\rm Ly\alpha} = %g$'%redshift)
	ax3.loglog(k_fid,-pk4_fid,linestyle='--',color = kitsune)

	ax3.set_ylabel(r'$|\tilde{P}_{\rm Ly\alpha}^4(k,z=%g)|$'%round(redshift,2))

	ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax3.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	
	ax3.set_xlim(1e-3,1e1)
	ax3.set_ylim(1e-1,1e6)
	
	##########################################

	if save_figure:		
		plt.savefig(save_fig_dir + 'pk_dm_1e25.pdf')

	plt.show()
		
	return 
