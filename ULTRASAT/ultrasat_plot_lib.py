# SL: last update 01/20/2023

from LIM_b7 import *
from LIM_b7.ultrasat_fiducial_pars import astrocosmo_dict


save_fig_dir = './results_ultrasat/' 

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


# Input pars: parameter to change, value to use
# Output: plot Pm(k), dndM, PT, Bi
def plot_vid(save_figure = False):

	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

		
	##########################################
	# 1) compute PT
	##########################################

	print('\nComputing fiducial PT...')
	pT_fid = fid_model.PT	
	T_fid = fid_model.T + fid_model.Tmean
	
	##########################################
	# 2) compute VID per shell
	##########################################
	
	print('\n...computing fiducial VID...')

	Ti_fid, Bi_fid = fid_model.get_VID(\
		Tbin = False, T = fid_model.T, Tmean = fid_model.Tmean,\
	    PT = pT_fid,\
		PT_zero=fid_model.PT_zero,\
	    Nbin = Npt_i, minBi = min_Bi)

	# print('... different models ...')
	# # this is the modified model 
	# test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

	# test_data['model_par'] = dict(A_lya=A_Lyalpha_fid/3.,B_lya=B_Lyalpha_fid, D_lya=D_Lyalpha_fid, dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)
	
	# test_model = update_VID(\
	# 	detector_params(redshift), \
	#     test_data)[0]

	# pT_A = test_model.PT	
	# T_A = test_model.T + test_model.Tmean
	
	# Ti_A, Bi_A = test_model.get_VID(\
	# 	Tbin = Ti_fid, T = test_model.T, Tmean = test_model.Tmean,\
	#     PT = pT_A,\
	# 	PT_zero=test_model.PT_zero,\
	#     Nbin = Npt_i, minBi = min_Bi)

	# test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

	# test_data['model_par'] = dict(A_lya=A_Lyalpha_fid,B_lya=B_Lyalpha_fid/3., D_lya=D_Lyalpha_fid, dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)
	
	# test_model = update_VID(\
	# 	detector_params(redshift), \
	#     test_data)[0]

	# pT_B = test_model.PT	
	# T_B = test_model.T + test_model.Tmean
	
	# Ti_B, Bi_B = test_model.get_VID(\
	# 	Tbin = Ti_fid, T = test_model.T, Tmean = test_model.Tmean,\
	#     PT = pT_B,\
	# 	PT_zero=test_model.PT_zero,\
	#     Nbin = Npt_i, minBi = min_Bi)

	# test_data = {**astrocosmo_dict(model_data['developer'],redshift), **model_data}

	# test_data['model_par'] = dict(A_lya=A_Lyalpha_fid,B_lya=B_Lyalpha_fid, D_lya=D_Lyalpha_fid/3., dndL_Lcut=astro_Lcut_fid, SFR_file=SFR_file_fid)
	
	# test_model = update_VID(\
	# 	detector_params(redshift), \
	#     test_data)[0]

	# pT_D = test_model.PT	
	# T_D = test_model.T + test_model.Tmean
	
	# Ti_D, Bi_D = test_model.get_VID(\
	# 	Tbin = Ti_fid, T = test_model.T, Tmean = test_model.Tmean,\
	#     PT = pT_D,\
	# 	PT_zero=test_model.PT_zero,\
	#     Nbin = Npt_i, minBi = min_Bi)


	print('\n-------------------------------')
	print('\nThe plot is on its way!')

	##########################################
	# plot procedure
	##########################################

	plt.figure(figsize=(15,7))

	outer = gridspec.GridSpec(\
		1, 2,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	ax1 = plt.subplot(outer[0])
	ax2 = plt.subplot(outer[1])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)

	##########################################
	# 1) plot PT
	##########################################

	ax1.loglog(T_fid, pT_fid,color = kitsune,linestyle='-',label=r'${\rm Ly}\alpha$')
	
	# ax1.loglog(T_fid, pT_A,color = suoko,linestyle='-',label=r'${\mathcal{A}/3}$')
	
	# ax1.loglog(T_fid, pT_B,color = suoko,linestyle='--',label=r'${\mathcal{B}/3}$')
	
	# ax1.loglog(T_fid, pT_D,color = seiheki,linestyle='-',label=r'${\mathcal{D}/3}$')
	
	ax1.set_ylabel(r'$\mathcal{P}(I,z= %g)$'%round(redshift,2))
	ax1.set_xticks([])
	ax1.legend(loc=1)

	ax1.set_xlabel(r'$I\, {\rm [Jy/sr]}$')
	ax1.set_xscale('log')
	ax1.set_xlim(3e-1,2e2)
	ax1.set_ylim(1e-6,1e-1)

	##########################################
	# 2) plot VID
	##########################################

	ax2.loglog(Ti_fid,Bi_fid,color = kitsune,label=r'${\rm Ly}\alpha$')

	# ax2.loglog(Ti_fid, Bi_A,color = suoko,linestyle='-',label=r'${\mathcal{A}/3}$')
	
	# ax2.loglog(Ti_fid, Bi_B,color = suoko,linestyle='--',label=r'${\mathcal{B}/3}$')
	
	# ax2.loglog(Ti_fid, Bi_D,color = seiheki,linestyle='-',label=r'${\mathcal{D}/3}$')

	ax2.vlines(fid_model.sigma_N.value,1e6,2e7,color='gray',linestyle=':')

	ax2.set_ylabel(r'$B_i(z= %g)$'%round(redshift,2))
	ax2.set_xticks([])
	
	ax2.legend(loc=1)

	ax2.set_xscale('log')	
	ax2.set_xlabel(r'$I_i\, {\rm [Jy/sr]}$')
		
	ax2.text(10,1e7,r'${\rm Noise\, dominated}$',fontsize=20)

	ax2.set_xlim(2e0,2e2)
	ax2.set_ylim(1e6,2e7)
	
	##########################################

	if save_figure:		
		plt.savefig(save_fig_dir + 'VID.pdf')

	plt.show()

	return Ti_fid, pT_fid, Bi_fid


def plot_pk(save_figure = False):

	if save_figure and not os.path.exists(save_fig_dir): 
		os.makedirs(save_fig_dir,exist_ok=True)

		
	##########################################
	# 1) compute PT
	##########################################

	print('\nComputing monopole...')
	k_fid = fid_model_pk.k	
	pk0_fid = fid_model_pk.Pk_0
	
	##########################################
	# 2) compute VID per shell 
	##########################################
	
	print('\n...computing quadrupole...')

	pk2_fid = fid_model_pk.Pk_2

	print('\n...computing execapole...')

	pk4_fid = fid_model_pk.Pk_4

	print('\n-------------------------------')
	print('\nThe plot is on its way!')

	##########################################
	# plot procedure
	##########################################

	plt.figure(figsize=(20,6))

	outer = gridspec.GridSpec(\
		1, 3,left=0.08,bottom=0.1,right= 0.98, 
		top=0.98,wspace=.25,hspace=.25)

	ax1 = plt.subplot(outer[0])
	ax2 = plt.subplot(outer[1])
	ax3 = plt.subplot(outer[2])

	ax1.get_yaxis().set_label_coords(-0.12,0.5)
	ax2.get_yaxis().set_label_coords(-0.12,0.5)
	ax3.get_yaxis().set_label_coords(-0.12,0.5)

	##########################################
	# 1) plot P0
	##########################################

	ax1.loglog(k_fid, pk0_fid,color = kitsune,linestyle='-',label=r'${\rm Ly}\alpha$')
	
	ax1.set_ylabel(r'$\tilde{P}_{\rm Ly\alpha}^0(k,z=%g)$'%round(redshift,2))
	ax1.legend(loc=1)

	ax1.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	ax1.set_xlim(1e-3,1e1)
	ax1.set_ylim(1e-1,1e5)

	##########################################
	# 2) plot P2
	##########################################

	ax2.loglog(k_fid,pk2_fid,color = kitsune,label=r'${\rm Ly}\alpha$')
	ax2.loglog(k_fid,-pk2_fid,linestyle='--',color = kitsune)

	ax2.set_ylabel(r'$|\tilde{P}_{\rm Ly\alpha}^2(k,z=%g)|$'%round(redshift,2))
	
	ax2.legend(loc=1)

	ax2.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
		
	ax2.set_xlim(1e-3,1e1)
	ax2.set_ylim(1e-1,1e5)
	
	##########################################
	# 3) plot P4
	##########################################

	ax3.loglog(k_fid,pk4_fid,color = kitsune,label=r'${\rm Ly}\alpha$')
	ax3.loglog(k_fid,-pk4_fid,linestyle='--',color = kitsune)

	ax3.set_ylabel(r'$|\tilde{P}_{\rm Ly\alpha}^4(k,z=%g)|$'%round(redshift,2))

	ax3.legend(loc=1)
	ax3.set_xlabel(r'$k\, {\rm [Mpc^{-1}]}$')
	
	ax3.set_xlim(1e-3,1e1)
	ax3.set_ylim(1e-1,1e5)
	
	##########################################

	if save_figure:		
		plt.savefig(save_fig_dir + 'pk.pdf')

	plt.show()

	return k_fid, pk0_fid, pk2_fid

