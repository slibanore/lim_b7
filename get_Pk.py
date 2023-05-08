from LIM_b7 import *
from LIM_b7.fiducial_pars import astrocosmo_dict


save_dir = './results/nu_mass/' 
save_fig_dir = './PLOTS/nu_mass/' 

def run_pk(get_matter_pk,
           get_line_pk,
           model_data,
           survey,
           redshift,
           save_flag,
           plot_flag):

    if save_flag:
        create_dir('./results')
        create_dir('./PLOTS')
        create_dir(save_dir)
        create_dir(save_fig_dir)

    if plot_flag:
        plt.figure(figsize=(10,5))

    data = lambda developer, z:{**astrocosmo_dict(developer,z), **model_data}

    detector_params = lambda z, s: obs_params_lowz(z, s) if z < 4 else obs_params_highz(z, s)

    fid_model = update_Pkline(\
                detector_params(redshift, survey),\
                data(model_data['developer'],redshift))[0]

    k_fid = fid_model.k	

    if get_matter_pk:

        pk_m = fid_model.Pm[0]

        if save_flag:
            np.savetxt(save_dir + 'matter_pk',(k_fid.value,pk_m.value),header='k [Mpc-1], Pk [Mpc3]')

        if plot_flag:

            plt.loglog(k_fid, pk_m,
                       label=r'$\rm Dark Matter$')

    if get_line_pk:

        pk_line_clust = fid_model.Pk_clust
        pk_line = fid_model.Pk
        pk_line_monopole = fid_model.Pk_0

        if save_flag:
            np.savetxt(save_dir + 'line_pk',np.array((k_fid,pk_line_monopole)),header='k [Mpc-1], Pk_line_monopole [Mpc3]')

        if plot_flag:

            plt.loglog(k_fid,       
                       pk_line_monopole,
                       label=r'$\rm Line\, monopole$')
    
    if plot_flag:

        plt.xlabel(r'$k\,[{\rm Mpc^{-1}}]$',fontsize=20)
        plt.ylabel(r'$P(k,z =%g)\,[{\rm Mpc^{3}}]$'%round(redshift,2),fontsize=20)
        plt.legend(loc=3,fontsize=20)

        plt.ylim(1e-10,1e4) 

        plt.tight_layout()
        plt.savefig(save_dir+'pk.png')
        plt.show()


    return