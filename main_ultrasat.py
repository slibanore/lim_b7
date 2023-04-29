# SL: last update 01/28/2023

import ultrasat_plot_lib as vp 
import ultrasat_VID_fisher_lib as vf
import ultrasat_Pkline_fisher as pf
import ultrasat_Pkline_CBR as cbr
from LIM_b7.ultrasat_fiducial_pars import *


def main():

    ###########################################
    # CHOOSE MODULE
    ###########################################

    # 2) choose which module to run 

    VID_plot = False
    Pk_plot = False
    VID_fisher = False
    Pk_fisher = True
    Pk_CBR = False
    Pk_cross_CBR = False
    Pk_gal_CBR = False

    ###########################################
    # SETUP ANALYSIS
    ###########################################

    # 3) set up specific analysis parameters
    # NOTE that each module has inside more setting that can be tuned 

    if VID_plot:

        save_figure = False
        vp.plot_vid(save_figure)

    if Pk_plot:

        save_figure = True
        vp.plot_pk(save_figure)



    if VID_fisher:

        z = (0.89+1.38)/2. 
        developer = 'CDS'
        derive_pars = ['fNL','A_lya','B_lya','D_lya','astro_sig_scatter','astro_Lcut']
        #'z0_Lyalpha','f0_Lyalpha','t_Lyalpha', 'astro_sig_scatter','astro_Lcut']
        save_VID_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True

        vf.fisher_VID(z, 
            developer, 
            derive_pars,
            save_VID_flag, 
            save_fisher_flag, 
            plot_sigma_flag, 
            import_allowed)

        #vf.error_dndL(z, developer)
        vf.error_dndL_new(z, developer)
        
    if Pk_fisher:

        z = (0.89+1.38)/2. 
        developer = 'CDS'
        derive_model = ['fNL','A_lya','B_lya','D_lya']#,'B_lya','D_lya',] #['z0_Lyalpha','f0_Lyalpha','t_Lyalpha']
        derive_Pkline = ['Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0,2,4]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True

        pf.fisher_Pkline(z, 
                developer, 
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed)

    if Pk_CBR:

        z = (0.89+1.38)/2. 
        developer = 'CDS'
        derive_model = ['w0','wa','A_lya','B_lya','D_lya']#,'z0_Lyalpha','f0_Lyalpha','t_Lyalpha']
        derive_Pkline = ['Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True

        cbr.fisher_Pkline(z, 
                developer, 
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed)

    if Pk_cross_CBR:

        zmin = 0.9
        zmax = 1.4
        dz = 0.1
        Nbin = int(round((zmax-zmin)/dz,0))
        z_arr = np.zeros(Nbin)
        zi = zmin 
        for i in range(Nbin):
            z_arr[i] = zi 
            zi += dz

        bg = 2

        developer = 'CDS'
        derive_model = ['w0','wa','A_lya','B_lya','D_lya']#,'z0_Lyalpha','f0_Lyalpha','t_Lyalpha']
        derive_Pkline = ['Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True

        cbr.fisher_Pkcross(z_arr, Nbin,
                developer, bg,
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed)

    if Pk_gal_CBR:

        zmin = 0.9
        zmax = 1.4
        dz = 0.1
        Nbin = int(round((zmax-zmin)/dz,0))
        z_arr = np.zeros(Nbin)
        zi = zmin 
        for i in range(Nbin):
            z_arr[i] = zi 
            zi += dz

        bg = 2

        developer = 'CDS'
        derive_model = ['w0','wa']
        derive_Pkline = ['Pshot','sNL','alpha_par','alpha_perp'] 
        multipole = [0]
        save_Pkline_flag = True
        save_fisher_flag = True
        plot_sigma_flag = True
        import_allowed = True

        cbr.fisher_Pkgal(z_arr, Nbin,
                developer, bg,
                derive_model,
                derive_Pkline,	
                multipole,
                save_Pkline_flag, \
                save_fisher_flag, \
                plot_sigma_flag, \
                import_allowed)


if __name__ == "__main__":
    main()