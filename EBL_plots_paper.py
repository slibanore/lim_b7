from EBL_forecast import * 


def plot_collective():

    fig, ax = plt.subplots(4,1, figsize=(9,19))

    fontsize_label = 23
    fontsize_legend = 21
    fontsize_tick = 19
###########################
    wavelenght, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T

    #plt.subplot(311)
    z = [0.,0.5,1.,2]

    rest_wave = np.linspace(700,3500,400)
    for zv in range(len(z)):
        s = np.zeros(len(rest_wave))
        wave = np.zeros((len(rest_wave)))
        for i in range(len(rest_wave)):
            wave[i] = rest_wave[i] * (1.+z[zv])
            if z[zv] <= 1:
                detector = 'GALEX_NUV' # same for galex_fuv in EW
            else:
                detector = 'ULTRASAT' # same for hetdex and spherex
            s[i] = signal(rest_wave[i]*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value 
        s_1500 = signal(1500*u.AA,z[zv],detector,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,to_plot=True).value 
        if z[zv] == 1.:
            ax[0].plot(wave, s, color= 'k')
            ax[0].plot(1500*(1+z[zv]), s_1500, marker = 'o',color= 'k', zorder= len(z)-zv,markersize=7)
        ax[1].plot(wave, s, color= colors[zv], label=r'$z = %g$'%z[zv], zorder= len(z)-zv)
        ax[1].plot(1500*(1+z[zv]), s_1500, marker = 'o',color= colors[zv], zorder= len(z)-zv,markersize=7)

    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_legend)
    ax[0].set_xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize_label)
    ax[0].set_xlim(700,4900)
    ax[0].set_ylim(5e24,1e28)

    ax[0].text(1500*(1+1), 3e26, r'$\alpha_{1500}(z)$', rotation=-1.5, fontsize=fontsize_legend, ha='center', va='center')
    ax[0].annotate('', xy=(0.55, 0.45), xytext=(0.55, 0.01),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='<->'), ha='center', va='center')
    ax[0].text(1570*(1+1), 3e25, r'$\epsilon_{1500}(z)$', rotation=90, fontsize=fontsize_legend, ha='center', va='center')

    ax[0].text(1030*(1+1), 1.8e26, r'$\alpha_{1100}(z)$', rotation=25, fontsize=fontsize_legend, ha='center', va='center')

    ax[0].text(750*(1+1), 2.6e25, r'$\alpha_{900}(z)$', rotation=18, fontsize=fontsize_legend, ha='center', va='center')

    ax[0].text(990*(1+1), 2.5e25, r'$f_{\rm LyC}(z)$', rotation=90, fontsize=fontsize_legend, ha='center', va='center')

    ax[0].text(1216*(1+1), 5e27, r'$\rm EW(z)$', rotation=0, fontsize=fontsize_legend, ha='center', va='center')


    ax[0].set_xticks([]) 
    ax[0].set_yticks([]) 
    ax[0].set_title(r'$\rm Comoving\,Volume\,Emissivity$',fontsize=fontsize_label,y=1.01)
#    ax[0].set_title(r'$\rm Comoving\,volume\,emissivity$',fontsize=fontsize_label, y=.8, x=.71)

    ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize_label)
    ax[1].set_ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_label)
    ax[1].legend(loc=2,ncol=2,fontsize=fontsize_legend)
    ax[1].set_xlim(700,4900)
    ax[1].set_ylim(1e25,5e28)

    ax[1].set_title(r'$\rm UV\,-\,EBL$',fontsize=fontsize_label,y=1.01)

    ax[1].tick_params(axis='y', labelsize=fontsize_tick) 
    ax[1].tick_params(axis='x', labelsize=fontsize_tick) 
###########################

    #plt.subplot(312)
    wavelenght_F = np.linspace(wavelenght_min('GALEX_FUV'),wavelenght_max('GALEX_FUV'))
    wavelenght_N = np.linspace(wavelenght_min('GALEX_NUV'),wavelenght_max('GALEX_NUV'))
    wavelenght = np.linspace(wavelenght_min('ULTRASAT'),wavelenght_max('ULTRASAT'))

    R = np.zeros(len(wavelenght))
    RN = np.zeros(len(wavelenght))
    RF = np.zeros(len(wavelenght))

    for i in range(len(wavelenght)):
        R[i] = Response(wavelenght[i],'ULTRASAT')
        RN[i] = Response(wavelenght_N[i],'GALEX_NUV')
        RF[i] = Response(wavelenght_F[i],'GALEX_FUV')

    ax[2].plot(wavelenght_F,RF-min(RF),label=r'$\rm GALEX\, FUV$',color=color_FUV,zorder=6)
    ax[2].fill_between(wavelenght_F.value, 0, RF - min(RF), color=color_FUV, alpha=0.3,zorder=5)    
    ax[2].plot(wavelenght_N,RN-min(RN),label=r'$\rm GALEX\, NUV$',color=color_NUV,zorder=2)
    ax[2].fill_between(wavelenght_N.value, 0, RN - min(RN), color=color_NUV, alpha=0.3,zorder=1)
    ax[2].plot(wavelenght,R-min(R),label=r'$\rm ULTRASAT$',color=color_ULTRASAT,zorder=4)
    ax[2].fill_between(wavelenght.value, 0, R - min(R), color=color_ULTRASAT, alpha=0.3,zorder=3)

    ax[2].set_yticks([])
    ax[2].tick_params(axis='x', labelsize=fontsize_tick) 

    ax[2].set_xlim(700,4900)
    ax[2].set_ylim(0,5.5)

    ax[2].set_title(r'$\rm Broadband\,Observation$',fontsize=fontsize_label,y=1.01)

    ax[2].legend(loc=1,ncol=1,fontsize=fontsize_legend)
    ax[2].set_xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize_label)
    ax[2].set_ylabel(r'$R(\lambda_{\rm obs})$',fontsize=fontsize_label)
########################

    dJ_U =  dJdz(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename='results/EBL/dJdz_ULTRASAT_reduced.dat')

    dJ_N =  dJdz(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename='results/EBL/dJdz_GALEX_NUV_reduced.dat')

    dJ_F =  dJdz(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename='results/EBL/dJdz_GALEX_FUV_reduced.dat')

    bJ_U = bJ_z(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_ULTRASAT_reduced.dat')

    bJ_N = bJ_z(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_GALEX_NUV_reduced.dat')

    bJ_F = bJ_z(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_GALEX_FUV_reduced.dat')


    window_size = 70  # Adjust this value based on your preference
    bJU_smoothed = moving_average(dJ_U*bJ_U, window_size)
    bJGN_smoothed = moving_average(dJ_N*bJ_N, window_size)
    bJGF_smoothed = moving_average(dJ_F*bJ_F, window_size)

    #plt.subplot(313)
    ax[3].plot(z_gals('SDSS')[:len(bJGF_smoothed)],bJGF_smoothed,label=r'$\rm GALEX\, FUV$',color=color_FUV)
    ax[3].plot(z_gals('SDSS')[:len(bJGN_smoothed)],bJGN_smoothed,label=r'$\rm GALEX\, NUV$',color=color_NUV)
    ax[3].plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_smoothed,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    ax[3].set_xlabel(r'$z$',fontsize=fontsize_label)
    ax[3].set_ylabel(r'$b_JdJ_{\nu_{\rm obs}}/dz\,[{\rm Jy/sr}]$',fontsize=fontsize_label)
    ax[3].legend(fontsize=fontsize_legend,loc=1)

    ax[3].tick_params(axis='y', labelsize=fontsize_tick) 
    ax[3].tick_params(axis='x', labelsize=fontsize_tick) 

    ax[3].set_title(r'$\rm CBR\,Reconstruction$',fontsize=fontsize_label,y=1.01)

    ax[3].set_xlim(0,2.3)
###############################


    ax[0].xaxis.set_label_coords(0.5, -.06)
    ax[0].yaxis.set_label_coords(-.07, 0.5)
    ax[1].yaxis.set_label_coords(-.07, 0.5)
    ax[2].yaxis.set_label_coords(-.07, 0.5)
    ax[3].yaxis.set_label_coords(-.07, 0.5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.37)
    plt.savefig('results/PLOTS/EBL/collective_plot.png',bbox_inches='tight')

    plt.show()

    return 



def plot_dNgdz():

    plt.figure(figsize = (9, 8))
    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    z_sdss = np.linspace(zmin_gal,zmax_gal,80) 
    z_desi = np.linspace(zmin_gal,zmax_gal,80) 
    N = np.zeros(len(z_sdss))
    Ndesi =  np.zeros(len(z_desi))

    for i in range(len(z_sdss)):
        N[i] = dNgdz(z_sdss[i],'SDSS')*delta_zi_interp#delta_zi('SDSS')
    for i in range(len(z_desi)):
        Ndesi[i] = dNgdz(z_desi[i],'DESI')*delta_zi_interp #delta_zi('DESI') 
        
    print(np.trapz(N/delta_zi('SDSS'),z_sdss))
    print(np.trapz(Ndesi/delta_zi('DESI'),z_desi))

    window_size = 7  # Adjust this value based on your preference
    Nsdss_smoothed = moving_average(N, window_size)
    Ndesi_smoothed = moving_average(Ndesi, window_size)

    plt.plot(z_desi[:len(Ndesi_smoothed)], Ndesi_smoothed,color=color_ULTRASAT,label=r'$\rm 5yr\,DESI$',)
    plt.plot(z_sdss[:len(Nsdss_smoothed)],Nsdss_smoothed,label=r'$\rm SDSS$',color=color_FUV)

    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.ylabel(r'$N_{g,i}\,(\Delta z_i = %g)$'%delta_zi_interp,fontsize=fontsize_label)
    plt.legend(loc=4,fontsize=fontsize_legend)
    plt.tight_layout()

    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    plt.xlim(min(z_gals_interp),min(max(z_gals_interp),z_desi[:len(Ndesi_smoothed)][-1]))

    plt.savefig('results/PLOTS/EBL/Ng.png',bbox_inches='tight')
    plt.show()

    return 


def plot_noise_grouped():

    plt.figure(figsize = (9, 8))
    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    groups = np.concatenate((np.linspace(0.001,1,50),np.linspace(1,50,50)))
    sigma_GAL_NUV = np.zeros(len(groups))
    sigma_GAL_FUV = np.zeros(len(groups))
    sigma_GAL_NUV_D = np.zeros(len(groups))
    sigma_GAL_FUV_D = np.zeros(len(groups))
    sigma_ULT = np.zeros(len(groups))
    
    for i in tqdm(range(len(groups))):
        sigma_GAL_NUV[i] = sigma_wz(0.5,'GALEX_NUV','SDSS',groups[i]).value
        sigma_GAL_FUV[i] = sigma_wz(0.5,'GALEX_FUV','SDSS',groups[i]).value
        sigma_GAL_NUV_D[i] = sigma_wz(0.5,'GALEX_NUV','DESI',groups[i]).value
        sigma_GAL_FUV_D[i] = sigma_wz(0.5,'GALEX_FUV','DESI',groups[i]).value
        sigma_ULT[i] = sigma_wz(1,'ULTRASAT','DESI',groups[i]).value

    plt.plot(groups,sigma_GAL_FUV,color_FUV,label=r'${\rm GALEX\, FUV}\times{\rm SDSS}$',linestyle='--')
    plt.plot([5],[sigma_wz(0.5,'GALEX_FUV','SDSS',False).value],color_FUV,marker='o',markersize=10)
    plt.plot([50],[sigma_wz(0.5,'GALEX_FUV','SDSS',True).value],color_FUV,marker='x',markersize=20)
    plt.plot(groups,sigma_GAL_FUV_D,color_FUV,label=r'${\rm GALEX\, FUV}\times{\rm DESI}$')
    plt.plot([5],[sigma_wz(0.5,'GALEX_FUV','DESI',False).value],color_FUV,marker='o',markersize=10)
    plt.plot([50],[sigma_wz(0.5,'GALEX_FUV','DESI',True).value],color_FUV,marker='x',markersize=20)

    plt.plot(groups,sigma_GAL_NUV,color_NUV,label=r'${\rm GALEX\, NUV}\times{\rm SDSS}$',linestyle='--')
    plt.plot(groups,sigma_GAL_NUV_D,color_NUV,label=r'${\rm GALEX\, NUV}\times{\rm DESI}$')
    plt.plot([5],[sigma_wz(0.5,'GALEX_NUV','SDSS',False).value],color_NUV,marker='o',markersize=10)
    plt.plot([50],[sigma_wz(0.5,'GALEX_NUV','SDSS',True).value],color_NUV,marker='x',markersize=20)
    plt.plot([5],[sigma_wz(0.5,'GALEX_NUV','DESI',False).value],color_NUV,marker='o',markersize=10)
    plt.plot([50],[sigma_wz(0.5,'GALEX_NUV','DESI',True).value],color_NUV,marker='x',markersize=20)

    plt.plot(groups,sigma_ULT,color_ULTRASAT,label=r'${\rm ULTRASAT}\times{\rm DESI}$')
    plt.plot([5.4],[sigma_wz(1,'ULTRASAT','DESI',False).value],color_ULTRASAT,marker='o',markersize=10)
    plt.plot([5.4],[sigma_wz(1,'ULTRASAT','DESI',True).value],color_ULTRASAT,marker='x',markersize=20)

    plt.yscale('log')
    plt.xlabel(r'$L_{\rm pix}\,[{\rm arcsec/pix}]$',fontsize=fontsize_label)
    plt.ylabel(r'$\mathcal{N}_{\rm CBR}$',fontsize=fontsize_label)
    plt.legend(loc=1,ncol=1,fontsize=fontsize_legend)

    filename = 'results/PLOTS/EBL/noise_voxgroup.png'

    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 


def plot_noise():

    plt.figure(figsize = (9, 8))
    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    use_z = z_gals_interp
    use_z_ultrasat_desi = z_gals_interp
    reduced_label = '_reduced'

    sn = np.zeros(len(use_z))
    sf = np.zeros(len(use_z))
    snD = np.zeros(len(use_z))
    sfD = np.zeros(len(use_z))
    sUD = np.zeros(len(use_z_ultrasat_desi))
    nG = np.zeros(len(use_z))
    nGD = np.zeros(len(use_z))
    nUD = np.zeros(len(use_z_ultrasat_desi))

    filename_ULTDESI = 'results/EBL/wJg_ULTRASAT,DESI' + reduced_label + '.dat'
        
    for i in tqdm(range(len(use_z))):
        sn[i] = abs(wJgz(use_z[i],'GALEX_NUV','SDSS',False,filename='results/EBL/wJg_GALEX_NUV,SDSS' + reduced_label + '.dat'))
        sf[i] = abs(wJgz(use_z[i],'GALEX_FUV','SDSS',False,filename='results/EBL/wJg_GALEX_FUV,SDSS'+ reduced_label+'.dat'))
        nG[i] = sigma_wz(use_z[i],'GALEX_NUV','SDSS',True).value
        snD[i] = abs(wJgz(use_z[i],'GALEX_NUV','DESI',False,filename='results/EBL/wJg_GALEX_NUV,DESI' + reduced_label + '.dat'))
        sfD[i] = abs(wJgz(use_z[i],'GALEX_FUV','DESI',False,filename='results/EBL/wJg_GALEX_FUV,DESI'+ reduced_label+'.dat'))
        nGD[i] = sigma_wz(use_z[i],'GALEX_NUV','DESI',True).value
    for i in tqdm(range(len(use_z_ultrasat_desi))):
        sUD[i] = abs(wJgz(use_z_ultrasat_desi[i],'ULTRASAT','DESI',False,filename=filename_ULTDESI))
        nUD[i] = sigma_wz(use_z_ultrasat_desi[i],'ULTRASAT','DESI',True).value

    window_size = 2  # Adjust this value based on your preference
    U_smoothed = moving_average(sUD/nUD, window_size)
    GN_smoothed = moving_average(sn/nG, window_size)
    GF_smoothed = moving_average(sf/nG ,window_size)
    GND_smoothed = moving_average(snD/nGD, window_size)
    GFD_smoothed = moving_average(sfD/nGD ,window_size)
    #

    plt.plot(use_z[:len(GF_smoothed)],GF_smoothed,label=r'$\rm GALEX\,FUV\times SDSS$',color=color_FUV,linestyle='--')
    plt.plot(use_z[:len(GFD_smoothed)],GFD_smoothed,label=r'$\rm GALEX\,FUV\times DESI$',color=color_FUV)
    plt.plot(use_z[:len(GN_smoothed)],GN_smoothed,label=r'$\rm GALEX\,NUV\times SDSS$',color=color_NUV,linestyle='--')
    plt.plot(use_z[:len(GND_smoothed)],GND_smoothed,label=r'$\rm GALEX\,NUV\times DESI$',color=color_NUV)

    plt.plot(use_z_ultrasat_desi[:len(U_smoothed)],U_smoothed,label=r'$\rm ULTRASAT\times DESI$',color=color_ULTRASAT)
    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.ylabel(r'${\bar{\omega}_{\tilde{J}{\rm g}}(z)}/\mathcal{N}_{\rm CBR}$',fontsize=fontsize_label)
    plt.legend(loc=1,ncol=1,fontsize=fontsize_legend)

    plt.ylim(1,1e6)
    plt.xlim(use_z[0],min(use_z[-1],use_z_ultrasat_desi[:len(U_smoothed)][-1]))
    plt.tight_layout()

    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    filename = 'results/PLOTS/EBL/wJg_noise' + reduced_label + '.png'
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 


def plot_multi_line_PAPER():

    plt.figure(figsize = (11, 10))
    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    z = np.linspace(zmin_gal,zmax_gal+.5,200)

    sigmas_nonion, fid_nonion = plot_err_noninonizing_cont_wb(z, 1500, use_pars_fid = pars_original_c18_fid,group_vox=True,run=False,galex_detector='SDSS',prior='conservative')

    sigmas_nonion_both = plot_err_noninonizing_cont_wb(z, 1500, use_pars_fid = pars_fid,prior = 'conservative',group_vox=True,run=False,galex_detector='DESI')[0]
    sigmas_nonion_both_opt = plot_err_noninonizing_cont_wb(z, 1500, use_pars_fid = pars_fid,prior = 'optimistic',group_vox=True,run=False,galex_detector='DESI')[0]

    s_g = sigmas_nonion[0] * (70/67.67)**3
    s_b = sigmas_nonion_both[2] * (70/67.67)**3
    s_b_opt = sigmas_nonion_both_opt[2] * (70/67.67)**3


    z_data = np.asarray([0.3,0.5,0.7,1,2,2.9,3]) # schimonovich 2005

    logrho_1500_data = np.asarray([25.86,25.97,26.16,26.11,26.45,26.52,25.58])

    log_error_up = np.asarray([0.05,0.15,0.31,0.31,0.25,0.17,0.31])

    log_error_down = np.asarray([0.05,0.08,0.13,0.13,0.09,0.07,0.17])

    eps_1500_data = pow(10,logrho_1500_data)

    eps_err_up = pow(10,log_error_up)

    eps_err_down = pow(10,log_error_down)


    z_data_hst = (np.asarray([0.92,1.62,2.08]) + np.asarray([1.33,1.88,2.37]))/2.

    log_rho_hst = np.asarray([26.19,26.46,36.34])

    log_error_hst = np.asarray([0.08,0.12,0.09])

    eps_hst = pow(10,log_rho_hst)

    eps_err_hst = pow(10,log_error_hst)


    z_d11, eps_d11 = np.genfromtxt('./dat/d11.dat').T
    z_gil, eps_gil = np.genfromtxt('./dat/gilmore.dat').T
    z_gil2, eps_gil2 = np.genfromtxt('./dat/gilmore_mod.dat').T

    
    plt.errorbar([(1.9+2.7)/2], [np.log10(3.63e26)], yerr=[0.4], fmt='*', capsize=5, markersize=10, color='k', label=r'$\rm Reddy\, et\, al.,\,2008$',alpha=.6,markerfacecolor='k')
    
    plt.errorbar(z_data, np.log10(eps_1500_data), yerr=[np.log10(eps_err_down), np.log10(eps_err_up)], fmt='o', capsize=5, markersize=13, color='k', label=r'$\rm Schiminovich\, et\, al.,\,2005$',alpha=.6,markerfacecolor='k')
    
    plt.errorbar(z_data_hst, np.log10(eps_hst), yerr=np.log10(eps_err_hst), fmt='D', capsize=5, markersize=13, color='k', label=r'$\rm Dahlen\, et\, al.,\,2006$',alpha=.6,markerfacecolor='none',)

    plt.fill_between(z, np.log10(fid_nonion - s_b_opt), np.log10(fid_nonion + s_b_opt), color=color_ULTRASAT, alpha = 0.9,label=r'$\rm Optimistic\,bias\,prior$')

    plt.fill_between(z, np.log10(fid_nonion - s_b), np.log10(fid_nonion + s_b), color=color_ULTRASAT, alpha = 0.2,label=r'$\rm Conservative\,bias\,prior$')

    plt.plot(z, np.log10(fid_nonion - s_g),color_FUV,alpha=0.6,label=r'$\rm C18\,(reproduced)$',linestyle=(0,(5,1)),linewidth=1.8)
    
    plt.plot(z, np.log10(fid_nonion + s_g),color_FUV,alpha=0.6,linestyle=(0,(5,1)),linewidth=1.8)
    
    plt.plot(z_d11, np.log10(eps_d11), 'g', label=r'$\rm Dominguez\, et\, al.,\,2011$',alpha=.6,linewidth=1.4)
    
    plt.plot(z_gil, np.log10(eps_gil), 'k:', label=r'$\rm Gilmore\, et\, al.,\,2011$',linewidth=1.8)
    
    plt.plot(z_gil2, np.log10(eps_gil2), 'k:',linewidth=1.8)
    
    plt.plot(z, np.log10(fid_nonion),'k',zorder=2)

    print(np.max(s_b/fid_nonion))
    print(np.max(s_b_opt/fid_nonion))

    plt.ylabel(r'$\log_{10}\epsilon_{1500}(z)[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_label)
    plt.legend(loc=2,ncol=2,fontsize=fontsize_legend)
    plt.ylim(25.7,27.2)

    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.xlim(z[0],z[-1])

    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    filename = 'results/PLOTS/EBL/eps_err_multiwave_PAPER.png'
 
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()


    return 

pars_original_c18_fid = ['eps_1500_0', 'bias_1500_0','gamma_1500', 'alpha_1500_0', 'C_alpha_1500', 'alpha_1100_0', 'C_alpha_1100', 'EW_z1', 'EW_z2', 'gamma_bv', 'gamma_bz']

def plot_err_emissivity_PAPER():

    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    z = [0.,0.5,1.,2.]
    wave = np.linspace(700,3000,500)

    for zi in z:
        signal_val = np.zeros(len(wave))

        sigma_galex = np.zeros(len(wave))
        sigma_both = np.zeros(len(wave))
        sigma_both_opt = np.zeros(len(wave))
        for w in range(len(wave)):

            signal_val[w] =  signal(wave[w]*u.AA,zi,'ULTRASAT',vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False).value

            if wave[w] <= 912:
                sigma_galex[w] = plot_err_ioncont_wb([zi],wave[w], use_pars_fid = pars_original_c18_fid,group_vox=True,run=False, prior='conservative',galex_detector='SDSS')[0][0][0]
                
                sigma_both[w] = plot_err_ioncont_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False, prior='conservative',galex_detector='DESI')[0][2][0]

                sigma_both_opt[w] = plot_err_ioncont_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False, prior='optimistic',galex_detector='DESI')[0][2][0]

            elif 912 < wave[w] < 1216*(1-0.005):
                sigma_galex[w] = plot_err_line_surrounding_wb([zi],wave[w], use_pars_fid = pars_original_c18_fid,group_vox=True,run=False,prior='conservative',galex_detector='SDSS')[0][0][0]

                sigma_both[w] = plot_err_line_surrounding_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='conservative',galex_detector='DESI')[0][2][0]

                sigma_both_opt[w] = plot_err_line_surrounding_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='optimistic',galex_detector='DESI')[0][2][0]

            elif 1216*(1-0.005) < wave[w] <= 1216:
                sigma_galex[w] = plot_err_line_wb([zi],wave[w], use_pars_fid = pars_original_c18_fid,group_vox=True,run=False,prior='conservative',galex_detector='SDSS')[0][0][0]
    
                sigma_both[w] = plot_err_line_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='conservative',galex_detector='DESI')[0][2][0]
    
                sigma_both_opt[w] = plot_err_line_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='optimistic',galex_detector='DESI')[0][2][0]
    
            else:
                sigma_galex[w] = plot_err_noninonizing_cont_wb([zi],wave[w], use_pars_fid = pars_original_c18_fid,group_vox=True,run=False,prior='conservative',galex_detector='SDSS')[0][0][0]
                
                sigma_both[w] = plot_err_noninonizing_cont_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='conservative',galex_detector='DESI')[0][2][0]

                sigma_both_opt[w] = plot_err_noninonizing_cont_wb([zi],wave[w], use_pars_fid = pars_fid,group_vox=True,run=False,prior='optimistic',galex_detector='DESI')[0][2][0]


        plt.figure(figsize = (10, 7))

        plt.axvline(lambda_from_nu(nu_min_gNUV).value/(1+zi), color='k',linewidth=1.5,linestyle=':')

        plt.fill_betweenx([1e25,5e28],lambda_from_nu(nu_min_US).value/(1+zi),3000,color='k',alpha=0.1)
        plt.fill_betweenx([1e25,5e28],lambda_from_nu(nu_max_gFUV).value/(1+zi),color='k',alpha=0.1)

        plt.plot(wave, signal_val,'k',zorder=10)
        plt.fill_between(wave, signal_val - sigma_both_opt, signal_val + sigma_both_opt, color=color_ULTRASAT, alpha = 0.9, label=r'$\rm Optimistic\,bias\,prior$')
        plt.fill_between(wave, signal_val - sigma_both, signal_val + sigma_both, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm Conservative\,bias\,prior$')
        plt.plot(wave, signal_val - sigma_galex,color_FUV,alpha=0.6,label=r'$\rm C18\,(reproduced)$',linestyle=(0,(5,1)),linewidth=1.8)
        plt.plot(wave, signal_val + sigma_galex,color_FUV,alpha=0.6,linestyle=(0,(5,1)),linewidth=1.8)

        plt.ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize_label)
        plt.yscale('log')

        plt.tick_params(axis='y', labelsize=fontsize_tick) 
        plt.tick_params(axis='x', labelsize=fontsize_tick) 

        if zi == 2.:
            plt.legend(loc=1,fontsize=fontsize_legend)
            plt.ylim(1e25,3e28)
        elif zi == 1. or zi == 1.5:
            plt.legend(loc=1,fontsize=fontsize_legend)
            plt.ylim(1e25,3e28)
        else:
            plt.legend(loc=1,)
            plt.ylim(1e25,1e27)

        plt.xlim(wave[0],wave[-1])

        plt.title(r'$z=%g$'%zi,fontsize=fontsize_label,y=1.01)
        plt.xlabel(r'$\lambda_{\rm obs}/(1+z)\,[{\rm \AA}]$',fontsize=fontsize_label)

        filename = 'results/PLOTS/EBL/err_final_emissivity_z' + str(zi) + '.png'
        
        plt.tight_layout()
        plt.savefig(filename,bbox_inches='tight')
        plt.show()


    return 


def plot_err_escape(group_vox=True,run=False,plot_flag = True):

    plt.figure(figsize = (11, 10))
    fontsize_label = 27
    fontsize_legend = 24
    fontsize_tick = 23

    z = np.linspace(zmin_gal,zmax_gal,200)

    required_pars = ['log_fLyC_1','log_fLyC_2']

    use_pars_fid = pars_fid
    F_N_c18 = Fisher_change_var(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F_c18 = Fisher_change_var(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV','DESI',group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV','DESI',group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)
    F_N_opt = Fisher_change_var(use_pars_fid,'GALEX_NUV','DESI',group_vox,run)
    F_F_opt = Fisher_change_var(use_pars_fid,'GALEX_FUV','DESI',group_vox,run)
    F_U_opt = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    sigma_b1500_opt = 0.01
    sigma_gammabnu_opt = 0.3
    sigma_gammabz_opt = 0.1

    sigma_b1500 = 0.05
    sigma_gammabnu = 1.3
    sigma_gammabz = 0.3

    for j in range(len(use_pars_fid)):
        if use_pars_fid[j] == 'C_alpha_1500':
            F_N_c18[j,j] += 1/1.5**2
            F_F_c18[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_1100':
            F_N_c18[j,j] += 1/1.5**2
            F_F_c18[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'gamma_bv':
            F_N_c18[j,j] += 1/sigma_gammabnu**2
            F_F_c18[j,j] += 1/sigma_gammabnu**2
            F_N[j,j] += 1/sigma_gammabnu**2
            F_F[j,j] += 1/sigma_gammabnu**2
            F_U[j,j] += 1/sigma_gammabnu**2
        if use_pars_fid[j] == 'gamma_bz':
            F_N_c18[j,j] += 1/sigma_gammabz**2
            F_F_c18[j,j] += 1/sigma_gammabz**2
            F_N[j,j] += 1/sigma_gammabz**2
            F_F[j,j] += 1/sigma_gammabz**2
            F_U[j,j] += 1/sigma_gammabz**2
        if use_pars_fid[j] == 'bias_1500_0':
            F_N_c18[j,j] += 1/sigma_b1500**2
            F_F_c18[j,j] += 1/sigma_b1500**2
            F_N[j,j] += 1/sigma_b1500**2
            F_F[j,j] += 1/sigma_b1500**2
            F_U[j,j] += 1/sigma_b1500**2

        if use_pars_fid[j] == 'gamma_bv':
            F_N_opt[j,j] += 1/sigma_gammabnu_opt**2
            F_F_opt[j,j] += 1/sigma_gammabnu_opt**2
            F_U_opt[j,j] += 1/sigma_gammabnu_opt**2
        if use_pars_fid[j] == 'gamma_bz':
            F_N_opt[j,j] += 1/sigma_gammabz_opt**2
            F_F_opt[j,j] += 1/sigma_gammabz_opt**2
            F_U_opt[j,j] += 1/sigma_gammabz_opt**2
        if use_pars_fid[j] == 'bias_1500_0':
            F_N_opt[j,j] += 1/sigma_b1500_opt**2
            F_F_opt[j,j] += 1/sigma_b1500_opt**2
            F_U_opt[j,j] += 1/sigma_b1500_opt**2


    F_galex = F_N_c18 + F_F_c18
    F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)
    F_both_opt = sum_galex_ultrasat_desi(use_pars_fid, F_N_opt, F_F_opt, F_U_opt)

    all_inv_F_G = np.linalg.inv(F_galex)
    all_inv_F_b = np.linalg.inv(F_both)
    all_inv_F_b_opt = np.linalg.inv(F_both_opt)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b_opt = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
       id_i = use_pars_fid.index(required_pars[i])
       for j in range(len(required_pars)):
           id_j = use_pars_fid.index(required_pars[j])
           inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
           inv_F_b[i,j] = all_inv_F_b[id_i,id_j]
           inv_F_b_opt[i,j] = all_inv_F_b_opt[id_i,id_j]
#
    fid_fLy = np.zeros(len(z))

    sigma_eps_G = np.zeros(len(z))
    sigma_eps_b = np.zeros(len(z))
    sigma_eps_b_opt = np.zeros(len(z))

    for i in range(len(z)):

        fid_fLy[i] =  fLyC(z[i], False) 

        deps_dlogfLy1 = (1 -  np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) )
        deps_dlogfLy2 = np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 

        J = np.array((deps_dlogfLy1,deps_dlogfLy2))

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))
        sigma_eps_b_opt[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b_opt,J]))


    plt.plot(z, fid_fLy,color='k')
    #plt.fill_between(z, fid_fLy - sigma_eps_b_opt, fid_fLy + sigma_eps_b_opt, color=color_ULTRASAT, alpha = 0.9, label=r'$\rm Optimistic\,bias\,prior$')
    plt.fill_between(z, fid_fLy - sigma_eps_b, fid_fLy + sigma_eps_b, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm Conservative\,bias\,prior$')
    #plt.plot(z, fid_fLy - sigma_eps_G,color_FUV,alpha=0.6,linestyle=(0,(5,1)),label=r'$\rm C18\,)(reproduced)$')
    #plt.plot(z, fid_fLy + sigma_eps_G,color_FUV,alpha=0.6,linestyle=(0,(5,1)))


    plt.ylabel(r'$f_{\rm LyC}$',fontsize=fontsize_label)
    plt.legend(loc=7,fontsize=fontsize_legend)
    plt.xlabel(r'$z$',fontsize=fontsize_label)
    plt.xlim(z[0],z[-1])
    plt.ylim(0.1,1)
    plt.tick_params(axis='y', labelsize=fontsize_tick) 
    plt.tick_params(axis='x', labelsize=fontsize_tick) 

    print(max(sigma_eps_b/fid_fLy))
    print(max(sigma_eps_b_opt/fid_fLy))

    filename = 'results/PLOTS/EBL/errfLyC.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return




#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################

def plot_response():

    wavelenght, QE_dat = np.genfromtxt('dat/ULTRASAT.dat').T
    
    plt.figure()
    plt.plot(wavelenght,QE_dat,label=r'$ULTRASAT$',color=color_ULTRASAT)
    plt.legend(loc=1)
    plt.xlabel(r'$\lambda\,[{\rm A}]$',fontsize=fontsize)
    plt.ylabel(r'$\rm QE$',fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/QE_ultrasat.png',bbox_inches='tight')

    wavelenght_F = np.linspace(wavelenght_min('GALEX_FUV'),wavelenght_max('GALEX_FUV'))
    wavelenght_N = np.linspace(wavelenght_min('GALEX_NUV'),wavelenght_max('GALEX_NUV'))
    wavelenght = np.linspace(wavelenght_min('ULTRASAT'),wavelenght_max('ULTRASAT'))

    R = np.zeros(len(wavelenght))
    RN = np.zeros(len(wavelenght))
    RF = np.zeros(len(wavelenght))

    for i in range(len(wavelenght)):
        R[i] = Response(wavelenght[i],'ULTRASAT')
        RN[i] = Response(wavelenght_N[i],'GALEX_NUV')
        RF[i] = Response(wavelenght_F[i],'GALEX_FUV')

    plt.figure()
    plt.plot(wavelenght,R-min(R),label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    plt.plot(wavelenght_N,RN-min(RN),label=r'$\rm GALEX\, NUV$',color=color_NUV)
    plt.plot(wavelenght_F,RF-min(RF),label=r'$\rm GALEX\, FUV$',color=color_FUV)

    plt.yticks([])

    plt.legend(loc=1,ncol=1)
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize)
    plt.ylabel(r'$R(\lambda_{\rm obs})$',fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/response.png',bbox_inches='tight')

    plt.show()

    return 


def plot_z_signal(vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False):

    # reproduce fig 2, mid panel and fig 8

    z = [0,0.3,1.,2]#,5,10]

    rest_wave = np.linspace(700,3500,400)
    for zv in range(len(z)):
        s = np.zeros(len(rest_wave))
        wave = np.zeros((len(rest_wave)))
        for i in range(len(rest_wave)):
            wave[i] = rest_wave[i] * (1.+z[zv])
            if z[zv] <= 1:
                detector = 'GALEX_NUV' # same for galex_fuv in EW
            else:
                detector = 'ULTRASAT' # same for hetdex and spherex
            s[i] = signal(rest_wave[i]*u.AA,z[zv],detector,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,True).value 

        plt.plot(wave, s, color= colors[zv], label=r'$z = %g$'%z[zv])

    plt.yscale('log')
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm \AA}]$',fontsize=fontsize)
    plt.ylabel(r'$\epsilon(\nu,z)\,[{\rm erg\,s^{-1}Hz^{-1}Mpc^{-3}}]$',fontsize=fontsize)
    plt.legend(loc=4,ncol=2)
    plt.xlim(700,8000)
    plt.ylim(3e24,2e28)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/emissivity.png',bbox_inches='tight')
    plt.show()

    return 



def plot_signal(vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False):

    dJ_U =  dJdz(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_ULTRASAT_reduced.dat')

    dJ_N =  dJdz(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_GALEX_NUV_reduced.dat')

    dJ_F =  dJdz(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,filename='results/EBL/dJdz_GALEX_FUV_reduced.dat')

    bJ_U = bJ_z(z_gals('DESI'),detector='ULTRASAT',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = val_bias,filename='results/EBL/bJ_ULTRASAT_reduced.dat')

    bJ_N = bJ_z(z_gals('SDSS'),detector='GALEX_NUV',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = val_bias,filename='results/EBL/bJ_GALEX_NUV_reduced.dat')

    bJ_F = bJ_z(z_gals('SDSS'),detector='GALEX_FUV',run=False,vals_eps1500=vals_eps1500,vals_alpha1500=vals_alpha1500,vals_alpha1100=vals_alpha1100,val_EW=val_EW,val_flyc=val_flyc,val_alpha900=val_alpha900,val_bias = val_bias,filename='results/EBL/bJ_GALEX_FUV_reduced.dat')


    window_size = 70  # Adjust this value based on your preference
    bJU_smoothed = moving_average(dJ_U*bJ_U, window_size)
    bJGN_smoothed = moving_average(dJ_N*bJ_N, window_size)
    bJGF_smoothed = moving_average(dJ_F*bJ_F, window_size)

    plt.figure(figsize=(20,12))
    plt.plot(z_gals('SDSS')[:len(bJGF_smoothed)],bJGF_smoothed,label=r'$\rm GALEX\, FUV$',color=color_FUV)
    plt.plot(z_gals('SDSS')[:len(bJGN_smoothed)],bJGN_smoothed,label=r'$\rm GALEX\, NUV$',color=color_NUV)
    plt.plot(z_gals('DESI')[:len(bJU_smoothed)],bJU_smoothed,label=r'$\rm ULTRASAT$',color=color_ULTRASAT)
    #plt.ylim(-10,200)
    plt.xlabel(r'$z$',fontsize=fontsize*1.2)
    plt.ylabel(r'$b_JdJ_{\nu_{\rm obs}}/dz\,[{\rm Jy/sr}]$',fontsize=fontsize*1.2)
    plt.legend(bbox_to_anchor = (.95,-0.1), ncol=3,fontsize=fontsize*1.2)
    plt.xticks(fontsize=1.2*.8*fontsize)
    plt.yticks(fontsize=1.2*.8*fontsize)
    
    plt.xlim(0,2.3)

    plt.tight_layout()
    plt.savefig('results/PLOTS/EBL/bJdJdz.png',bbox_inches='tight')


    plt.show()
    return 



def plot_err_ioncont_wb(z, lambda_val, use_pars_fid = pars_all,group_vox=True,run=False, prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)
    nu_912 = nu_from_lambda(912*u.AA)

    if use_pars_fid == pars_original_c18_fid:
        required_pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100']

    else:
        required_pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','alpha_900','log_fLyC_1','log_fLyC_2']

    fid_eps_1500 = fiducials['eps1500'][0]
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]
    fid_alpha900 = fiducials['alpha900']

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    if prior == 'optimistic':
        sigma_b1500 = 0.01
        sigma_gammabnu = 0.3
        sigma_gammabz = 0.1
    elif prior == 'conservative':
        sigma_b1500 = 0.05
        sigma_gammabnu = 1.3
        sigma_gammabz = 0.3

    if prior:
        for j in range(len(use_pars_fid)):
            if galex_detector == 'SDSS':
                if use_pars_fid[j] == 'C_alpha_1500':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
                if use_pars_fid[j] == 'C_alpha_1100':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'gamma_bv':
                F_N[j,j] += 1/sigma_gammabnu**2
                F_F[j,j] += 1/sigma_gammabnu**2
                F_U[j,j] += 1/sigma_gammabnu**2
            if use_pars_fid[j] == 'gamma_bz':
                F_N[j,j] += 1/sigma_gammabz**2
                F_F[j,j] += 1/sigma_gammabz**2
                F_U[j,j] += 1/sigma_gammabz**2
            if use_pars_fid[j] == 'bias_1500_0':
                F_N[j,j] += 1/sigma_b1500**2
                F_F[j,j] += 1/sigma_b1500**2
                F_U[j,j] += 1/sigma_b1500**2

    F_galex = F_N + F_F
    F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)

    all_inv_F_G = np.linalg.inv(F_galex)
    all_inv_F_U = np.linalg.inv(F_U)
    all_inv_F_b = np.linalg.inv(F_both)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    inv_F_U = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
       id_i = use_pars_fid.index(required_pars[i])
       for j in range(len(required_pars)):
           id_j = use_pars_fid.index(required_pars[j])
           inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
           inv_F_U[i,j] = all_inv_F_U[id_i,id_j]
           inv_F_b[i,j] = all_inv_F_b[id_i,id_j]
#

    fid_eps = np.zeros(len(z))

    sigma_eps_G = np.zeros(len(z))
    sigma_eps_U = np.zeros(len(z))
    sigma_eps_b = np.zeros(len(z))

    for i in range(len(z)):

        A =  fLyC(z[i], False) 
        B =  (nu_912 / nu_1216)**(fid_alpha1100_0 + fid_C1100*np.log10(1+z[i])) * (use_nu / nu_912)**fid_alpha900

        C = fid_eps_1500 * (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp = C * ( A * B )
        fid_eps[i] = temp 

        deps_deps1500 = temp / fid_eps_1500

        deps_dgamma = temp * np.log(1+z[i])
        deps_dalpha1500 = temp * np.log(nu_1216 / nu_1500)
        deps_dC1500 = temp * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100 = temp * np.log(nu_912 / nu_1216)
        deps_dC1100 = temp * np.log(nu_912 / nu_1216) * np.log10(1+z[i])
        if not use_pars_fid == pars_original_c18_fid:
            deps_dlogfLy1 = C * B * (1 -  np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) )
            deps_dlogfLy2 = C * B * np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 
            deps_dalpha900 = temp * np.log(use_nu / nu_912)

            J = np.array((deps_deps1500,deps_dgamma,deps_dalpha1500,deps_dC1500,deps_dalpha1100,deps_dC1100,deps_dalpha900,deps_dlogfLy1,deps_dlogfLy2))
        else:
            J = np.array((deps_deps1500,deps_dgamma,deps_dalpha1500,deps_dC1500,deps_dalpha1100,deps_dC1100))

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))

    return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps




def plot_err_line_surrounding_wb(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False,prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)


    required_pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100']

    fid_eps_1500 = fiducials['eps1500'][0]
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)


    if prior == 'optimistic':
        sigma_b1500 = 0.01
        sigma_gammabnu = 0.3
        sigma_gammabz = 0.1
    elif prior == 'conservative':
        sigma_b1500 = 0.05
        sigma_gammabnu = 1.3
        sigma_gammabz = 0.3

    if prior:
        for j in range(len(use_pars_fid)):
            if galex_detector == 'SDSS':
                if use_pars_fid[j] == 'C_alpha_1500':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
                if use_pars_fid[j] == 'C_alpha_1100':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'gamma_bv':
                F_N[j,j] += 1/sigma_gammabnu**2
                F_F[j,j] += 1/sigma_gammabnu**2
                F_U[j,j] += 1/sigma_gammabnu**2
            if use_pars_fid[j] == 'gamma_bz':
                F_N[j,j] += 1/sigma_gammabz**2
                F_F[j,j] += 1/sigma_gammabz**2
                F_U[j,j] += 1/sigma_gammabz**2
            if use_pars_fid[j] == 'bias_1500_0':
                F_N[j,j] += 1/sigma_b1500**2
                F_F[j,j] += 1/sigma_b1500**2
                F_U[j,j] += 1/sigma_b1500**2

    F_galex = F_N + F_F
    F_both = np.zeros((len(F_N),len(F_N)))

    # since EW is parametrized independently using (0.3,1) and (1,2) , we decide to rely on the (1,2) choice and ignore the EW_0.3 parameter
    for i in range(len(F_N)):
        if use_pars_fid[i][:2] != 'EW':
            for j in range(len(F_N)):
                if use_pars_fid[j][:2] != 'EW':
                    F_both[i,j] = F_galex[i,j] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z1':
                    F_both[i,j] = F_galex[i,use_pars_fid.index('EW_z2')] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z2':
                    F_both[i,j] = F_U[i,j]
        elif use_pars_fid[i] == 'EW_z1':
            for j in range(len(F_N)):
                if use_pars_fid[j][:2] != 'EW':
                    F_both[i,j] = F_galex[use_pars_fid.index('EW_z2'),j] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z1':
                    F_both[i,j] = F_galex[use_pars_fid.index('EW_z2'),use_pars_fid.index('EW_z2')] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z2':
                    F_both[i,j] = F_U[i,j]
        elif use_pars_fid[i] == 'EW_z2':
            for j in range(len(F_N)):
                F_both[i,j] = F_U[i,j]

    #F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)

    all_inv_F_G = np.linalg.inv(F_galex)
    all_inv_F_U = np.linalg.inv(F_U)
    all_inv_F_b = np.linalg.inv(F_both)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    inv_F_U = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
       id_i = use_pars_fid.index(required_pars[i])
       for j in range(len(required_pars)):
           id_j = use_pars_fid.index(required_pars[j])
           inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
           inv_F_U[i,j] = all_inv_F_U[id_i,id_j]
           inv_F_b[i,j] = all_inv_F_b[id_i,id_j]
#
    fid_eps_galex = np.zeros(len(z))
    fid_eps_ultrasat = np.zeros(len(z))

    sigma_eps_G = np.zeros(len(z))
    sigma_eps_U = np.zeros(len(z))
    sigma_eps_b = np.zeros(len(z))

    for i in range(len(z)):

        B =  (use_nu / nu_1216)**(fid_alpha1100_0 + fid_C1100*np.log10(1+z[i]))

        C = fid_eps_1500 * (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp_galex = C * ( B )

        temp_ultrasat = C * ( B )

        fid_eps_galex[i] = temp_galex
        fid_eps_ultrasat[i] = temp_ultrasat

        deps_deps_galex = temp_galex / fid_eps_1500
        deps_dgamma_galex = temp_galex * np.log(1+z[i])
        deps_dalpha1500_galex = temp_galex * np.log(nu_1216 / nu_1500)
        deps_dC1500_galex = temp_galex * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_galex = C  * np.log(use_nu / nu_1216)
        deps_dC1100_galex = C  * np.log(use_nu / nu_1216) * np.log10(1+z[i])

        deps_deps_ultrasat = temp_ultrasat / fid_eps_1500
        deps_dgamma_ultrasat = temp_ultrasat * np.log(1+z[i])
        deps_dalpha1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500)
        deps_dC1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_ultrasat = C * np.log(use_nu / nu_1216)
        deps_dC1100_ultrasat = C  * np.log(use_nu / nu_1216) * np.log10(1+z[i])


        J_galex = np.array((deps_deps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,))

        J_ultrasat = np.array((deps_deps_ultrasat,deps_dgamma_ultrasat,deps_dalpha1500_ultrasat,deps_dC1500_ultrasat,deps_dalpha1100_ultrasat,deps_dC1100_ultrasat))

        J_both = np.array((deps_deps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex)) 

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J_galex,inv_F_G,J_galex]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J_ultrasat,inv_F_U,J_ultrasat]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J_both,inv_F_b,J_both]))

    return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps_ultrasat




def plot_err_line_wb(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False,prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)


    required_pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2']

    fid_eps_1500 = fiducials['eps1500'][0]
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)


    if prior == 'optimistic':
        sigma_b1500 = 0.01
        sigma_gammabnu = 0.3
        sigma_gammabz = 0.1
    elif prior == 'conservative':
        sigma_b1500 = 0.05
        sigma_gammabnu = 1.3
        sigma_gammabz = 0.3

    if prior:
        for j in range(len(use_pars_fid)):
            if galex_detector == 'SDSS':
                if use_pars_fid[j] == 'C_alpha_1500':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
                if use_pars_fid[j] == 'C_alpha_1100':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'gamma_bv':
                F_N[j,j] += 1/sigma_gammabnu**2
                F_F[j,j] += 1/sigma_gammabnu**2
                F_U[j,j] += 1/sigma_gammabnu**2
            if use_pars_fid[j] == 'gamma_bz':
                F_N[j,j] += 1/sigma_gammabz**2
                F_F[j,j] += 1/sigma_gammabz**2
                F_U[j,j] += 1/sigma_gammabz**2
            if use_pars_fid[j] == 'bias_1500_0':
                F_N[j,j] += 1/sigma_b1500**2
                F_F[j,j] += 1/sigma_b1500**2
                F_U[j,j] += 1/sigma_b1500**2


    F_galex = F_N + F_F
    F_both = np.zeros((len(F_N),len(F_N)))


    # since EW is parametrized independently using (0.3,1) and (1,2) , we decide to rely on the (1,2) choice and ignore the EW_0.3 parameter
    for i in range(len(F_N)):
        if use_pars_fid[i][:2] != 'EW':
            for j in range(len(F_N)):
                if use_pars_fid[j][:2] != 'EW':
                    F_both[i,j] = F_galex[i,j] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z1':
                    F_both[i,j] = F_galex[i,use_pars_fid.index('EW_z2')] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z2':
                    F_both[i,j] = F_U[i,j]
        elif use_pars_fid[i] == 'EW_z1':
            for j in range(len(F_N)):
                if use_pars_fid[j][:2] != 'EW':
                    F_both[i,j] = F_galex[use_pars_fid.index('EW_z2'),j] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z1':
                    F_both[i,j] = F_galex[use_pars_fid.index('EW_z2'),use_pars_fid.index('EW_z2')] + F_U[i,j]
                elif use_pars_fid[j] == 'EW_z2':
                    F_both[i,j] = F_U[i,j]
        elif use_pars_fid[i] == 'EW_z2':
            for j in range(len(F_N)):
                F_both[i,j] = F_U[i,j]

    #F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)

    all_inv_F_G = np.linalg.inv(F_galex)
    all_inv_F_U = np.linalg.inv(F_U)
    all_inv_F_b = np.linalg.inv(F_both)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    inv_F_U = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
       id_i = use_pars_fid.index(required_pars[i])
       for j in range(len(required_pars)):
           id_j = use_pars_fid.index(required_pars[j])
           inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
           inv_F_U[i,j] = all_inv_F_U[id_i,id_j]
           inv_F_b[i,j] = all_inv_F_b[id_i,id_j]
#
    fid_eps_galex = np.zeros(len(z))
    fid_eps_ultrasat = np.zeros(len(z))

    sigma_eps_G = np.zeros(len(z))
    sigma_eps_U = np.zeros(len(z))
    sigma_eps_b = np.zeros(len(z))

    for i in range(len(z)):

        A_galex =  EW_val(z[i],'GALEX_NUV', False) / (0.005*nu_1216) * (use_nu**2 / cu.c.to(u.AA/u.s))
        B =  (use_nu / nu_1216)**(fid_alpha1100_0 + fid_C1100*np.log10(1+z[i]))

        C = fid_eps_1500* (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp_galex = C * ( A_galex + B )

        A_ultrasat =  EW_val(z[i],'ULTRASAT', False) / (0.005*nu_1216) * (use_nu**2 / cu.c.to(u.AA/u.s))
        temp_ultrasat = C * ( A_ultrasat + B )

        fid_eps_galex[i] = temp_galex
        fid_eps_ultrasat[i] = temp_ultrasat

        deps_deps_galex = temp_galex / fid_eps_1500
        deps_dgamma_galex = temp_galex * np.log(1+z[i])
        deps_dalpha1500_galex = temp_galex * np.log(nu_1216 / nu_1500)
        deps_dC1500_galex = temp_galex * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_galex = C * A_galex * np.log(use_nu / nu_1216)
        deps_dC1100_galex = C * A_galex * np.log(use_nu / nu_1216) * np.log10(1+z[i])
        deps_dEW03_galex = C * B * (1 -  np.log10((1+z[i])/(1+0.3)) / np.log10((1+1)/(1+0.3)) )
        deps_dEW1_galex = C * B * np.log10((1+z[i])/(1+0.3)) / np.log10((1+1)/(1+0.3)) 

        deps_deps_ultrasat = temp_ultrasat / fid_eps_1500
        deps_dgamma_ultrasat = temp_ultrasat * np.log(1+z[i])
        deps_dalpha1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500)
        deps_dC1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_ultrasat = C * A_ultrasat * np.log(use_nu / nu_1216)
        deps_dC1100_ultrasat = C * A_ultrasat * np.log(use_nu / nu_1216) * np.log10(1+z[i])
        deps_dEW1_ultrasat = C * B * (1 - np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)))
        deps_dEW2_ultrasat = C * B * np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 


        J_galex = np.array((deps_deps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,deps_dEW03_galex,deps_dEW1_galex))

        J_ultrasat = np.array((deps_deps_ultrasat,deps_dgamma_ultrasat,deps_dalpha1500_ultrasat,deps_dC1500_ultrasat,deps_dalpha1100_ultrasat,deps_dC1100_ultrasat,deps_dEW1_ultrasat,deps_dEW2_ultrasat))

        J_both = np.array((deps_deps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,deps_dEW1_ultrasat,deps_dEW2_ultrasat)) 

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J_galex,inv_F_G,J_galex]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J_ultrasat,inv_F_U,J_ultrasat]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J_both,inv_F_b,J_both]))


    return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps_ultrasat



def plot_err_noninonizing_cont_wb(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False, prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)


    # required vals: log10(b1500eps1500), gamma1500, alpha1500, C1500

    required_pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500']

    fid_log10_beps = fiducials['eps1500'][0]
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)
    

    if prior == 'optimistic':
        sigma_b1500 = 0.01
        sigma_gammabnu = 0.3
        sigma_gammabz = 0.1
    elif prior == 'conservative':
        sigma_b1500 = 0.05
        sigma_gammabnu = 1.3
        sigma_gammabz = 0.3

    if prior:
        for j in range(len(use_pars_fid)):
            if galex_detector == 'SDSS':
                if use_pars_fid[j] == 'C_alpha_1500':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
                if use_pars_fid[j] == 'C_alpha_1100':
                    F_N[j,j] += 1/1.5**2
                    F_F[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'gamma_bv':
                F_N[j,j] += 1/sigma_gammabnu**2
                F_F[j,j] += 1/sigma_gammabnu**2
                F_U[j,j] += 1/sigma_gammabnu**2
            if use_pars_fid[j] == 'gamma_bz':
                F_N[j,j] += 1/sigma_gammabz**2
                F_F[j,j] += 1/sigma_gammabz**2
                F_U[j,j] += 1/sigma_gammabz**2
            if use_pars_fid[j] == 'bias_1500_0':
                F_N[j,j] += 1/sigma_b1500**2
                F_F[j,j] += 1/sigma_b1500**2
                F_U[j,j] += 1/sigma_b1500**2

    F_galex = F_N + F_F
    F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)

    all_inv_F_G = np.linalg.inv(F_galex)
    all_inv_F_U = np.linalg.inv(F_U)
    all_inv_F_b = np.linalg.inv(F_both)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    inv_F_U = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
        id_i = use_pars_fid.index(required_pars[i])
        for j in range(len(required_pars)):
            id_j = use_pars_fid.index(required_pars[j])
            inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
            inv_F_U[i,j] = all_inv_F_U[id_i,id_j]
            inv_F_b[i,j] = all_inv_F_b[id_i,id_j]


    fid_eps_1500 = np.zeros(len(z))
    sigma_eps1500_G = np.zeros(len(z))
    sigma_eps1500_U = np.zeros(len(z))
    sigma_eps1500_b = np.zeros(len(z))

    for i in range(len(z)):

        temp = (fid_log10_beps) * (1+z[i])**fid_gamma1500 * (use_nu / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        fid_eps_1500[i] = temp

        deps_deps = temp / fid_log10_beps
        deps_dgamma = temp * np.log(1+z[i])
        deps_dalpha = temp * np.log(use_nu / nu_1500)
        deps_dC = temp * np.log(use_nu / nu_1500) * np.log10(1+z[i])

        J = np.array((deps_deps,deps_dgamma,deps_dalpha,deps_dC))

        sigma_eps1500_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_eps1500_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_eps1500_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))


    return [sigma_eps1500_G, sigma_eps1500_U, sigma_eps1500_b], fid_eps_1500


