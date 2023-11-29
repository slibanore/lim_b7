from EBL_observable import * 


def sigma_N(detector, group_vox=True):

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV':
        if not group_vox:
            px_size = 5*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 50*u.arcsec

        mag_limit = 20.5

    elif detector == 'ULTRASAT':
        if not group_vox:
            px_size = 5.4*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 5.4*u.arcsec
    
        mag_limit = 23.5

    else: 
        print('Detector not recognized!')
        return -1 

    fnu = 10**((mag_limit - 8.90)/(-2.5))*u.Jy
    fnu_px = fnu/(((px_size)**2).to(u.steradian))

    noise_per_voxel = fnu_px/5

    return noise_per_voxel



def sigma_wz(z, detector, gal_survey, group_vox):

    delta_zc = 1e-2

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV':
        if not group_vox:
            px_size = 5*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 50.*u.arcsec
        #sky_LIM = (5500*u.deg**2).to(u.steradian)

    elif detector == 'ULTRASAT':
        if not group_vox:
            px_size = 5.4*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 5.4*u.arcsec

        #sky_LIM = (40000*u.deg**2).to(u.steradian)

    #Nvox = sky_LIM / (px_size**2).to(u.steradian)
    #Nvox_std = 1 / (px_size**2).to(u.steradian) #Nvox / sky_LIM
    #Nvox_inthetamax =  np.pi* use_thetamax**2 /  (px_size**2).to(u.steradian)  #Nvox_std * np.pi* use_thetamax**2

    Avox = (((px_size)**2).to(u.steradian))
#    Ngal_ref = dNgdz(z,gal_survey) * delta_zi(gal_survey)

#    noise = delta_zi(gal_survey) / delta_zc * use_thetamax * np.sqrt(np.pi) * np.sqrt(Jnu_monopole(detector)**2 + sigma_N(detector, group_vox).value**2) / np.sqrt((Ngal_ref*Nvox_inthetamax))
# 
    dzi = delta_zi(gal_survey)
    #if gal_survey == 'DESI' and z >= 1.6:
    #    dzi *= 7

    noise = np.sqrt(dzi) / delta_zc * np.sqrt(Avox * (Jnu_monopole(detector)**2 + sigma_N(detector, group_vox).value**2)) / np.sqrt(dNgdz(z,gal_survey)) 

    return noise



def plot_noise_grouped():

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
    plt.xlabel(r'$L_{\rm pix}\,[{\rm arcsec/pix}]$',fontsize=fontsize)
    plt.ylabel(r'$\mathcal{N}_{\rm CBR}$',fontsize=fontsize)
    plt.legend(loc=1,ncol=1)

    plt.tight_layout()

    if scale_physical_max == 300*u.Mpc: 
        filename = 'results/PLOTS/EBL/noise_voxgroup_thetamax.png'
    else:
        filename = 'results/PLOTS/EBL/noise_voxgroup.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 


def plot_signal_and_noise(detector,gal_survey,reduced_z = False):

    if reduced_z:
        use_z = z_gals_interp if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1
        reduced_label = '_reduced'
    else:
        use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1
        reduced_label = ''

    signal = np.zeros(len(use_z))

    filename = 'results/EBL/wJg_' + detector +',' + gal_survey + reduced_label + '.dat'

    for i in tqdm(range(len(use_z))):
        signal[i] = abs(wJgz(use_z[i],detector,gal_survey,False,filename=filename))

    nsingle = np.zeros(len(use_z))
    ngrouped = np.zeros(len(use_z))
    if detector == 'ULTRASAT':
        ngrouped50 = np.zeros(len(use_z))
    for i in tqdm(range(len(use_z))):
        nsingle[i] = sigma_wz(use_z[i],detector,gal_survey, False).value
        ngrouped[i] = sigma_wz(use_z[i],detector,gal_survey, True).value
        if detector == 'ULTRASAT':
            ngrouped50[i] = sigma_wz(use_z[i],detector,gal_survey, 50).value

    color = color_ULTRASAT if detector == 'ULTRASAT' else color_NUV  if detector == 'GALEX_NUV' else color_FUV if detector == 'GALEX_FUV' else -1

    plt.plot(use_z,signal,label=r'$\rm %s\times$'%detector + r'$\rm %s$'%gal_survey, color=color)

    #label_grouped = r'$3\times 3 \,{\rm voxels}$' if detector == 'ULTRASAT' else r'$\rm Voxels = 50\,arcsec$'  

    #plt.plot(use_z,nsingle,'k:', label=r'$\rm Noise\,per\,voxel$')
    
    #label_grouped = r'${\rm Noise:}\, 3\times 3 \rm voxels$' if detector == 'ULTRASAT' else r'$\rm Voxels = 50\,arcsec$'  
    plt.plot(use_z,ngrouped,'k-',label=r'$\rm Noise$')#label_grouped)
    #if detector == 'ULTRASAT':
    #    plt.plot(use_z,ngrouped50,'k--',label= r'$\rm Noise:\, 50\,arcsec\, voxels$')

    plt.xlabel(r'$z$')
    plt.ylabel(r'$|\sigma_{w_{\tilde{J}_\nu{\rm g}}(z)}|$')
    plt.yscale('log')
    plt.legend(loc=1)

    plt.tight_layout()

    filefig = 'results/PLOTS/EBL/noise_' + detector + ',' + gal_survey + '.png'
    
    #plt.savefig(filefig,bbox_inches='tight')
    plt.show()

    return 

def plot_noise(reduced_z = False):

    if reduced_z:
        use_z = z_gals_interp
        #use_z_ultrasat = z_gals_interp
        use_z_ultrasat_desi = z_gals_interp
        reduced_label = '_reduced'
    else:
        use_z = z_gals('SDSS')
        reduced_label = '_reduced'
        #use_z_ultrasat = z_gals('SPHEREx')
        use_z_ultrasat_desi = z_gals('DESI')

    sn = np.zeros(len(use_z))
    sf = np.zeros(len(use_z))
    snD = np.zeros(len(use_z))
    sfD = np.zeros(len(use_z))
    #sU = np.zeros(len(use_z_ultrasat))
    sUD = np.zeros(len(use_z_ultrasat_desi))
    nG = np.zeros(len(use_z))
    nGD = np.zeros(len(use_z))
    #nU = np.zeros(len(use_z_ultrasat))
    nUD = np.zeros(len(use_z_ultrasat_desi))

    #filename_ULT = 'results/EBL/wJg_ULTRASAT,SPHEREx' + reduced_label + '.dat'
    plt.figure()    
    filename_ULTDESI = 'results/EBL/wJg_ULTRASAT,DESI' + reduced_label + '.dat'
        
    for i in tqdm(range(len(use_z))):
        sn[i] = abs(wJgz(use_z[i],'GALEX_NUV','SDSS',False,filename='results/EBL/wJg_GALEX_NUV,SDSS' + reduced_label + '.dat'))
        sf[i] = abs(wJgz(use_z[i],'GALEX_FUV','SDSS',False,filename='results/EBL/wJg_GALEX_FUV,SDSS'+ reduced_label+'.dat'))
        nG[i] = sigma_wz(use_z[i],'GALEX_NUV','SDSS',True).value
        snD[i] = abs(wJgz(use_z[i],'GALEX_NUV','DESI',False,filename='results/EBL/wJg_GALEX_NUV,DESI' + reduced_label + '.dat'))
        sfD[i] = abs(wJgz(use_z[i],'GALEX_FUV','DESI',False,filename='results/EBL/wJg_GALEX_FUV,DESI'+ reduced_label+'.dat'))
        nGD[i] = sigma_wz(use_z[i],'GALEX_NUV','DESI',True).value
    #for i in tqdm(range(len(use_z_ultrasat))):
    #    sU[i] = abs(wJgz(use_z_ultrasat[i],'ULTRASAT','SPHEREx',False,filename=filename_ULT))
    #    nU[i] = sigma_wz(use_z[i],'ULTRASAT','SPHEREx',True).value
    for i in tqdm(range(len(use_z_ultrasat_desi))):
        sUD[i] = abs(wJgz(use_z_ultrasat_desi[i],'ULTRASAT','DESI',False,filename=filename_ULTDESI))
        nUD[i] = sigma_wz(use_z_ultrasat_desi[i],'ULTRASAT','DESI',True).value

    window_size = 3  # Adjust this value based on your preference
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
    #plt.plot(use_z_ultrasat,sU/nU,label=r'$ULTRASAT\times SPHEREx$',color=color_ULTRASAT)
    plt.plot(use_z_ultrasat_desi[:len(U_smoothed)],U_smoothed,label=r'$\rm ULTRASAT\times DESI$',color=color_ULTRASAT)
    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'${\bar{\omega}_{\tilde{J}{\rm g}}(z)}/\mathcal{N}_{\rm CBR}$',fontsize=fontsize)
    plt.legend(loc=1,ncol=1,fontsize=fontsize*.8)

    plt.ylim(1,5e5)
#    plt.hlines(1,use_z[0],use_z[-1],linewidth=1.,color='k')
    plt.xlim(use_z[0],min(use_z[-1],use_z_ultrasat_desi[:len(U_smoothed)][-1]))
    plt.tight_layout()


    filename = 'results/PLOTS/EBL/wJg_noise' + reduced_label + '.png'
    
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 



pars_fid = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2','log_fLyC_1','log_fLyC_2','bias_1500_0','gamma_bv','gamma_bz','alpha_900']
#
fid_vals_det = lambda detector: [fiducials['eps1500'][0],
      fiducials['eps1500'][1],
      fiducials['alpha1500'][0],
      fiducials['alpha1500'][1],
      fiducials['alpha1100'][0],
      fiducials['alpha1100'][1],
      fiducials['EW'](detector)[0],
      fiducials['EW'](detector)[1],
      fiducials['fescape'][0],
      fiducials['fescape'][1],
      fiducials['bias'][0],
      fiducials['bias'][1],
      fiducials['bias'][2],
        fiducials['alpha900']
      ]


delta_der = 0.01

def ders(parameter,detector,gal_survey,reduced_z = False):

    reduced_label = '_reduced' if reduced_z else ''

    print('Doing der: ' + str(parameter))

    filename = 'results/EBL/der/d' + parameter + '_' + detector + ',' + gal_survey + reduced_label + '.dat'

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT':
        filename_wm = 'results/EBL/wmz_' + detector + reduced_label + '.dat'
    else:
        filename_wm = -1

    filename_dJ = 'results/EBL/der/dJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'
    filename_bJ = 'results/EBL/der/bJ_d' + parameter + '_' + detector + ',' + gal_survey  + reduced_label +  '.dat'

    filename_dJ_fid = 'results/EBL/dJdz_' + detector  + reduced_label +  '.dat'
    filename_bJ_fid = 'results/EBL/bJ_' + detector  + reduced_label +  '.dat' 

    fid_vals = fid_vals_det(detector)
    if os.path.exists(filename):
        use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

        z_arr, der_wz_arr = np.genfromtxt(filename)
        der_wz = interp1d(z_arr,der_wz_arr)(use_z)

    else:

        if reduced_z:
            use_z = z_gals_interp if detector == 'GALEX_NUV' or detector ==    'GALEX_FUV' or detector == 'ULTRASAT' else -1
        else:
            use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector ==    'GALEX_FUV' or detector == 'ULTRASAT' else -1

        use_pars_up = []
        for i in range(len(pars_fid)):
            
            if pars_fid[i] == parameter:
                par_up = fid_vals[i]*(1+delta_der)
                step = fid_vals[i]*delta_der
            else:
                par_up = fid_vals[i]

            use_pars_up.append(par_up)

        with Pool(6) as pool:

            print('Doing parameter up')
            dJ_up_f = partial(dJdz,detector=detector,run=True,\
                    vals_eps1500=[use_pars_up[pars_fid.index('eps_1500_0')],use_pars_up[pars_fid.index('gamma_1500')]],\
                    vals_alpha1500=[use_pars_up[pars_fid.index('alpha_1500_0')],use_pars_up[pars_fid.index('C_alpha_1500')]],\
                    vals_alpha1100=[use_pars_up[pars_fid.index('alpha_1100_0')],use_pars_up[pars_fid.index('C_alpha_1100')]],\
                    val_EW=[use_pars_up[pars_fid.index('EW_z1')],use_pars_up[pars_fid.index('EW_z2')]],\
                    val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],\
                    val_alpha900=use_pars_up[pars_fid.index('alpha_900')],\
                    filename=filename_dJ)
            
            dJ_up = pool.map(dJ_up_f, use_z)

            print('Doing bias up')
            bJ_up_f = partial(bJ_z,detector=detector,run=True,\
                    vals_eps1500=[use_pars_up[pars_fid.index('eps_1500_0')],use_pars_up[pars_fid.index('gamma_1500')]],\
                    vals_alpha1500=[use_pars_up[pars_fid.index('alpha_1500_0')],use_pars_up[pars_fid.index('C_alpha_1500')]],\
                    vals_alpha1100=[use_pars_up[pars_fid.index('alpha_1100_0')],use_pars_up[pars_fid.index('C_alpha_1100')]],\
                    val_EW=[use_pars_up[pars_fid.index('EW_z1')],use_pars_up[pars_fid.index('EW_z2')]],\
                    val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],\
                    val_alpha900=use_pars_up[pars_fid.index('alpha_900')],\
                    val_bias=[use_pars_up[pars_fid.index('bias_1500_0')],use_pars_up[pars_fid.index('gamma_bv')],use_pars_up[pars_fid.index('gamma_bz')]],\
                    filename=filename_bJ)
            
            bJ_up = pool.map(bJ_up_f, use_z)

        np.savetxt(filename_dJ,(use_z,dJ_up))
        np.savetxt(filename_bJ,(use_z,bJ_up))


        filename_fid = 'results/EBL/wJg_' + detector + ',' + gal_survey + reduced_label + '.dat'
       
        if os.path.exists(filename_fid):
                run_fid = False 
        else:
                run_fid = True       
       
        der_wz = np.zeros(len(use_z))
        for z in tqdm(range(len(use_z))):
            wz_up = wJgz(use_z[z],detector, gal_survey,
                True,
                vals_eps1500=[use_pars_up[pars_fid.index('eps_1500_0')],use_pars_up[pars_fid.index('gamma_1500')]],
                vals_alpha1500=[use_pars_up[pars_fid.index('alpha_1500_0')],use_pars_up[pars_fid.index('C_alpha_1500')]],
                vals_alpha1100=[use_pars_up[pars_fid.index('alpha_1100_0')],use_pars_up[pars_fid.index('C_alpha_1100')]],
                val_EW=[use_pars_up[pars_fid.index('EW_z1')],use_pars_up[pars_fid.index('EW_z2')]],
                val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],
                val_alpha900=use_pars_up[pars_fid.index('alpha_900')],
                val_bias=[use_pars_up[pars_fid.index('bias_1500_0')],use_pars_up[pars_fid.index('gamma_bv')],use_pars_up[pars_fid.index('gamma_bz')]],\
                #filename = filename, 
                filename_wm = filename_wm,
                filename_dJ = filename_dJ,
                filename_bJ =filename_bJ)

            wz = wJgz(use_z[z],detector,gal_survey,
                run=run_fid,
                vals_eps1500=False,
                vals_alpha1500=False,
                vals_alpha1100=False,
                val_EW=False,
                val_flyc=False,
                val_alpha900=False,
                val_bias=False,
                filename =filename_fid, 
                filename_wm = filename_wm,
                filename_dJ = filename_dJ_fid,
                filename_bJ =filename_bJ_fid)

            try:
                der_wz[z] = ((wz_up - wz)/step)
            except:
                der_wz[z] = ((wz_up - wz)/step).value

        np.savetxt(filename,(use_z, der_wz))
        if run_fid:
            np.savetxt(filename_fid,(use_z, (use_z,wz)))

    return der_wz


def run_ders(detector,gal_survey,reduced_z=False):

    for parameter in pars_fid:
        ders(parameter,detector,gal_survey,reduced_z)

    return


def plot_ders(pars, reduced_z = False):

    wn = np.zeros(len(z_gals('SDSS')))
    reduced_label = '_reduced' if reduced_z else ''

    use_filename_nuv_J = 'results/EBL/dJdz_GALEX_NUV' + reduced_label + '.dat'
    use_filename_nuv_b = 'results/EBL/bJ_GALEX_NUV' + reduced_label + '.dat'

    for i in (range(len(z_gals('SDSS')))):
        print('\nDoing z = ' + str(z_gals('SDSS')[i]))
        print('NUV')
        wn[i] = dJdz(z_gals('SDSS')[i], 'GALEX_NUV',run = False,\
                           vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = use_filename_nuv_J) * bJ_z(z_gals('SDSS')[i], 'GALEX_NUV', run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = use_filename_nuv_b)

    plt.figure()
    plt.plot(z_gals('SDSS'),(wn),color_NUV,label=r'$\rm NUV\times SDSS$')

    for p in pars:
        wder = np.zeros(len(z_gals('SDSS')))

        use_filename_wder_J = 'results/EBL/der/dJ_d' + p + '_' + 'GALEX_NUV,SDSS' + reduced_label +  '.dat'
        use_filename_wder_b = 'results/EBL/der/bJ_d' + p + '_' + 'GALEX_NUV,SDSS' + reduced_label +  '.dat'

        for i in (range(len(z_gals('SDSS')))):
            print('\nDoing z = ' + str(z_gals('SDSS')[i]))
            wder[i] = dJdz(z_gals('SDSS')[i], 'GALEX_NUV',run = False,\
                           vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = use_filename_wder_J) * bJ_z(z_gals('SDSS')[i], 'GALEX_NUV', run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = use_filename_wder_b)

        plt.plot(z_gals('SDSS'),wder,label=r'$\rm d%s$'%p)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\bar{\omega}_{J_\nu{\rm g}}(z)$',fontsize=fontsize)
    plt.legend(loc=1,ncol=2)

    plt.tight_layout()
    plt.show()

    return 


def plot_ders_ultr(pars, reduced_z = False):

    wn = np.zeros(len(z_gals('DESI')))
    reduced_label = '_reduced' if reduced_z else ''

    use_filename_nuv_J = 'results/EBL/dJdz_ULTRASAT' + reduced_label + '.dat'
    use_filename_nuv_b = 'results/EBL/bJ_ULTRASAT' + reduced_label + '.dat'

    for i in (range(len(z_gals('DESI')))):
        print('\nDoing z = ' + str(z_gals('DESI')[i]))
        print('NUV')
        wn[i] = dJdz(z_gals('DESI')[i], 'ULTRASAT',run = False,\
                           vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = use_filename_nuv_J) * bJ_z(z_gals('DESI')[i], 'ULTRASAT', run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = use_filename_nuv_b)

    plt.figure()
    plt.plot(z_gals('DESI'),(wn),color_NUV,label=r'$\rm ULTRASAT\times DESI$')

    for p in pars:
        wder = np.zeros(len(z_gals('DESI')))

        use_filename_wder_J = 'results/EBL/der/dJ_d' + p + '_' + 'ULTRASAT,DESI' + reduced_label +  '.dat'
        use_filename_wder_b = 'results/EBL/der/bJ_d' + p + '_' + 'ULTRASAT,DESI' + reduced_label +  '.dat'

        for i in (range(len(z_gals('DESI')))):
            print('\nDoing z = ' + str(z_gals('DESI')[i]))
            wder[i] = dJdz(z_gals('DESI')[i], 'ULTRASAT',run = False,\
                           vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = use_filename_wder_J) * bJ_z(z_gals('DESI')[i], 'ULTRASAT', run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = use_filename_wder_b)

        plt.plot(z_gals('DESI'),wder,label=r'$\rm d%s$'%p)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\bar{\omega}_{J_\nu{\rm g}}(z)$',fontsize=fontsize)
    plt.legend(loc=1,ncol=2)

    plt.tight_layout()
    plt.show()

    return 


def fisher_matrix(pars,detector,gal_survey,group_vox,run = False):

    if group_vox: 
        filename = 'results/EBL/FISHER_' + detector + ',' + gal_survey + '.dat'
    else:
        filename = 'results/EBL/FISHER_' + detector + ',' + gal_survey + '_singleVOX.dat'

    if not run and os.path.exists(filename):
        Fisher_start = np.genfromtxt(filename)

        Fisher = np.zeros((len(pars),len(pars)))
        for i in range(len(pars)):
            for j in range(len(pars)):
                Fisher[i,j] = Fisher_start[pars_fid.index(pars[i])][pars_fid.index(pars[j])]
        return Fisher

    Fisher = np.zeros((len(pars),len(pars)))
    
    use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

    der_pars = {}
    for p in range(len(pars)):
        #temp = ders(pars[p],detector,gal_survey,reduced_z = True)
        #all_zbin_der = interp1d(use_z,temp)
        all_zbin_der = ders(pars[p],detector,gal_survey,reduced_z = True)

        der_pars[pars[p]] = all_zbin_der 

    sigma2 = np.zeros(len(use_z))
    print('Doing sigma2')
    for zv in tqdm(range(len(use_z))):

        sigma2[zv] = sigma_wz(use_z[zv],detector,gal_survey,group_vox).value**2

    Fisher = np.zeros((len(pars),len(pars)))
    print('Doing Fisher')
    for p1 in tqdm(range(len(pars))):
        for p2 in range(len(pars)):
            if p1 <= p2:
                der_p1 = np.asarray(der_pars[pars[p1]])                
                der_p2 = np.asarray(der_pars[pars[p2]])        
                        
                temp = sum(der_p1*der_p2/sigma2)
                    
                Fisher[p1,p2] = temp
            else: 
                Fisher[p1,p2] = Fisher[p2,p1]

    if group_vox: 
        np.savetxt(filename,Fisher)
    else:
        np.savetxt(filename,Fisher)

    return Fisher


def Fisher_change_var(pars,detector,gal_survey,group_vox,run=False):

    F = fisher_matrix(pars_fid,detector,gal_survey,group_vox, run)

    deps_dP = np.log(10) * fid_vals_det(detector)[pars_fid.index('eps_1500_0')]
    db_dP = np.log(10) * fid_vals_det(detector)[pars_fid.index('bias_1500_0')]

    Fisher_prime = np.zeros((len(pars),len(pars)))
    
    for i in range(len(pars)):

        if pars[i] == 'log10_epsbias_1500_0':

            Fisher_prime[i,i] = F[pars_fid.index('eps_1500_0')][pars_fid.index('eps_1500_0')] * deps_dP**2 + 2*F[pars_fid.index('eps_1500_0')][pars_fid.index('bias_1500_0')] * deps_dP*db_dP + F[pars_fid.index('bias_1500_0')][pars_fid.index('bias_1500_0')] * db_dP**2

            for j in range(len(pars)):
                
                if i < j:
                    Fisher_prime[i,j] = F[pars_fid.index('eps_1500_0')][pars_fid.index(pars[j])] * deps_dP + F[pars_fid.index('bias_1500_0')][pars_fid.index(pars[j])] * db_dP

        else:
            for j in range(len(pars)):
                if i <= j:
                    Fisher_prime[i,j] = F[pars_fid.index(pars[i])][pars_fid.index(pars[j])]

    for i in range(len(pars)):
        for j in range(len(pars)):
            if i > j:
                Fisher_prime[i,j] = Fisher_prime[j,i]

    return Fisher_prime




def sigma_epsbias(group_vox=True, fix_bias = False):

    print('Doing contour plots')

    if fix_bias:
        pars = ['eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','log_fLyC_1','log_fLyC_2', 'gamma_bv','gamma_bz','EW_z1','EW_z2']

    else:
        pars = ['eps_1500_0','bias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','log_fLyC_1','log_fLyC_2', 'gamma_bv','gamma_bz','EW_z1','EW_z2']
#
    F_N = fisher_matrix(pars,'GALEX_NUV','SDSS',group_vox,False)
    F_F = fisher_matrix(pars,'GALEX_FUV','SDSS',group_vox,False)
    F_GALEX = F_N + F_F

    F_ULTRASAT = fisher_matrix(pars,'ULTRASAT','SPHEREx',group_vox,False)

    F_ULTRASAT_DESI = fisher_matrix(pars,'ULTRASAT','DESI',group_vox,False)

    F_both = np.zeros((len(pars)+1,len(pars)+1))
    F_both_DESI = np.zeros((len(pars)+1,len(pars)+1))
    for i in range(len(pars)):
        if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
            for j in range(len(pars)):
                if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                    F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                    F_both_DESI[i,j] = F_ULTRASAT_DESI[i,j] + F_GALEX[i,j]
                elif pars[j] == 'EW_z1':
                    F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                    F_both_DESI[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                elif pars[j] == 'EW_z2':
                    F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
                    F_both_DESI[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT_DESI[i,pars.index('EW_z1')]
            F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
            F_both_DESI[i,-1] = F_ULTRASAT_DESI[i,pars.index('EW_z2')]
        elif pars[i] == 'EW_z1':
            for j in range(len(pars)):
                F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
                F_both_DESI[i,j] = F_GALEX[pars.index('EW_z1'),j]
        elif pars[i] == 'EW_z2':
            for j in range(len(pars)):
                if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                    F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                    F_both_DESI[i,j] = F_ULTRASAT_DESI[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                elif pars[j] == 'EW_z1':
                    F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                    F_both_DESI[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                elif pars[j] == 'EW_z2':
                    F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
                    F_both_DESI[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT_DESI[pars.index('EW_z1'),pars.index('EW_z1')]
            F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
            F_both_DESI[i,-1] = F_ULTRASAT_DESI[pars.index('EW_z1'),pars.index('EW_z2')]
    for j in range(len(pars)):
        if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
            F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
            F_both_DESI[-1,j] = F_ULTRASAT_DESI[pars.index('EW_z2'),j]
        elif pars[j] == 'EW_z2':
            F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
            F_both_DESI[-1,j] = F_ULTRASAT_DESI[pars.index('EW_z2'),pars.index('EW_z1')]
    F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]
    F_both_DESI[-1,-1] = F_ULTRASAT_DESI[pars.index('EW_z2'),pars.index('EW_z2')]


    for j in range(len(F_GALEX)):
        if pars[j] == 'gamma_1500':
            F_GALEX[j,j] += 1/.3**2
            F_ULTRASAT[j,j] += 1/.3**2
            F_both[j,j] += 1/.3**2
            F_ULTRASAT_DESI[j,j] += 1/.3**2
            F_both_DESI[j,j] += 1/.3**2
        if pars[j] == 'C_alpha_1500':
            F_GALEX[j,j] += 1/1.5**2
            F_ULTRASAT[j,j] += 1/1.5**2
            F_both[j,j] += 1/1.5**2
            F_ULTRASAT_DESI[j,j] += 1/1.5**2
            F_both_DESI[j,j] += 1/1.5**2
        if pars[j] == 'C_alpha_1100':
            F_GALEX[j,j] += 1/1.5**2
            F_ULTRASAT[j,j] += 1/1.5**2
            F_both[j,j] += 1/1.5**2
            F_ULTRASAT_DESI[j,j] += 1/1.5**2
            F_both_DESI[j,j] += 1/1.5**2

    
    inv_F_GALEX = np.linalg.inv(F_GALEX)
    inv_F_ULTRASAT = np.linalg.inv(F_ULTRASAT)
    inv_F_both = np.linalg.inv(F_both)

    inv_F_ULTRASAT_DESI = np.linalg.inv(F_ULTRASAT_DESI)
    inv_F_both_DESI = np.linalg.inv(F_both_DESI)

    print('\nEPS_1500 GALEX = ' + str(round(100 * np.sqrt(inv_F_GALEX[0][0])/fid_vals_det('GALEX_NUV')[pars_fid.index('eps_1500_0')],2)) + ' %')
    print('EPS_1500 ULTRASAT = ' + str(round(100 * np.sqrt(inv_F_ULTRASAT[0][0])/fid_vals_det('ULTRASAT')[pars_fid.index('eps_1500_0')],2)) + ' %')
    print('EPS_1500 both = ' + str(round(100 * np.sqrt(inv_F_both[0][0])/fid_vals_det('GALEX_NUV')[pars_fid.index('eps_1500_0')],2)) + ' %\n')

    print('EPS_1500 ULTRASAT x DESI = ' + str(round(100 * np.sqrt(inv_F_ULTRASAT_DESI[0][0])/fid_vals_det('ULTRASAT')[pars_fid.index('eps_1500_0')],2)) + ' %')
    print('EPS_1500 both DESI = ' + str(round(100 * np.sqrt(inv_F_both_DESI[0][0])/fid_vals_det('GALEX_NUV')[pars_fid.index('eps_1500_0')],2)) + ' %\n')

    if not fix_bias:
        print('b_1500 GALEX = ' + str(round(np.sqrt(inv_F_GALEX[1][1])/fid_vals_det('GALEX_NUV')[pars_fid.index('bias_1500_0')],2)) + ' %')
        print('b_1500 ULTRASAT = ' + str(round(np.sqrt(inv_F_ULTRASAT[1][1])/fid_vals_det('ULTRASAT')[pars_fid.index('bias_1500_0')],2)) + ' %')
        print('b_1500 both = ' + str(round(np.sqrt(inv_F_both[1][1])/fid_vals_det('GALEX_NUV')[pars_fid.index('bias_1500_0')],2)) + ' %\n')

        print('b_1500 ULTRASAT x DESI = ' + str(round(np.sqrt(inv_F_ULTRASAT_DESI[1][1])/fid_vals_det('ULTRASAT')[pars_fid.index('bias_1500_0')],2)) + ' %')
        print('b_1500 both DESI = ' + str(round(np.sqrt(inv_F_both_DESI[1][1])/fid_vals_det('GALEX_NUV')[pars_fid.index('bias_1500_0')],2)) + ' %')

    return


pars_original_c18 = ['log10_epsbias_1500_0', 'gamma_1500', 'alpha_1500_0', 'C_alpha_1500', 'alpha_1100_0', 'C_alpha_1100', 'EW_z1', 'EW_z2', 'gamma_bv', 'gamma_bz']
pars_all = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2','log_fLyC_1','log_fLyC_2','gamma_bv','gamma_bz','alpha_900']


def compare_surveys(detector = 'GALEX_ULTRASAT', pars =pars_original_c18, prior = False, plot_flag = True, group_vox = True):

    #if pars == pars_original_c18:
    #    prior = True

    print('Doing contour plots')

    if detector == 'GALEX':

        F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
        F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
        temp = F_N + F_F

    elif detector == 'GALEX_DESI':

        F_N = Fisher_change_var(pars,'GALEX_NUV','DESI',group_vox,False)
        F_F = Fisher_change_var(pars,'GALEX_FUV','DESI',group_vox,False)
        temp = F_N + F_F

    elif detector == 'GALEX_NUV':

        F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
        temp = F_N

    elif detector == 'GALEX_FUV':

        F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
        temp = F_F


    elif detector == 'ULTRASAT':

        F = Fisher_change_var(pars,'ULTRASAT','SPHEREx',group_vox,False)
        temp = F

    elif detector == 'ULTRASAT_DESI':

        F = Fisher_change_var(pars,'ULTRASAT','DESI',group_vox,False)
        temp = F

    elif detector == 'GALEX_ULTRASAT':

        F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
        F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
        F_GALEX = F_N + F_F
        F_ULTRASAT = Fisher_change_var(pars,'ULTRASAT','SPHEREx',group_vox,False)

        F_both = np.zeros((len(pars)+1,len(pars)+1))
        for i in range(len(pars)):
            if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
            elif pars[i] == 'EW_z1':
                for j in range(len(pars)):
                    F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
            elif pars[i] == 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
        for j in range(len(pars)):
            if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
            elif pars[j] == 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
        F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]

        temp = F_both

    elif detector == 'GALEX_ULTRASAT_DESI':

        F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
        F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
        F_GALEX = F_N + F_F
        F_ULTRASAT = Fisher_change_var(pars,'ULTRASAT','DESI',group_vox,False)

        F_both = np.zeros((len(pars)+1,len(pars)+1))
        for i in range(len(pars)):
            if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
            elif pars[i] == 'EW_z1':
                for j in range(len(pars)):
                    F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
            elif pars[i] == 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
        for j in range(len(pars)):
            if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
            elif pars[j] == 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
        F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]

        temp = F_both

    elif detector == 'GALEX_ULTRASAT_DESIDESI':

        F_N = Fisher_change_var(pars,'GALEX_NUV','DESI',group_vox,False)
        F_F = Fisher_change_var(pars,'GALEX_FUV','DESI',group_vox,False)
        F_GALEX = F_N + F_F
        F_ULTRASAT = Fisher_change_var(pars,'ULTRASAT','DESI',group_vox,False)

        F_both = np.zeros((len(pars)+1,len(pars)+1))
        for i in range(len(pars)):
            if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
            elif pars[i] == 'EW_z1':
                for j in range(len(pars)):
                    F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
            elif pars[i] == 'EW_z2':
                for j in range(len(pars)):
                    if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                        F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                    elif pars[j] == 'EW_z1':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                    elif pars[j] == 'EW_z2':
                        F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
                F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
        for j in range(len(pars)):
            if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
            elif pars[j] == 'EW_z2':
                F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
        F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]

        temp = F_both

    if prior: 
        for j in range(len(pars)):
            if pars[j] == 'gamma_1500':
                temp[j,j] += 1/.3**2
            if pars[j] == 'C_alpha_1500':
                temp[j,j] += 1/1.5**2
            if pars[j] == 'C_alpha_1100':
                temp[j,j] += 1/1.5**2

    inv_F = np.linalg.inv(temp)


    names = []
    fiducials_pars = []    

    if detector == 'GALEX_ULTRASAT' or detector=='GALEX_ULTRASAT_DESI' or detector=='GALEX_ULTRASAT_DESIDESI':
        for i in pars:
            if i == 'log10_epsbias_1500_0':
                fiducials_pars.append(np.log10(fid_vals_det('GALEX_NUV')[pars_fid.index('eps_1500_0')]*fid_vals_det('GALEX_NUV')[pars_fid.index('bias_1500_0')]))
            elif i == 'EW_z1' or i == 'EW_z2':
                fiducials_pars.append(fid_vals_det('GALEX_NUV')[pars_fid.index(i)].value)
            else:
                fiducials_pars.append(fid_vals_det('GALEX_NUV')[pars_fid.index(i)])

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else  r'$\alpha_{900}$' if i == 'alpha_900' else -1
            names.append(name)

        fiducials_pars.append(fid_vals_det('ULTRASAT')[pars_fid.index('EW_z2')].value)
        names.append(r'$EW^{z=2}$')

    elif detector == 'ULTRASAT_DESI':
        for i in pars:
            if i == 'log10_epsbias_1500_0':
                fiducials_pars.append(np.log10(fid_vals_det('ULTRASAT')[pars_fid.index('eps_1500_0')]*fid_vals_det('ULTRASAT')[pars_fid.index('bias_1500_0')]))
            else:
                try:
                    fiducials_pars.append(fid_vals_det('ULTRASAT')[pars_fid.index(i)].value)
                except:
                    fiducials_pars.append(fid_vals_det('ULTRASAT')[pars_fid.index(i)])

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=1}$' if i == 'EW_z1' else r'$EW^{z=2}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else  r'$\alpha_{900}$' if i == 'alpha_900' else -1

            names.append(name)
    else:
        for i in pars:
            if i == 'log10_epsbias_1500_0':
                fiducials_pars.append(np.log10(fid_vals_det(detector)[pars_fid.index('eps_1500_0')]*fid_vals_det(detector)[pars_fid.index('bias_1500_0')]))
            else:
                try:
                    fiducials_pars.append(fid_vals_det(detector)[pars_fid.index(i)].value)
                except:
                    fiducials_pars.append(fid_vals_det(detector)[pars_fid.index(i)])

            name = r'$\log_{10}(\epsilon_{1500}^{z=0}b_{1500}^{z=0})$' if i == 'log10_epsbias_1500_0' else r'$\epsilon_{1500}^{z=0}$' if i == 'eps_1500_0' else r'$\gamma_{1500}$' if i == 'gamma_1500' else r'$\alpha_{1500}^{z=0}$' if i == 'alpha_1500_0' else r'$C_{\alpha_{1500}}$' if i == 'C_alpha_1500' else r'$\alpha_{1100}^{z=0}$' if i == 'alpha_1100_0' else r'$C_{\alpha_{1100}}$' if i == 'C_alpha_1100' else r'$EW^{z=0.3}$' if i == 'EW_z1' else r'$EW^{z=1}$' if i == 'EW_z2' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else r'$\alpha_{900}$' if i == 'alpha_900' else -1

            names.append(name)

    print('DETECTOR: ' + detector)
    for i in range(len(names)):
        print('--- ' + names[i] + ': ' + str(fiducials_pars[i]) + ' +- ' + str(round(np.sqrt(inv_F[i][i]),6)))            

    if plot_flag:

        label_det = r'$\rm GALEX\times SDSS$' if detector == 'GALEX' else r'$\rm NUV\times SDSS$' if detector == 'GALEX_NUV' else r'$\rm FUV\times SDSS$' if detector == 'GALEX_FUV' else  r'$\rm ULTRASAT\times SPHEREx$' if detector == 'ULTRASAT' else  r'$\rm ULTRASAT\times DESI$' if detector == 'ULTRASAT_DESI' else r'$\rm Combined$' if detector == 'GALEX_ULTRASAT' else r'$\rm Combined\, with\,DESI$' if detector == 'GALEX_ULTRASAT_DESI' else -1

        color = color_FUV if detector == 'GALEX' else color_ULTRASAT if (detector == 'ULTRASAT' or detector == 'ULTRASAT_DESI') else color_ULTRASAT if (detector == 'GALEX_ULTRASAT' or detector == 'GALEX_ULTRASAT_DESI') else color_FUV if detector == 'GALEX_FUV' else color_NUV if detector == 'GALEX_NUV' else -1

        line_arg = {'ls':'-', 'color':color} 

        plot_dist_list = []
        plot_dist_list.append(GaussianND(fiducials_pars, inv_F, names=names))

        settings = gp.GetDistPlotSettings()

        settings.norm_prob_label = False
        settings.axes_fontsize = fontsize*.8
        settings.axes_labelsize = fontsize*.8
        settings.legend_fontsize = fontsize*.9
        settings.fontsize = fontsize
        settings.fig_width_inch = 12

        g = gp.get_subplot_plotter(settings=settings)

        g.settings.figure_legend_frame = False
        g.settings.alpha_filled_add=0.4
        g.settings.title_limit_fontsize = 14
        g.settings.linewidth= linewidth

        g.triangle_plot(plot_dist_list, names,  
        colors=[color],line_args=line_arg,\
        legend_labels=[label_det],legend_loc='upper right',\
            markers={'x2':0})

        for i in range(len(pars)-1):
            ax = g.subplots[i+1, 0]
            ax.yaxis.set_label_coords(-0.7, 0.5)

        plt.tight_layout()

        if detector == 'GALEX_ULTRASAT':
            name_plot = 'ellipse_combined' 
        elif detector == 'GALEX_ULTRASAT_DESI':
            name_plot = 'ellipse_combined_DESI' 
        else:
            name_plot = 'ellipse_' + detector
        filename = 'results/PLOTS/EBL/' + name_plot + '.png'
        plt.savefig(filename,bbox_inches='tight')

        plt.show()

    return


def sum_galex_ultrasat_desi(pars, F_N, F_F, F_ULTRASAT):

    F_GALEX = F_N + F_F

    F_both = np.zeros((len(pars)+1,len(pars)+1))
    for i in range(len(pars)):
        if pars[i] != 'EW_z1' and pars[i] != 'EW_z2':
            for j in range(len(pars)):
                if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                    F_both[i,j] = F_ULTRASAT[i,j] + F_GALEX[i,j]
                elif pars[j] == 'EW_z1':
                    F_both[i,j] = F_GALEX[i,pars.index('EW_z1')]                    
                elif pars[j] == 'EW_z2':
                    F_both[i,j] = F_GALEX[i,pars.index('EW_z2')]  +  F_ULTRASAT[i,pars.index('EW_z1')]
            F_both[i,-1] = F_ULTRASAT[i,pars.index('EW_z2')]
        elif pars[i] == 'EW_z1':
            for j in range(len(pars)):
                F_both[i,j] = F_GALEX[pars.index('EW_z1'),j]
        elif pars[i] == 'EW_z2':
            for j in range(len(pars)):
                if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
                    F_both[i,j] = F_ULTRASAT[pars.index('EW_z1'),j] + F_GALEX[pars.index('EW_z2'),j]
                elif pars[j] == 'EW_z1':
                    F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z1')]                    
                elif pars[j] == 'EW_z2':
                    F_both[i,j] = F_GALEX[pars.index('EW_z2'),pars.index('EW_z2')]  +  F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z1')]
            F_both[i,-1] = F_ULTRASAT[pars.index('EW_z1'),pars.index('EW_z2')]
    for j in range(len(pars)):
        if pars[j] != 'EW_z1' and pars[j] != 'EW_z2':
            F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),j]
        elif pars[j] == 'EW_z2':
            F_both[-1,j] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z1')]
    F_both[-1,-1] = F_ULTRASAT[pars.index('EW_z2'),pars.index('EW_z2')]

    return F_both 


###################################################



def plot_err_noninonizing_cont(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = True, prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)


    # required vals: log10(b1500eps1500), gamma1500, alpha1500, C1500

    required_pars = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500']

    fid_log10_beps = np.log10(fiducials['eps1500'][0]*fiducials['bias'][0])
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)
    
    if prior: 
        for j in range(len(use_pars_fid)):
            if use_pars_fid[j] == 'gamma_1500':
                F_N[j,j] += 1/.3**2
                F_F[j,j] += 1/.3**2
                F_U[j,j] += 1/.3**2
            if use_pars_fid[j] == 'C_alpha_1500':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'C_alpha_1100':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2

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

        temp = 10**(fid_log10_beps) * (1+z[i])**fid_gamma1500 * (use_nu / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        fid_eps_1500[i] = temp

        deps_dlog10_beps = temp * np.log(10)
        deps_dgamma = temp * np.log(1+z[i])
        deps_dalpha = temp * np.log(use_nu / nu_1500)
        deps_dC = temp * np.log(use_nu / nu_1500) * np.log10(1+z[i])

        J = np.array((deps_dlog10_beps,deps_dgamma,deps_dalpha,deps_dC))

        sigma_eps1500_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_eps1500_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_eps1500_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))

    if plot_flag:
        plt.plot(z, fid_eps_1500,color='k')
        plt.fill_between(z, fid_eps_1500 - sigma_eps1500_b, fid_eps_1500 + sigma_eps1500_b, color='k', alpha = 0.2, label=r'$\rm GALEX + ULTRASAT$')
        plt.fill_between(z, fid_eps_1500 + sigma_eps1500_b, fid_eps_1500 + sigma_eps1500_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')
        plt.fill_between(z, fid_eps_1500 + sigma_eps1500_U, fid_eps_1500 + sigma_eps1500_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX$')
        plt.fill_between(z, fid_eps_1500 - sigma_eps1500_U, fid_eps_1500 - sigma_eps1500_b, color=color_ULTRASAT, alpha = 0.2,)
        plt.fill_between(z, fid_eps_1500 - sigma_eps1500_G, fid_eps_1500 - sigma_eps1500_U, color=color_FUV, alpha = 0.2, )

        #plt.ylim(1e25,1e27)
        plt.ylabel(r'$b_{1500}^{z=0}\,\epsilon_{%g}$'%lambda_val)
        plt.yscale('log')
        plt.legend(loc=4)
        plt.xlabel(r'$z$')
        #plt.ylim(1e25,9e26)
        plt.xlim(z[0],z[-1])


        if use_pars_fid == pars_original_c18:
            filename = 'results/PLOTS/EBL/eps_nonion_err_wave' + str(1500) + '.png'
        elif use_pars_fid == pars_all:
            filename = 'results/PLOTS/EBL/eps_nonion_err_wave' + str(1500) + '_allpar.png'

        # !!! bias as function of z and nu ??? !!!

        plt.savefig(filename)
        plt.show()

    else:
        return [sigma_eps1500_G, sigma_eps1500_U, sigma_eps1500_b], fid_eps_1500

    return



def plot_err_line(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = True, prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)


    required_pars = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_z1','EW_z2']

    fid_log10_beps = np.log10(fiducials['eps1500'][0]*fiducials['bias'][0])
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    if prior: 
        for j in range(len(use_pars_fid)):
            if use_pars_fid[j] == 'gamma_1500':
                F_N[j,j] += 1/.3**2
                F_F[j,j] += 1/.3**2
                F_U[j,j] += 1/.3**2
            if use_pars_fid[j] == 'C_alpha_1500':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'C_alpha_1100':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2

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

        C = 10**(fid_log10_beps) * (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp_galex = C * ( A_galex + B )

        A_ultrasat =  EW_val(z[i],'ULTRASAT', False) / (0.005*nu_1216) * (use_nu**2 / cu.c.to(u.AA/u.s))
        temp_ultrasat = C * ( A_ultrasat + B )

        fid_eps_galex[i] = temp_galex
        fid_eps_ultrasat[i] = temp_ultrasat

        deps_dlog10_beps_galex = temp_galex * np.log(10)
        deps_dgamma_galex = temp_galex * np.log(1+z[i])
        deps_dalpha1500_galex = temp_galex * np.log(nu_1216 / nu_1500)
        deps_dC1500_galex = temp_galex * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_galex = C * A_galex * np.log(use_nu / nu_1216)
        deps_dC1100_galex = C * A_galex * np.log(use_nu / nu_1216) * np.log10(1+z[i])
        deps_dEW03_galex = C * B * (1 -  np.log10((1+z[i])/(1+0.3)) / np.log10((1+1)/(1+0.3)) )
        deps_dEW1_galex = C * B * np.log10((1+z[i])/(1+0.3)) / np.log10((1+1)/(1+0.3)) 

        deps_dlog10_beps_ultrasat = temp_ultrasat * np.log(10)
        deps_dgamma_ultrasat = temp_ultrasat * np.log(1+z[i])
        deps_dalpha1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500)
        deps_dC1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_ultrasat = C * A_ultrasat * np.log(use_nu / nu_1216)
        deps_dC1100_ultrasat = C * A_ultrasat * np.log(use_nu / nu_1216) * np.log10(1+z[i])
        deps_dEW1_ultrasat = C * B * (1 - np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)))
        deps_dEW2_ultrasat = C * B * np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 


        J_galex = np.array((deps_dlog10_beps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,deps_dEW03_galex,deps_dEW1_galex))

        J_ultrasat = np.array((deps_dlog10_beps_ultrasat,deps_dgamma_ultrasat,deps_dalpha1500_ultrasat,deps_dC1500_ultrasat,deps_dalpha1100_ultrasat,deps_dC1100_ultrasat,deps_dEW1_ultrasat,deps_dEW2_ultrasat))

        J_both = np.array((deps_dlog10_beps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,deps_dEW1_ultrasat,deps_dEW2_ultrasat)) 

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J_galex,inv_F_G,J_galex]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J_ultrasat,inv_F_U,J_ultrasat]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J_both,inv_F_b,J_both]))

    if plot_flag:
        plt.plot(z, fid_eps_ultrasat,color='k')
        plt.fill_between(z, fid_eps_ultrasat - sigma_eps_b, fid_eps_ultrasat + sigma_eps_b, color='k', alpha = 0.2, label=r'$\rm GALEX + ULTRASAT$')
        plt.fill_between(z, fid_eps_ultrasat + sigma_eps_b, fid_eps_ultrasat + sigma_eps_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')
        plt.fill_between(z, fid_eps_galex + sigma_eps_U, fid_eps_galex + sigma_eps_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX$')
        plt.fill_between(z, fid_eps_ultrasat - sigma_eps_U, fid_eps_ultrasat - sigma_eps_b, color=color_ULTRASAT, alpha = 0.2,)
        plt.fill_between(z, fid_eps_galex - sigma_eps_G, fid_eps_ultrasat - sigma_eps_U, color=color_FUV, alpha = 0.2, )

        # !!! the bias now has to be frequency dependent !!!

        #plt.ylim(1e25,1e27)
        plt.ylabel(r'$b_{1500}^{z=0}\,\epsilon_\lambda$')#_{%g}$'%lambda_val)
        plt.yscale('log')
        plt.legend(loc=4)
        plt.xlabel(r'$z$')
        #plt.ylim(1e25,9e26)
        plt.xlim(z[0],z[-1])

        if use_pars_fid == pars_original_c18:
            filename = 'results/PLOTS/EBL/eps_line_err_wave' + str(1216) + '.png'
        elif use_pars_fid == pars_all:
            filename = 'results/PLOTS/EBL/eps_line_err_wave' + str(1216) + '_allpar.png'

        plt.savefig(filename)
        plt.show()

    else:
        return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps_ultrasat

    return


def plot_err_line_surrounding(z, lambda_val, use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = True,prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)


    required_pars = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100']

    fid_log10_beps = np.log10(fiducials['eps1500'][0]*fiducials['bias'][0])
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    if prior: 
        for j in range(len(use_pars_fid)):
            if use_pars_fid[j] == 'gamma_1500':
                F_N[j,j] += 1/.3**2
                F_F[j,j] += 1/.3**2
                F_U[j,j] += 1/.3**2
            if use_pars_fid[j] == 'C_alpha_1500':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'C_alpha_1100':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2

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

        C = 10**(fid_log10_beps) * (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp_galex = C * ( B )

        temp_ultrasat = C * ( B )

        fid_eps_galex[i] = temp_galex
        fid_eps_ultrasat[i] = temp_ultrasat

        deps_dlog10_beps_galex = temp_galex * np.log(10)
        deps_dgamma_galex = temp_galex * np.log(1+z[i])
        deps_dalpha1500_galex = temp_galex * np.log(nu_1216 / nu_1500)
        deps_dC1500_galex = temp_galex * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_galex = C  * np.log(use_nu / nu_1216)
        deps_dC1100_galex = C  * np.log(use_nu / nu_1216) * np.log10(1+z[i])

        deps_dlog10_beps_ultrasat = temp_ultrasat * np.log(10)
        deps_dgamma_ultrasat = temp_ultrasat * np.log(1+z[i])
        deps_dalpha1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500)
        deps_dC1500_ultrasat = temp_ultrasat * np.log(use_nu / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100_ultrasat = C * np.log(use_nu / nu_1216)
        deps_dC1100_ultrasat = C  * np.log(use_nu / nu_1216) * np.log10(1+z[i])


        J_galex = np.array((deps_dlog10_beps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex,))

        J_ultrasat = np.array((deps_dlog10_beps_ultrasat,deps_dgamma_ultrasat,deps_dalpha1500_ultrasat,deps_dC1500_ultrasat,deps_dalpha1100_ultrasat,deps_dC1100_ultrasat))

        J_both = np.array((deps_dlog10_beps_galex,deps_dgamma_galex,deps_dalpha1500_galex,deps_dC1500_galex,deps_dalpha1100_galex,deps_dC1100_galex)) 

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J_galex,inv_F_G,J_galex]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J_ultrasat,inv_F_U,J_ultrasat]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J_both,inv_F_b,J_both]))

    if plot_flag:
        plt.plot(z, fid_eps_ultrasat,color='k')
        plt.fill_between(z, fid_eps_ultrasat - sigma_eps_b, fid_eps_ultrasat + sigma_eps_b, color='k', alpha = 0.2, label=r'$\rm GALEX + ULTRASAT$')
        plt.fill_between(z, fid_eps_ultrasat + sigma_eps_b, fid_eps_ultrasat + sigma_eps_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')
        plt.fill_between(z, fid_eps_galex + sigma_eps_U, fid_eps_galex + sigma_eps_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX$')
        plt.fill_between(z, fid_eps_ultrasat - sigma_eps_U, fid_eps_ultrasat - sigma_eps_b, color=color_ULTRASAT, alpha = 0.2,)
        plt.fill_between(z, fid_eps_galex - sigma_eps_G, fid_eps_ultrasat - sigma_eps_U, color=color_FUV, alpha = 0.2, )

        # !!! the bias now has to be frequency dependent !!!

        #plt.ylim(1e25,1e27)
        plt.ylabel(r'$b_{1500}^{z=0}\,\epsilon_\lambda$')#_{%g}$'%lambda_val)
        plt.yscale('log')
        plt.legend(loc=4)
        plt.xlabel(r'$z$')
        #plt.ylim(1e25,9e26)
        plt.xlim(z[0],z[-1])

        if use_pars_fid == pars_original_c18:
            filename = 'results/PLOTS/EBL/eps_line_err_wave' + str(1216) + '.png'
        elif use_pars_fid == pars_all:
            filename = 'results/PLOTS/EBL/eps_line_err_wave' + str(1216) + '_allpar.png'

        plt.savefig(filename)
        plt.show()

    else:
        return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps_ultrasat

    return


def plot_err_ioncont(z, lambda_val, use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = True, prior=False,galex_detector='SDSS'):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)
    nu_1216 = nu_from_lambda(1216*u.AA)
    nu_912 = nu_from_lambda(912*u.AA)

    if use_pars_fid == pars_original_c18:
        required_pars = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100']

    else:
        required_pars = ['log10_epsbias_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','alpha_900','log_fLyC_1','log_fLyC_2']

    fid_log10_beps = np.log10(fiducials['eps1500'][0]*fiducials['bias'][0])
    fid_gamma1500 = fiducials['eps1500'][1]
    fid_alpha1500_0 = fiducials['alpha1500'][0]
    fid_C1500 = fiducials['alpha1500'][1]
    fid_alpha1100_0 = fiducials['alpha1100'][0]
    fid_C1100 = fiducials['alpha1100'][1]
    fid_alpha900 = fiducials['alpha900']

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV',galex_detector,group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV',galex_detector,group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    if prior: 
        for j in range(len(use_pars_fid)):
            if use_pars_fid[j] == 'gamma_1500':
                F_N[j,j] += 1/.3**2
                F_F[j,j] += 1/.3**2
                F_U[j,j] += 1/.3**2
            if use_pars_fid[j] == 'C_alpha_1500':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'C_alpha_1100':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2

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

        C = 10**(fid_log10_beps) * (1+z[i])**fid_gamma1500 * (nu_1216 / nu_1500)**(fid_alpha1500_0 + fid_C1500*np.log10(1+z[i]))

        temp = C * ( A * B )
        fid_eps[i] = temp 

        deps_dlog10_beps = temp * np.log(10)
        deps_dgamma = temp * np.log(1+z[i])
        deps_dalpha1500 = temp * np.log(nu_1216 / nu_1500)
        deps_dC1500 = temp * np.log(nu_1216 / nu_1500) * np.log10(1+z[i])
        deps_dalpha1100 = temp * np.log(nu_912 / nu_1216)
        deps_dC1100 = temp * np.log(nu_912 / nu_1216) * np.log10(1+z[i])
        if not use_pars_fid == pars_original_c18:
            deps_dlogfLy1 = C * B * (1 -  np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) )
            deps_dlogfLy2 = C * B * np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 
            deps_dalpha900 = temp * np.log(use_nu / nu_912)

            J = np.array((deps_dlog10_beps,deps_dgamma,deps_dalpha1500,deps_dC1500,deps_dalpha1100,deps_dC1100,deps_dalpha900,deps_dlogfLy1,deps_dlogfLy2))
        else:
            J = np.array((deps_dlog10_beps,deps_dgamma,deps_dalpha1500,deps_dC1500,deps_dalpha1100,deps_dC1100))

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))

    if plot_flag:
        plt.plot(z, fid_eps,color='k')
        plt.fill_between(z, fid_eps - sigma_eps_b, fid_eps + sigma_eps_b, color='k', alpha = 0.2, label=r'$\rm GALEX + ULTRASAT$')
        plt.fill_between(z, fid_eps + sigma_eps_b, fid_eps + sigma_eps_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')
        plt.fill_between(z, fid_eps + sigma_eps_U, fid_eps + sigma_eps_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX$')
        plt.fill_between(z, fid_eps - sigma_eps_U, fid_eps - sigma_eps_b, color=color_ULTRASAT, alpha = 0.2,)
        plt.fill_between(z, fid_eps - sigma_eps_G, fid_eps - sigma_eps_U, color=color_FUV, alpha = 0.2, )

        # !!! the bias now has to be frequency dependent !!!

        #plt.ylim(1e25,1e27)
        plt.ylabel(r'$b_{1500}^{z=0}\,\epsilon_\lambda$')#_{%g}$'%lambda_val)
        plt.yscale('log')
        plt.legend(loc=4)
        plt.xlabel(r'$z$')
        #plt.ylim(1e25,9e26)
        plt.xlim(z[0],z[-1])

        if use_pars_fid == pars_original_c18:
            filename = 'results/PLOTS/EBL/eps_ion_err_wave' + str(912) + '.png'
        elif use_pars_fid == pars_all:
            filename = 'results/PLOTS/EBL/eps_ion_err_wave' + str(912) + '_allpar.png'

        plt.savefig(filename)
        plt.show()

    else:
        return [sigma_eps_G, sigma_eps_U, sigma_eps_b], fid_eps

    return





def plot_err_escape(group_vox=True,run=False,plot_flag = True):

    z = np.linspace(zmin_gal,zmax_gal,200)

    required_pars = ['log_fLyC_1','log_fLyC_2']

    use_pars_fid = pars_all
    F_N_c18 = Fisher_change_var(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F_c18 = Fisher_change_var(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV','DESI',group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV','DESI',group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    for j in range(len(use_pars_fid)):
        if use_pars_fid[j] == 'gamma_1500':
            F_N_c18[j,j] += 1/.3**2
            F_F_c18[j,j] += 1/.3**2
        if use_pars_fid[j] == 'C_alpha_1500':
            F_N_c18[j,j] += 1/1.5**2
            F_F_c18[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_1100':
            F_N_c18[j,j] += 1/1.5**2
            F_F_c18[j,j] += 1/1.5**2


    F_galex = F_N_c18 + F_F_c18
    F_both = sum_galex_ultrasat_desi(use_pars_fid, F_N, F_F, F_U)

    all_inv_F_G = np.linalg.inv(F_galex)
    #all_inv_F_U = np.linalg.inv(F_U)
    all_inv_F_b = np.linalg.inv(F_both)

    inv_F_G = np.zeros((len(required_pars),len(required_pars)))
    #inv_F_U = np.zeros((len(required_pars),len(required_pars)))
    inv_F_b = np.zeros((len(required_pars),len(required_pars)))

    for i in range(len(required_pars)):
       id_i = use_pars_fid.index(required_pars[i])
       for j in range(len(required_pars)):
           id_j = use_pars_fid.index(required_pars[j])
           inv_F_G[i,j] = all_inv_F_G[id_i,id_j]
           #inv_F_U[i,j] = all_inv_F_U[id_i,id_j]
           inv_F_b[i,j] = all_inv_F_b[id_i,id_j]
#
    fid_fLy = np.zeros(len(z))

    sigma_eps_G = np.zeros(len(z))
    #sigma_eps_U = np.zeros(len(z))
    sigma_eps_b = np.zeros(len(z))

    for i in range(len(z)):

        fid_fLy[i] =  fLyC(z[i], False) 

        deps_dlogfLy1 = (1 -  np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) )
        deps_dlogfLy2 = np.log10((1+z[i])/(1+1)) / np.log10((1+2)/(1+1)) 

        J = np.array((deps_dlogfLy1,deps_dlogfLy2))

        sigma_eps_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        #sigma_eps_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_eps_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))


    plt.plot(z, fid_fLy,color='k')
    plt.fill_between(z, fid_fLy - sigma_eps_b, fid_fLy + sigma_eps_b, color=color_ULTRASAT, alpha = 0.5, label=r'$\rm This\,work,\,1\sigma$')
    plt.fill_between(z, fid_fLy - 3*sigma_eps_b, fid_fLy + 3*sigma_eps_b, color=color_ULTRASAT, alpha = 0.1, label=r'$\rm This\,work,\,3\sigma$')

    plt.plot(z, fid_fLy - sigma_eps_G,color_FUV,alpha=0.6,label=r'$\rm GALEX\times SDSS,\,1\sigma\,with\,prior$')
    plt.plot(z, fid_fLy + sigma_eps_G,color_FUV,alpha=0.6,)

    plt.ylabel(r'$f_{\rm LyC}$',fontsize=fontsize)#_{%g}$'%lambda_val)
    plt.legend(loc=1)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])
    plt.ylim(0.1,0.9)

    print(max(3*sigma_eps_b/fid_fLy))

    filename = 'results/PLOTS/EBL/errfLyC.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return


def plot_err_bias(lambda_val, use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = True,prior=False):

    use_nu = nu_from_lambda(lambda_val*u.AA)
    nu_1500 = nu_from_lambda(1500*u.AA)

    z = np.linspace(zmin_gal,zmax_gal,200)

    required_pars = ['log10_epsbias_1500_0','gamma_bv','gamma_bz']

    fid_log10_beps = np.log10(fiducials['eps1500'][0]*fiducials['bias'][0])
    fid_gamma_bnu = fiducials['bias'][1]
    fid_gamma_bz = fiducials['bias'][2]

    F_N = Fisher_change_var(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F = Fisher_change_var(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_U = Fisher_change_var(use_pars_fid,'ULTRASAT','DESI',group_vox,run)

    if prior: 
        for j in range(len(use_pars_fid)):
            if use_pars_fid[j] == 'gamma_1500':
                F_N[j,j] += 1/.3**2
                F_F[j,j] += 1/.3**2
                F_U[j,j] += 1/.3**2
            if use_pars_fid[j] == 'C_alpha_1500':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2
            if use_pars_fid[j] == 'C_alpha_1100':
                F_N[j,j] += 1/1.5**2
                F_F[j,j] += 1/1.5**2
                F_U[j,j] += 1/1.5**2


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
    fid_bias = np.zeros(len(z))

    sigma_bias_G = np.zeros(len(z))
    sigma_bias_U = np.zeros(len(z))
    sigma_bias_b = np.zeros(len(z))

    for i in range(len(z)):

        fid_bias[i] =  pow(10,fid_log10_beps) / 1e25 *(1+z[i])**fid_gamma_bz * (use_nu / nu_1500)**fid_gamma_bnu 

        db_dlogepsb = fid_bias[i] * np.log(10)
        db_dgamma_bnu = fid_bias[i] * np.log(use_nu / nu_1500)
        db_dgamma_bz = fid_bias[i] * np.log((1+z[i]))

        J = np.array((db_dlogepsb,db_dgamma_bnu,db_dgamma_bz))

        sigma_bias_G[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_G,J]))
        sigma_bias_U[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_U,J]))
        sigma_bias_b[i] = np.sqrt(np.linalg.multi_dot([J,inv_F_b,J]))

    plt.plot(z, fid_bias,color='k')
    plt.fill_between(z, fid_bias - sigma_bias_b, fid_bias + sigma_bias_b, color='k', alpha = 0.2, label=r'$\rm GALEX + ULTRASAT$')
    plt.fill_between(z, fid_bias + sigma_bias_b, fid_bias + sigma_bias_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')
    plt.fill_between(z, fid_bias + sigma_bias_U, fid_bias + sigma_bias_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX$')
    plt.fill_between(z, fid_bias - sigma_bias_U, fid_bias - sigma_bias_b, color=color_ULTRASAT, alpha = 0.2,)
    plt.fill_between(z, fid_bias - sigma_bias_G, fid_bias - sigma_bias_U, color=color_FUV, alpha = 0.2, )

    # !!! the bias now has to be frequency dependent !!!

    #plt.ylim(1e25,1e27)
    plt.ylabel(r'$b(\nu,z)\epsilon^{z=0}_{1500}/10^{25}$')#_{%g}$'%lambda_val)
    #plt.yscale('log')
    plt.legend(loc=1)
    plt.xlabel(r'$z$')
    plt.ylim(0,5)
    plt.xlim(z[0],z[-1])

    if use_pars_fid == pars_original_c18:
        filename = 'results/PLOTS/EBL/errbias_' + str(lambda_val)+ '.png'
    elif use_pars_fid == pars_all:
        filename = 'results/PLOTS/EBL/errbias_' + str(lambda_val)+ '_allpar.png'

    plt.savefig(filename)
    plt.show()

    return [sigma_bias_G, sigma_bias_U, sigma_bias_b], fid_bias


def plot_multi_line_PAPER():

    z = np.linspace(zmin_gal,zmax_gal+.5,200)

    sigmas_nonion, fid_nonion = plot_err_noninonizing_cont(z, 1500, use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False,galex_detector='SDSS',prior=True)

    sigmas_nonion_nonp = plot_err_noninonizing_cont(z, 1500, use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False,galex_detector='SDSS')[0]

    sigmas_nonion_both = plot_err_noninonizing_cont(z, 1500, use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = False,galex_detector='DESI')[0]

    s_G = sigmas_nonion[0] / fiducials['bias'][0] #* (70/67.67)**3
    s_G_np = sigmas_nonion_nonp[0] / fiducials['bias'][0]#* (70/67.67)**3
    s_b = sigmas_nonion_both[2] / fiducials['bias'][0]#* (70/67.67)**3
    s_3b = 3*sigmas_nonion_both[2] / fiducials['bias'][0]#* (70/67.67)**3

    fid_nonion /= fiducials['bias'][0] 
    #fid_nonion *= (70/67.67)**3

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

    plt.figure(figsize=(19,10))

    #plt.errorbar([0.055], [25.54], yerr=[[0.02],[0.09]], fmt='o', capsize=5, markersize=7, color='k', label=r'$\rm Wyder\, et\, al.,\,2005$',alpha=.5,markerfacecolor='none')
    
    plt.errorbar([(1.9+2.7)/2], [np.log10(3.63e26)], yerr=[0.4], fmt='*', capsize=5, markersize=7, color='g', label=r'$\rm Reddy\, et\, al.,\,2008$',alpha=.7,markerfacecolor='none')
    
    plt.errorbar(z_data, np.log10(eps_1500_data), yerr=[np.log10(eps_err_down), np.log10(eps_err_up)], fmt='o', capsize=5, markersize=15, color='indigo', label=r'$\rm Schiminovich\, et\, al.,\,2005$',alpha=.7,markerfacecolor='none')
    
    plt.errorbar(z_data_hst, np.log10(eps_hst), yerr=np.log10(eps_err_hst), fmt='D', capsize=5, markersize=15, color='b', label=r'$\rm Dahlen\, et\, al.,\,2006$',alpha=.7,markerfacecolor='none')

    plt.fill_between(z, np.log10(fid_nonion - s_b), np.log10(fid_nonion + s_b), color=color_ULTRASAT, alpha = 0.5,label=r'$\rm This\,work,\,1\sigma$')

    plt.fill_between(z, np.log10(fid_nonion - s_3b), np.log10(fid_nonion + s_3b), color=color_ULTRASAT, alpha = 0.1,label=r'$\rm This\,work,\,3\sigma$')

    plt.plot(z_d11, np.log10(eps_d11), 'k-.', label=r'$\rm Dominguez\, et\, al.,\,2011$',alpha=.6,linewidth=1.9)
    
    plt.plot(z_gil, np.log10(eps_gil), 'k--', label=r'$\rm Gilmore\, et\, al.,\,2011$',alpha=.6,linewidth=1.9)
    
    plt.plot(z_gil2, np.log10(eps_gil2), 'k--',alpha=.6,linewidth=1.9)
    
    plt.plot(z, np.log10(fid_nonion),'k')

    #plt.plot(z, np.log10(fid_nonion - s_G),color_FUV,alpha=0.6,label=r'$\rm C18\,results,\,(1\sigma)$',linewidth=2)
    #plt.plot(z, np.log10(fid_nonion + s_G),color_FUV,alpha=0.6,linewidth=2)

    # plt.plot(z, np.log10(fid_nonion - s_G_np),color_FUV,alpha=0.6,label=r'$\rm C18\,results,\,(1\sigma\,without\,priors)$',linewidth=2)
    # plt.plot(z, np.log10(fid_nonion + s_G_np),color_FUV,alpha=0.6,linewidth=2)

    print(np.max(s_b/fid_nonion))

    plt.ylabel(r'$\log_{10}\epsilon_{1500}(z)$',fontsize=fontsize)
    #plt.yscale('log')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim(25.5,27.)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])


    filename = 'results/PLOTS/EBL/eps_err_multiwave_PAPER.png'
 
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()


    return 









def plot_multi_line(detector, use_pars_fid = pars_original_c18, line = True, nonion_cont = True, ion_cont = True, DM_contr = False,galex_detector='DESI'):

    z = np.linspace(zmin_gal,zmax_gal,200)

    if line:
        sigmas_line, fid_line = plot_err_line(z, 1216, use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,galex_detector=galex_detector)
    if nonion_cont:
        sigmas_nonion, fid_nonion = plot_err_noninonizing_cont(z, 1500, use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,galex_detector=galex_detector)
    if ion_cont and use_pars_fid == pars_all:
        sigmas_ion, fid_ion = plot_err_ioncont(z, 912, use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,galex_detector=galex_detector)

    if detector == 'GALEX':
        if line:
            use_sigma_line = sigmas_line[0]
        if nonion_cont:
            use_sigma_nonion = sigmas_nonion[0]
        if ion_cont and use_pars_fid == pars_all:
            use_sigma_ion = sigmas_ion[0]
        use_color = color_NUV
        label_detector = r'$\rm GALEX \times SDSS$'

        use_z_long = np.linspace(z[0],z[-1],500)
        z_line_detected_fuv = []
        z_line_detected_nuv = []
        for i in range(len(use_z_long)):
            if (1216*(1+use_z_long[i]) < lambda_from_nu(nu_min_gFUV).value) and (1216*(1+use_z_long[i]) > lambda_from_nu(nu_max_gFUV).value) :
                z_line_detected_fuv.append(use_z_long[i])
            if (1216*(1+use_z_long[i]) < lambda_from_nu(nu_min_gNUV).value) and (1216*(1+use_z_long[i]) > lambda_from_nu(nu_max_gNUV).value) :
               z_line_detected_nuv.append(use_z_long[i])

        z_line_detected_nuv = np.asarray(z_line_detected_nuv)
        z_line_detected_fuv = np.asarray(z_line_detected_fuv)

    elif detector == 'ULTRASAT':
        if line:
            use_sigma_line = sigmas_line[1]
        if nonion_cont:
            use_sigma_nonion = sigmas_nonion[1]
        if ion_cont and use_pars_fid == pars_all:
            use_sigma_ion = sigmas_ion[1]
        use_color = color_ULTRASAT
        label_detector = r'$\rm ULTRASAT \times DESI$'

    elif detector == 'GALEX_ULTRASAT':
        if line:
            use_sigma_line = sigmas_line[2]
        if nonion_cont:
            use_sigma_nonion = sigmas_nonion[2]
        if ion_cont and use_pars_fid == pars_all:
            use_sigma_ion = sigmas_ion[2]
        use_color = 'k'
        label_detector = r'$\rm (GALEX + ULTRASAT) \times DESI)$'

    if ion_cont and use_pars_fid == pars_all:
        plt.figure(figsize=(20,10))
        plt.plot(z, fid_ion,'k:',label=r'$\lambda = 912{\rm \AA}$')
    if line:
        plt.plot(z, fid_line,color='k',label=r'$\lambda = 1216{\rm \AA}$')
    if nonion_cont:
        plt.plot(z, fid_nonion,'k--',label=r'$\lambda = 1500{\rm \AA}$')

    if line:
        plt.fill_between(z, fid_line - use_sigma_line, fid_line + use_sigma_line, color=use_color, alpha = 0.4, )
    if nonion_cont:
        plt.fill_between(z, fid_nonion - use_sigma_nonion, fid_nonion + use_sigma_nonion, color=use_color, alpha = 0.4,label=label_detector)
    if ion_cont and use_pars_fid == pars_all:
        plt.fill_between(z, fid_ion - use_sigma_ion, fid_ion + use_sigma_ion, color=use_color, alpha = 0.4)
    
    #plt.fill_between(z[0],min(z_line_detected_fuv),'k',alpha=0.2) 
    #plt.fill_between(max(z_line_detected_fuv),z[-1],'k',alpha=0.2)
#
    #plt.fill_between(z[0],min(z_line_detected_nuv),color_NUV,alpha=0.2)
    #plt.fill_between(max(z_line_detected_nuv),z[-1],color_NUV,alpha=0.2)

    if DM_contr:
        eps_DM = np.zeros(len(z))
        Gamma = 7e-40*u.s**-1 #H(0.).to(u.s**-1)
        rho_crit = (3*H(0.)**2 / (8*np.pi*cu.G)).to(u.Msun/u.Mpc**3)
        Omega_c = (camb_pars.omch2 / (H(0.).value/100)**2)
        conf_time = lambda zv: (cosmo.conformal_time(zv)*u.Mpc/cu.c).to(u.s)
        proper_time = lambda zv: conf_time(zv) / (1+zv)
        intrinsic_eps_DM = lambda zv: (Gamma * rho_crit * Omega_c * cu.c**2 * (1+zv)**3 *np.exp(-Gamma* (proper_time(zv)-proper_time(0.)))*u.Hz**-1).to(u.erg/u.s/u.Mpc**3/u.Hz)

        for zi in range(len(z)):
            eps_DM[zi] = intrinsic_eps_DM(z[zi]).value

        print(eps_DM* fiducials['bias'][0],)
        plt.plot(z,eps_DM * fiducials['bias'][0],'r',label=r'$\rm Decaying\, DM,\,\Theta_{\rm DM} = 7\times 10^{-40}\,s$')



    plt.ylabel(r'$\varepsilon(\nu,z)$',fontsize=fontsize)
    plt.yscale('log')
    if ion_cont and use_pars_fid == pars_all:
        plt.legend(bbox_to_anchor=(1,0.7))
    else:
        plt.legend(loc=4)
    plt.ylim(1e24,1e27)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.xlim(z[0],z[-1])


    if use_pars_fid == pars_original_c18:
        filename = 'results/PLOTS/EBL/eps_err_multiwave' + detector + '.png'
    elif use_pars_fid == pars_all:
        filename = 'results/PLOTS/EBL/eps_err_multiwave' + detector + '_allpars.png'

    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()


    return 


##############################
##############################
##############################
##############################



def plot_err_emissivity(z = [0.], use_pars_fid = pars_all, prior=False,plot_flag=True):

    
    wave = np.linspace(700,3000,500)
    signal_val = np.zeros(len(wave))

    sigma_galex = np.zeros(len(wave))
    sigma_ultrasat = np.zeros(len(wave))
    sigma_both = np.zeros(len(wave))
    for w in range(len(wave)):

        signal_val[w] =  signal(wave[w]*u.AA,z[0],'ULTRASAT',vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False).value

        if wave[w] <= 912:
            temp = plot_err_ioncont(z,wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False, prior=prior)[0]
            
        elif 912 < wave[w] < 1216*(1-0.005):
            temp = plot_err_line_surrounding(z,wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior)[0]

        elif 1216*(1-0.005) < wave[w] <= 1216:
              temp = plot_err_line(z,wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior)[0]
  
        else:
            temp = plot_err_noninonizing_cont(z,wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior)[0]

        sigma_galex[w] = temp[0][0] / fiducials['bias'][0]
        sigma_ultrasat[w] = temp[1][0]  / fiducials['bias'][0]
        sigma_both[w] = temp[2][0]  / fiducials['bias'][0]

    if plot_flag:

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,gridspec_kw={'hspace': 0,'height_ratios': [2, 1]})

        ax1.plot(wave, signal_val,'k')


        ax1.fill_between(wave, signal_val - sigma_both, signal_val + sigma_both, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm This\,work$')
        ax1.fill_between(wave, signal_val + sigma_both, signal_val + sigma_galex, color=color_FUV, alpha = 0.2, label=r'$\rm Results\, from\,C18\,(reproduced)$')
        ax1.fill_between(wave, signal_val - sigma_galex, signal_val - sigma_both, color=color_FUV, alpha = 0.2, )

        ax1.set_ylabel(r'$\epsilon(\nu,z)$',fontsize=fontsize)
        ax1.set_yscale('log')
        # ax2.set_xscale('log')
        ax1.legend(loc=1,)
        ax2.set_xlim(wave[0],wave[-1])
        ax1.set_ylim(3e24,3e28)

        ax2.plot(wave, (sigma_ultrasat - sigma_galex)/sigma_ultrasat, color=color_ULTRASAT, label=r'$\sigma_{\rm ULTRASAT} - \sigma_{\rm GALEX}$')

        ax2.axhline(y=0,color='k')

        # Set labels and legend for the lower panel
        ax2.set_ylabel(r'$1-\sigma_{\rm G\times S}/\sigma_{\rm U\times D}$')

        ax1.set_title(r'$z=%g$'%z[0],)
        ax2.set_xlabel(r'$\lambda_{\rm obs}\,[{\rm A}]$',fontsize=fontsize)
        #plt.xticks(fontsize=2*fontsize*.8)
        #plt.yticks(fontsize=2*fontsize*.8)

        if use_pars_fid == pars_original_c18:
            filename = 'results/PLOTS/EBL/err_rest_emissivity_z' + str(z) + '.png'
        elif use_pars_fid == pars_all:
            filename = 'results/PLOTS/EBL/err_rest_emissivity_allpars_z' + str(z) + '.png'

        plt.tight_layout()
        plt.savefig(filename,bbox_inches='tight')
        plt.show()


    return [sigma_galex, sigma_ultrasat, sigma_both], wave, signal_val 


def plot_err_emissivity_PAPER():

    z = [0.,0.3,0.5,0.6,1.,1.5,2.]
    wave = np.linspace(700,3000,500)

    for zi in z:
        signal_val = np.zeros(len(wave))

        sigma_galex = np.zeros(len(wave))
        sigma_both = np.zeros(len(wave))
        for w in range(len(wave)):

            signal_val[w] =  signal(wave[w]*u.AA,zi,'ULTRASAT',vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False).value

            if wave[w] <= 912:
                sigma_galex[w] = plot_err_ioncont([zi],wave[w], use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False, prior=True,galex_detector='SDSS')[0][0][0]
                
                sigma_both[w] = plot_err_ioncont([zi],wave[w], use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = False, prior=False,galex_detector='DESI')[0][2][0]

            elif 912 < wave[w] < 1216*(1-0.005):
                sigma_galex[w] = plot_err_line_surrounding([zi],wave[w], use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False,prior=True,galex_detector='SDSS')[0][0][0]

                sigma_both[w] = plot_err_line_surrounding([zi],wave[w], use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = False,prior=False,galex_detector='DESI')[0][2][0]

            elif 1216*(1-0.005) < wave[w] <= 1216:
                sigma_galex[w] = plot_err_line([zi],wave[w], use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False,prior=True,galex_detector='SDSS')[0][0][0]
    
                sigma_both[w] = plot_err_line([zi],wave[w], use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = False,prior=False,galex_detector='DESI')[0][2][0]
    
            else:
                sigma_galex[w] = plot_err_noninonizing_cont([zi],wave[w], use_pars_fid = pars_original_c18,group_vox=True,run=False,plot_flag = False,prior=True,galex_detector='SDSS')[0][0][0]
                
                sigma_both[w] = plot_err_noninonizing_cont([zi],wave[w], use_pars_fid = pars_all,group_vox=True,run=False,plot_flag = False,prior=False,galex_detector='DESI')[0][2][0]

            sigma_galex[w] /= fiducials['bias'][0]
            sigma_both[w] /= fiducials['bias'][0]
            
            #sigma_galex[w] *= 3
            #sigma_both[w] *= 3


        plt.figure(figsize=(12,9))

        # plt.axvline(lambda_from_nu(nu_min_US).value/(1+zi), color='k',linewidth=.5)
        plt.axvline(lambda_from_nu(nu_min_gNUV).value/(1+zi), color='k',linewidth=1,linestyle=':')
        # plt.axvline(lambda_from_nu(nu_max_gFUV).value/(1+zi), color='k',linewidth=.5)

        plt.fill_betweenx([1e25,5e28],lambda_from_nu(nu_min_US).value/(1+zi),3000,color='k',alpha=0.1)
        plt.fill_betweenx([1e25,5e28],lambda_from_nu(nu_max_gFUV).value/(1+zi),color='k',alpha=0.1)

        plt.plot(wave, signal_val,'k')
        #plt.plot(wave, signal_val - sigma_both, color=color_ULTRASAT, linewidth=1., alpha=0.9, label=r'$\rm This\,work,\,1\sigma$')
        #plt.plot(wave, signal_val + sigma_both, color=color_ULTRASAT, alpha=0.9, linewidth=1.)
        plt.fill_between(wave, signal_val - 3*sigma_both, signal_val + 3*sigma_both, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm This\,work,\,3\sigma$')
        plt.plot(wave, signal_val - sigma_galex,color_FUV,alpha=0.6,label=r'$\rm C18\,results\,(1\sigma,\,reproduced)$',linestyle=(0,(5,1)),linewidth=1.5)
        plt.plot(wave, signal_val + sigma_galex,color_FUV,alpha=0.6,linestyle=(0,(5,1)),linewidth=1.5)

        #plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_max_US).value,color=color_ULTRASAT,alpha=0.1)
        #plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_min_US).value,wave[-2]*(1+zall[-1]),color=color_ULTRASAT,alpha=0.1)


        plt.ylabel(r'$\epsilon(\nu,z)$',fontsize=fontsize)
        plt.yscale('log')
        # ax2.set_xscale('log')
        if zi == 2.:
            plt.legend(loc=1,)
            plt.ylim(1e25,3e28)
        elif zi == 1. or zi == 1.5:
            plt.legend(loc=1,)
            plt.ylim(1e25,3e28)
        # elif zi == 0.:
        #     yticks = [1e25,2e25,3e25,4e25,5e25,6e25,7e25]
        #     yticklabels = [r'$1$',r'$2$',r'$3$',r'$4$',r'$5$',r'$6$',r'$7$']
        #     plt.yticks(yticks, yticklabels)  
        #     plt.ylabel(r'$10^{25}\times\epsilon(\nu,z)$',fontsize=fontsize)
        #     plt.legend(loc=4,)
        #     plt.ylim(2e25,7e25)
        else:
            plt.legend(loc=1,)
            plt.ylim(1e25,1e27)

        plt.xlim(wave[0],wave[-1])

        plt.title(r'$z=%g$'%zi,)
        plt.xlabel(r'$\lambda_{\rm obs}/(1+z)\,[{\rm \AA}]$',fontsize=fontsize)

        filename = 'results/PLOTS/EBL/err_final_emissivity_z' + str(zi) + '.png'
        
        plt.tight_layout()
        plt.savefig(filename,bbox_inches='tight')
        plt.show()


    return 



def plot_err_emissivity_combined(zall = [0.,0.5,1.,2.], use_pars_fid = pars_all, prior=False,plot_flag=True,galex_detector='SDSS'):

    linestyle=['-','--','-.',':']
    plt.figure(figsize=(15,9))
    for z in zall:
        wave = np.linspace(700,3000,500)
        signal_val = np.zeros(len(wave))

        sigma_both = np.zeros(len(wave))
        for w in range(len(wave)):

            signal_val[w] =  signal(wave[w]*u.AA,z,'ULTRASAT',vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False).value

            if wave[w] <= 912:
                temp = plot_err_ioncont([z],wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False, prior=prior,galex_detector=galex_detector)[0]
                
            elif 912 < wave[w] < 1216*(1-0.005):
                temp = plot_err_line_surrounding([z],wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior,galex_detector=galex_detector)[0]

            elif 1216*(1-0.005) < wave[w] <= 1216:
                temp = plot_err_line([z],wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior,galex_detector=galex_detector)[0]
    
            else:
                temp = plot_err_noninonizing_cont([z],wave[w], use_pars_fid = use_pars_fid,group_vox=True,run=False,plot_flag = False,prior=prior,galex_detector=galex_detector)[0]

            sigma_both[w] = temp[2][0] / fiducials['bias'][0]

        plt.plot(wave*(1+z), signal_val,'k',linestyle=linestyle[zall.index(z)],label=r'$z=%g$'%z)

        plt.fill_between(wave*(1+z), signal_val - sigma_both, signal_val + sigma_both, color='k', alpha = 0.2)

    plt.axvline(lambda_from_nu(nu_min_US).value, color=color_ULTRASAT)
    plt.axvline(lambda_from_nu(nu_max_US).value, color=color_ULTRASAT)

    plt.axvline(lambda_from_nu(nu_max_gFUV).value, color=color_FUV)
    plt.axvline(lambda_from_nu(nu_min_gNUV).value, color=color_FUV)

    plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_max_US).value,color=color_ULTRASAT,alpha=0.1)
    plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_min_US).value,wave[-2]*(1+zall[-1]),color=color_ULTRASAT,alpha=0.1)

    plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_min_gNUV).value,wave[-2]*(1+zall[-1]),color=color_FUV,alpha=0.2)
    plt.fill_betweenx([3e24,1e28],lambda_from_nu(nu_max_gFUV).value,color=color_FUV,alpha=0.2)

    plt.ylabel(r'$\epsilon_\nu(z=0)$',fontsize=fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(loc=1)
    plt.xlim(912,wave[-2]*(1+zall[-1]))
    plt.ylim(3e24,1e28)

    plt.title(r'$\rm (GALEX + ULTRASAT)\times DESI)$')
    plt.xlabel(r'$\lambda_{\rm obs}\,[{\rm A}]$',fontsize=fontsize)

    if use_pars_fid == pars_original_c18:
        filename = 'results/PLOTS/EBL/err_rest_emissivity_combined.png'
    elif use_pars_fid == pars_all:
        filename = 'results/PLOTS/EBL/err_rest_emissivity_allpars_combined.png'

    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.show()


    return 




###

def err_bdJ(use_pars_fid = pars_all, prior=False):

    z = np.linspace(zmin_gal,zmax_gal,200)
    sigmas_emissivity, wave, signals = plot_err_emissivity(z, use_pars_fid, prior, False)
    
    detector = ['GALEX']#,'ULTRASAT','GALEX_ULTRASAT']

    #sigmas_bias  = [[],[],[]] 
    #for i in range(len(wave)):     
    #    temp = plot_err_bias(wave[i], use_pars_fid ,group_vox=True,run=False,plot_flag = True)
    #    sigmas_bias[0].append(temp[0])
    #    sigmas_bias[1].append(temp[1])
    #    sigmas_bias[2].append(temp[2])

    for i in detector:

        #bias = bJ_z(z, i, run = False, vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias=False, filename = '')

        det_i = 'GALEX_NUV' if i == 'GALEX' else 'ULTRASAT'
        label_i = r'$\rm GALEX \times SDSS$' if i == 'GALEX' else r'$\rm ULTRASAT \times DESI$' if i == 'ULTRASAT' else r'$\rm (G\times SDSS)+(U\times DESI)$'
        dJ = dJdz(z, i,run = False,\
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False, filename = './results/EBL/dJdz_' + det_i + '_reduced.dat')

        #bdJ = bias * dJ

        sigma_eps = sigmas_emissivity[detector.index(i)]
        #sigma_bias = sigmas_bias[detector.index(i)]

        sigma_dJ = np.zeros(len(z))
        for j in range(len(z)):
            intg = np.zeros(len(wave))
            for w in range(len(wave)):
                intg[w] = (Response(wave[w]*(1+z[j])*u.AA, i) * sigma_eps[w] * np.exp(-tau_Lya(wave[w]*u.AA*(1+z[j]),z[j]))/(nu_from_lambda((wave[w]*(1+z[j])*u.AA))).value)
            sigma_dJ[j] = (cu.c.to(u.km/u.s)).value / (4*np.pi*H(z[j]).value*(1+z[j])) * np.trapz(intg,wave)

        plt.plot(z,dJ)
        plt.errorbar(z,dJ,yerr=sigma_dJ,label=label_i)

    plt.xlabel(r'$z$')
    plt.ylabel(r'$b_{1500}^{z=0}dJ/dz$')

    plt.show()
    return 












def error_dJ_biasfix(group_vox,run = False):

    # !!! WE FIX THE BIAS 

    z_val = list(z_gals('SDSS')) 

    rest_wave = np.concatenate((np.linspace(91.2,121.6,50), np.linspace(121.6,300,50)))*u.nm

    use_pars_fid = ['log10_eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_0.3','EW_1',
                    #'log_fLyC_1','log_fLyC_2','bias_1500_0',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps1500'][0]),
        fiducials['eps1500'][1],
        fiducials['alpha1500'][0],
        fiducials['alpha1500'][1],
        fiducials['alpha1100'][0],
        fiducials['alpha1100'][1],
        fiducials['EW'][0],
        fiducials['EW'][1],
        #fiducials['fescape'][0],
        #fiducials['fescape'][1],
        #fiducials['bias'][0],
        fiducials['bias'][1],
        fiducials['bias'][2]
        ]
    
#    if detector == 'ULTRASAT':
    F_U = fisher_matrix(use_pars_fid,'ULTRASAT','SPHEREx',group_vox,run)
#    else:
    F_N = fisher_matrix(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F = fisher_matrix(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_all = F_N + F_F
    F_UF = F_U + F_F + F_N 

    for j in range(len(F_UF)):
        if use_pars_fid[j] == 'gamma_1500':
            F_UF[j,j] += 1/.3**2
            F_all[j,j] += 1/.3**2
        if use_pars_fid[j] == 'C_alpha_1500':
            F_UF[j,j] += 1/1.5**2
            F_all[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_1100':
            F_UF[j,j] += 1/1.5**2
            F_all[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'EW_0.3':
            F_UF[j,j] += 1/1**2
            F_all[j,j] += 1/1**2
        if use_pars_fid[j] == 'EW_1':
            F_UF[j,j] += 1/1**2
            F_all[j,j] += 1/1**2

    inv_F_all = np.linalg.inv(F_all)
    inv_F_UF_all = np.linalg.inv(F_UF)
 
    # remove the 2 bias parameters
    inv_F = np.zeros((len(inv_F_all)-2,len(inv_F_all)-2))
    inv_F_UF = np.zeros((len(inv_F_UF_all)-2,len(inv_F_UF_all)-2))

    for i in range(len(inv_F)):
        for j in range(len(inv_F)):
            inv_F[i,j] = inv_F_all[i,j]

    for i in range(len(inv_F_UF)):
        for j in range(len(inv_F_UF)):
            inv_F_UF[i,j] = inv_F_UF_all[i,j]

#    if detector == 'GALEX_NUV':
#        nu_min = nu_min_gNUV
#        nu_max = nu_max_gNUV
#
#    elif detector == 'GALEX_FUV':
#       nu_min = nu_min_gFUV
    nu_max = nu_max_gFUV

#   elif detector == 'ULTRASAT':
    nu_min = nu_min_US
#       nu_max = nu_max_US
#    else:
#        print('Check detector!')
#        return 

    nu_obs_val = np.linspace(nu_min.value,nu_max.value)
    sigma_dJ_FN = np.zeros(len(z_val))
    sigma_dJ_FU = np.zeros(len(z_val))
    sigma_dJ_U = np.zeros(len(z_val))
    dJ_F = np.zeros(len(z_val))
    dJ_U = np.zeros(len(z_val))
    detector = 'GALEX_NUV'
    for use_z in z_val:

        use_signal = lambda w: signal(w,use_z,detector,False,False,False,False,False,False)

        depsnu_deps1500_0 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) / pow(10,use_fid_vals[use_pars_fid.index('log10_eps_1500_0')])

        depsnu_dgamma1500 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * np.log10(1+use_z)

        coef_alpha1500 = lambda w: np.log(nu_from_lambda(w) / nu_from_lambda(1500*u.nm)) if w.value > 121.6 else np.log(nu_from_lambda(121.6*u.nm) / nu_from_lambda(1500*u.nm))

        depsnu_dalpha1500_0 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * coef_alpha1500(w)

        depsnu_dC1500 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * coef_alpha1500(w) * np.log10(use_z)

        depsnu_dalpha1100_0 = lambda w: (signal(w,use_z,detector,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm))).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha1100(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))).value if w.value <= 121.6 else 0.


        depsnu_dC1100 = lambda w: (signal(w,use_z,detector,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha1100(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 121.6 else 0.


        dEW_dEW03 = 1 - np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        dEW_dEW1 = np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        depsnu_dEW03 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW03).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 

        depsnu_dEW1 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW1).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 
        
        Jacobian = lambda w: np.asarray((depsnu_deps1500_0(w).value, depsnu_dgamma1500(w).value, depsnu_dalpha1500_0(w).value,depsnu_dC1500(w).value,depsnu_dalpha1100_0(w),depsnu_dC1100(w),depsnu_dEW03(w),depsnu_dEW1(w)))

        signal_val = np.zeros(len(rest_wave))
        sigma_eps_nu = np.zeros(len(rest_wave))
        sigma_eps_nu_UF = np.zeros(len(rest_wave))
        for ww in tqdm(range(len(rest_wave))):
            
            signal_val[ww] = use_signal(rest_wave[ww]).value

            sigma_eps_nu[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F,Jacobian(rest_wave[ww])]))
            sigma_eps_nu_UF[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F_UF,Jacobian(rest_wave[ww])]))
        
        sigma_eps_nu_use= lambda w: np.sqrt(np.linalg.multi_dot([Jacobian(w),inv_F,Jacobian(w).T])) * (use_signal(w)).unit
        sigma_eps_nu_UF_use= lambda w: np.sqrt(np.linalg.multi_dot([Jacobian(w),inv_F_UF,Jacobian(w).T])) * (use_signal(w)).unit


        dlambda_dnu = lambda nu_obs: (cu.c.to(u.km/u.s) / nu_obs).to(u.nm)

        intg_FN = np.zeros(len(nu_obs_val))
        intg_FU = np.zeros(len(nu_obs_val))
        intg_U = np.zeros(len(nu_obs_val))
        for n in range(len(nu_obs_val)):
            intg_FN[n] =  Response(lambda_from_nu(nu_obs_val[n]*u.Hz), 'GALEX_FUV') * sigma_eps_nu_use(lambda_from_nu(nu_obs_val[n]*u.Hz*(1+use_z))).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs_val[n]*u.Hz),use_z))/nu_obs_val[n]
            intg_FU[n] =  Response(lambda_from_nu(nu_obs_val[n]*u.Hz), 'GALEX_FUV') * sigma_eps_nu_UF_use(lambda_from_nu(nu_obs_val[n]*u.Hz*(1+use_z))).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs_val[n]*u.Hz),use_z))/nu_obs_val[n]
            intg_U[n] =  Response(lambda_from_nu(nu_obs_val[n]*u.Hz), 'ULTRASAT') * sigma_eps_nu_UF_use(lambda_from_nu(nu_obs_val[n]*u.Hz*(1+use_z))).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs_val[n]*u.Hz),use_z))/nu_obs_val[n]

        unit = (dlambda_dnu(1*u.Hz)* (u.nm)**-1 * sigma_eps_nu_use(lambda_from_nu(1*u.Hz)) /(1*u.Hz)).unit

        sigma_dJ_FN[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_FN,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + 'GALEX_FUV' + '.dat')

        sigma_dJ_FU[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_FU,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + 'GALEX_FUV' + '.dat')

        dJ_F[z_val.index(use_z)] = dJdz(use_z,'GALEX_FUV',False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + 'GALEX_FUV' + '.dat') * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + 'GALEX_FUV' + '.dat')
                
        sigma_dJ_U[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_U,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + 'ULTRASAT' + '.dat')

        dJ_U[z_val.index(use_z)] = dJdz(use_z,'ULTRASAT',False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + 'ULTRASAT' + '.dat') * bJ_z(use_z,'ULTRASAT',run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + 'ULTRASAT' + '.dat')
        
    #detector_name = r'$GALEX\, NUV$' if detector == 'GALEX_NUV' else r'$GALEX\, FUV$' if detector == 'GALEX_FUV' else r'$ULTRASAT$' if detector == 'ULTRASAT' else -1
    plt.plot(z_val,dJ_F,label=r'$\rm GALEX\,FUV\,signal$',color=color_FUV)
    plt.fill_between(z_val,dJ_F-sigma_dJ_FN,dJ_F-sigma_dJ_FU,color=color_FUV,alpha=0.3)
    plt.fill_between(z_val,dJ_F+sigma_dJ_FN,dJ_F+sigma_dJ_FU,color=color_FUV,alpha=0.3,label=r'$\rm NUV+FUV$')
    plt.fill_between(z_val,dJ_F-sigma_dJ_FU,dJ_F+sigma_dJ_FU,color=color_ULTRASAT,alpha=0.3,label=r'$\rm GALEX+ULTRASAT$')#+{\it u}$')
    plt.plot(z_val,dJ_U,label=r'$\rm ULTRASAT\,signal$',color=color_ULTRASAT)
    plt.fill_between(z_val,dJ_U-sigma_dJ_FU,dJ_U+sigma_dJ_FU,color=color_ULTRASAT,alpha=0.3)
    plt.xlim(z_val[0],z_val[-1])
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$dJ/dz\,[{\rm Jy/sr}]$',fontsize=fontsize)
    plt.legend(loc=4)

    plt.tight_layout()

    #if scale_physical_max == 300*u.Mpc: 
    #    filename = 'PLOTS/ULTRASAT/workshop/err_bJdJ_thetamax.png'
    #else:
    #    filename = 'PLOTS/ULTRASAT/workshop/err_bJdJ.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return sigma_dJ_FU

# !!!!
# !!!! TO FIX 
def error_dJ_ULTRASAT(group_vox,run = False):

    detector_all = ['GALEX_NUV','GALEX_FUV','ULTRASAT']
    z_val = list(z_gals('SDSS'))

    rest_wave = np.concatenate((np.linspace(91.2,121.6,200), np.linspace(121.6,350,200)))*u.nm

    use_pars_fid = ['log10_eps_1500_0','gamma_1500','alpha_1500_0','C_alpha_1500','alpha_1100_0','C_alpha_1100','EW_0.3','EW_1',#'log_fLyC_1','log_fLyC_2',
                    'bias_1500_0',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps1500'][0]),
        fiducials['eps1500'][1],
        fiducials['alpha1500'][0],
        fiducials['alpha1500'][1],
        fiducials['alpha1100'][0],
        fiducials['alpha1100'][1],
        fiducials['EW'][0],
        fiducials['EW'][1],
        #fiducials['fescape'][0],
        #fiducials['fescape'][1],
        fiducials['bias'][0],
        fiducials['bias'][1],
        fiducials['bias'][2]
        ]
    
    F_N = fisher_matrix(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F = fisher_matrix(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_U = fisher_matrix(use_pars_fid,'ULTRASAT','SPHEREx',group_vox,run)
    F_all = F_N + F_F + F_U
    F_GAL = F_N + F_F

    #for i in range(len(F_all)):
    #    if use_pars_fid[i] == 'gamma_1500':
    #       F_all[i,i] += 1/.3**2
    inv_F_all = np.linalg.inv(F_all)
    inv_F_GAL_all = np.linalg.inv(F_GAL)
    inv_F_u_all = np.linalg.inv(F_U)
 
    # remove the 2 bias parameters
    inv_F = np.zeros((len(inv_F_all)-3,len(inv_F_all)-3))
    inv_F_GAL = np.zeros((len(inv_F_all)-3,len(inv_F_all)-3))
    inv_F_u = np.zeros((len(inv_F_all)-3,len(inv_F_all)-3))

    plt.figure(figsize=figsize)

    subplot = [131,132,133]
    for i in range(len(inv_F)):
        for j in range(len(inv_F)):
            inv_F[i,j] = inv_F_all[i,j]
            inv_F_GAL[i,j] = inv_F_GAL_all[i,j]
            inv_F_u[i,j] = inv_F_u_all[i,j]

    for detector in detector_all:
        if detector == 'GALEX_NUV':
            nu_min = nu_min_gNUV
            nu_max = nu_max_gNUV

        elif detector == 'GALEX_FUV':
            nu_min = nu_min_gFUV
            nu_max = nu_max_gFUV

        elif detector == 'ULTRASAT':
            nu_min = nu_min_US
            nu_max = nu_max_US
        else:
            print('Check detector!')
            return 

        eps_large = [3 *10**fid_vals[0],fid_vals[1]]
        bias_large = [3*fid_vals[-3],fid_vals[-2],fid_vals[-1]]

        nu_obs_val = np.linspace(nu_min.value,nu_max.value)
        sigma_dJ = np.zeros(len(z_val))
        sigma_dJ_det = np.zeros(len(z_val))
        dJ = np.zeros(len(z_val))
        dJ_eps = np.zeros(len(z_val))
        dJ_bias = np.zeros(len(z_val))
        for use_z in z_val:

            use_signal = lambda w: signal(w,use_z,detector,False,False,False,False,False,False)

            depsnu_deps1500_0 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) / pow(10,use_fid_vals[use_pars_fid.index('log10_eps_1500_0')])

            depsnu_dgamma1500 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * np.log10(1+use_z)

            coef_alpha1500 = lambda w: np.log(nu_from_lambda(w) / nu_from_lambda(1500*u.nm)) if w.value > 121.6 else np.log(nu_from_lambda(121.6*u.nm) / nu_from_lambda(1500*u.nm))

            depsnu_dalpha1500_0 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * coef_alpha1500(w)

            depsnu_dC1500 = lambda w: signal(w,use_z,detector,False,False,False,False,False,False) * coef_alpha1500(w) * np.log10(use_z)

            depsnu_dalpha1100_0 = lambda w: (signal(w,use_z,detector,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm))).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha1100(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))).value if w.value <= 121.6 else 0.


            depsnu_dC1100 = lambda w: (signal(w,use_z,detector,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha1100(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 121.6 else 0.


            dEW_dEW03 = 1 - np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
            dEW_dEW1 = np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
            depsnu_dEW03 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW03).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 

            depsnu_dEW1 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW1).value  if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 
            
            Jacobian = lambda w: np.asarray((depsnu_deps1500_0(w).value, depsnu_dgamma1500(w).value, depsnu_dalpha1500_0(w).value,depsnu_dC1500(w).value,depsnu_dalpha1100_0(w),depsnu_dC1100(w),depsnu_dEW03(w),depsnu_dEW1(w)))

            signal_val = np.zeros(len(rest_wave))
            sigma_eps_nu = np.zeros(len(rest_wave))
            sigma_eps_nu_det = np.zeros(len(rest_wave))
            for ww in tqdm(range(len(rest_wave))):
                
                signal_val[ww] = use_signal(rest_wave[ww]).value

                sigma_eps_nu[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F,Jacobian(rest_wave[ww])]))
                if detector == 'ULTRASAT':
                    sigma_eps_nu[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F_u,Jacobian(rest_wave[ww])]))

                else:
                    sigma_eps_nu[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F_GAL,Jacobian(rest_wave[ww])]))

            sigma_eps_nu_use= lambda w: np.sqrt(np.linalg.multi_dot([Jacobian(w),inv_F,Jacobian(w).T])) * (use_signal(w)).unit

            if detector == 'ULTRASAT':
                sigma_eps_nu_use_det= lambda w: np.sqrt(np.linalg.multi_dot([Jacobian(w),inv_F_u,Jacobian(w).T])) * (use_signal(w)).unit
            else:
                sigma_eps_nu_use_det= lambda w: np.sqrt(np.linalg.multi_dot([Jacobian(w),inv_F_GAL,Jacobian(w).T])) * (use_signal(w)).unit

            dlambda_dnu = lambda nu_obs: (cu.c.to(u.km/u.s) / nu_obs).to(u.nm)

            intg = np.zeros(len(nu_obs_val))
            intg_det = np.zeros(len(nu_obs_val))
            for n in range(len(nu_obs_val)):
                intg[n] = dlambda_dnu(nu_obs_val[n]*u.Hz).value * Response(lambda_from_nu(nu_obs_val[n]*u.Hz), detector) * sigma_eps_nu_use(lambda_from_nu(nu_obs_val[n]*u.Hz*(1+use_z))).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs_val[n]*u.Hz),use_z))/nu_obs_val[n]

                intg_det[n] = dlambda_dnu(nu_obs_val[n]*u.Hz).value * Response(lambda_from_nu(nu_obs_val[n]*u.Hz), detector) * sigma_eps_nu_use_det(lambda_from_nu(nu_obs_val[n]*u.Hz*(1+use_z))).value * np.exp(-tau_Lya(lambda_from_nu(nu_obs_val[n]*u.Hz),use_z))/nu_obs_val[n]

            unit = (dlambda_dnu(1*u.Hz)* (u.nm)**-1 * sigma_eps_nu_use(lambda_from_nu(1*u.Hz)) /(1*u.Hz)).unit

            if detector == 'ULTRASAT':
            
                sigma_dJ[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')

                sigma_dJ_det[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_det,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')

                dJ[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')
            
                dJ_eps[z_val.index(use_z)] = dJdz(use_z,detector,True,eps_large,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '_eps1500_0_large.dat') * bJ_z(use_z,detector,run=True,vals_eps1500=eps_large,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '_eps1500_0_large.dat')
            
                dJ_bias[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = bias_large,filename='results/EBL/bJ_' + detector + '_b1500_0_large.dat')
            
            else:
                sigma_dJ[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')

                sigma_dJ_det[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_det,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')

                dJ[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=False,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '.dat')

                dJ_eps[z_val.index(use_z)] = dJdz(use_z,detector,True,eps_large,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '_eps1500_0_large.dat') * bJ_z(use_z,detector,run=True,vals_eps1500=eps_large,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = False,filename='results/EBL/bJ_' + detector + '_eps1500_0_large.dat')
            
                dJ_bias[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=True,vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,val_bias = bias_large,filename='results/EBL/bJ_' + detector + '_b1500_0_large.dat')

        plt.subplot(subplot[detector_all.index(detector)])

        detector_name = r'$GALEX\, NUV$' if detector == 'GALEX_NUV' else r'$GALEX\, FUV$' if detector == 'GALEX_FUV' else r'$ULTRASAT$' if detector == 'ULTRASAT' else -1
        plt.plot(z_val,dJ,'k',linestyle='-',label=detector_name)
        plt.plot(z_val,dJ_eps,'k',linestyle='--',label=r'$3\epsilon_{1500}^0$')
        plt.plot(z_val,dJ_bias,'k',linestyle='-.',label=r'$3b_{1500}^0$')
        #plt.fill_between(z_val,dJ-sigma_dJ_det,dJ+sigma_dJ_det,color='r',alpha=0.1,label=detector_name)
        plt.fill_between(z_val,dJ-sigma_dJ,dJ+sigma_dJ,color=color_ULTRASAT,alpha=0.1,label=r'$Combined$')
        plt.xlim(z_val[0],z_val[-1])
        #plt.ylim(-10,400)
        #plt.ylim(-10,200)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$dJ/dz\,[{\rm Jy/sr}]$')
        plt.legend(loc=3,fontsize=20,ncol=2)

    plt.tight_layout()
    if scale_physical_max == 300*u.Mpc: 
        filename = 'results/PLOTS/EBL/forecast_dJ_all_beps_thetamax.png'
    else:
        filename = 'results/PLOTS/EBL/forecast_dJ_all_beps.png'
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return sigma_dJ
