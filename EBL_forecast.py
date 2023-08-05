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
            px_size = 5.45*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 5.45*u.arcsec
    
        mag_limit = 23

    elif detector == 'CASTOR_UV' or 'CASTOR_U':        
        if not group_vox:
            px_size = 0.1*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 0.7*u.arcsec

        mag_limit = 27.4
        
    elif detector == 'CASTOR_G':
        if not group_vox:
            px_size = 0.1*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 0.7*u.arcsec

        mag_limit = 27
        
    else: 
        print('Detector not recognized!')
        return -1 

    fnu = 10**((mag_limit - 8.90)/(-2.5))*u.Jy
    fnu_px = fnu#/(((px_size)**2).to(u.steradian))

    noise_per_voxel = fnu_px/5

    return noise_per_voxel



def sigma_wz(z, detector, gal_survey, group_vox):

    delta_zc = 1e-3

    if detector == 'GALEX_NUV' or detector == 'GALEX_FUV':
        if not group_vox:
            px_size = 5*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 50.*u.arcsec
        sky_LIM = (5500*u.deg**2).to(u.steradian)

    elif detector == 'ULTRASAT':
        if not group_vox:
            px_size = 5.45*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 5.45*u.arcsec

        sky_LIM = (40000*u.deg**2).to(u.steradian)

    elif detector == 'CASTOR_UV' or 'CASTOR_U' or 'CASTOR_G':
        if not group_vox:
            px_size = 0.1*u.arcsec
        else:
            if type(group_vox) == float or type(group_vox) == int or type(group_vox) == np.float64:
                px_size = group_vox*u.arcsec
            else:
                px_size = 0.7*u.arcsec
        sky_LIM = (7700*u.deg**2).to(u.steradian)

    if scale_physical_max == 300*u.Mpc:
        if detector == 'ULTRASAT': 
            use_thetamax = theta_max(cosmo.angular_diameter_distance(z))
        else:
            use_thetamax = np.arctan((5*u.Mpc).value / cosmo.angular_diameter_distance(z))
    else:
        use_thetamax = theta_max(cosmo.angular_diameter_distance(z))

    Nvox = sky_LIM / (px_size**2).to(u.steradian)
    Nvox_std = 1 / (px_size**2).to(u.steradian) #Nvox / sky_LIM
    Nvox_inthetamax =  np.pi* use_thetamax**2 /  (px_size**2).to(u.steradian)  #Nvox_std * np.pi* use_thetamax**2

    Ngal_ref = dNgdz(z,gal_survey) * delta_zi(gal_survey)
    
    #noise = delta_zi(gal_survey) / delta_zc * use_thetamax * np.sqrt(np.pi / (Ngal_ref*Nvox/sky_LIM)) * np.sqrt(Jnu_monopole(detector)**2 + sigma_N(detector, group_vox).value**2)  

    noise = delta_zi(gal_survey) / delta_zc * use_thetamax * np.sqrt(np.pi) * np.sqrt(Jnu_monopole(detector)**2 + sigma_N(detector, group_vox).value**2) / np.sqrt((Ngal_ref*Nvox_inthetamax)) 

    return noise



def plot_noise_grouped():

    groups = np.concatenate((np.linspace(0.001,1,50),np.linspace(1,50,50)))
    sigma_GAL = np.zeros(len(groups))
    sigma_ULT = np.zeros(len(groups))
    sigma_CAS = np.zeros(len(groups))
    
    for i in range(len(groups)):
        sigma_GAL[i] = sigma_wz(0.5,'GALEX_FUV','SDSS',groups[i]).value
        sigma_ULT[i] = sigma_wz(1,'ULTRASAT','SPHEREx',groups[i]).value
        sigma_CAS[i] = sigma_wz(1,'CASTOR_G','SDSS',groups[i]).value

    plt.plot(groups,sigma_GAL,color_FUV,label=r'${\rm FUV,} z = 0.5$')
    plt.plot([5],[sigma_wz(0.5,'GALEX_FUV','SDSS',False).value],color_FUV,marker='D')
    plt.plot([50],[sigma_wz(0.5,'GALEX_FUV','SDSS',True).value],color_FUV,marker='x',markersize=15)
    plt.plot(groups,sigma_ULT,color_ULTRASAT,label=r'${\rm ULTRASAT,} z = 1$')
    plt.plot([5.45],[sigma_wz(1,'ULTRASAT','SPHEREx',False).value],color_ULTRASAT,marker='D')
    plt.plot([5.45],[sigma_wz(1,'ULTRASAT','SPHEREx',True).value],color_ULTRASAT,marker='x',markersize=15)
    plt.plot(groups,sigma_CAS,color_CASTOR,label=r'$g, z = 1$')
    plt.plot([0.1],[sigma_wz(1,'CASTOR_G','SDSS',False).value],color_CASTOR,marker='D')
    plt.plot([0.7],[sigma_wz(1,'CASTOR_G','SDSS',True).value],color_CASTOR,marker='x',markersize=15)

    plt.yscale('log')
    plt.xlabel(r'$\rm Voxel size\,[arcsec]$',fontsize=fontsize)
    plt.ylabel(r'$\sigma_{\omega_{Jg}}$',fontsize=fontsize)
    plt.legend(loc=1)

    plt.tight_layout()

    if scale_physical_max == 300*u.Mpc: 
        filename = 'PLOTS/ULTRASAT/EBL/noise_voxgroup_thetamax.png'
    else:
        filename = 'PLOTS/ULTRASAT/EBL/noise_voxgroup.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 


def plot_signal_and_noise(detector,gal_survey):

    use_z = z_small if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else z_small_castor if detector == 'CASTOR_UV' or detector == 'CASTOR_U' or detector == 'CASTOR_G' else -1

    signal = np.zeros(len(use_z))

    if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
        filename = 'results/map_EBL/wJg_' + detector +',' + gal_survey + '_thetamax.dat'
    else:
        filename = 'results/map_EBL/wJg_' + detector +',' + gal_survey + '.dat'

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

    color = color_ULTRASAT if detector == 'ULTRASAT' else color_NUV  if detector == 'GALEX_NUV' else color_FUV if detector == 'GALEX_FUV' else color_CASTOR

    plt.plot(use_z,signal,label=r'$\rm %s\times$'%detector + r'$\rm %s$'%gal_survey, color=color)

    #label_grouped = r'$3\times 3 \,{\rm voxels}$' if detector == 'ULTRASAT' else r'$\rm Voxels = 50\,arcsec$'  

    #plt.plot(use_z,nsingle,'k:', label=r'$\rm Noise\,per\,voxel$')
    
    #label_grouped = r'${\rm Noise:}\, 3\times 3 \rm voxels$' if detector == 'ULTRASAT' else r'$\rm Voxels = 50\,arcsec$'  
    plt.plot(use_z,ngrouped,'k-',label=r'$\rm Noise$')#label_grouped)
    #if detector == 'ULTRASAT':
    #    plt.plot(use_z,ngrouped50,'k--',label= r'$\rm Noise:\, 50\,arcsec\, voxels$')

    plt.xlabel(r'$z$')
    plt.ylabel(r'$|\sigma_{w_{J_\nu{\rm g}}(z)}|$')
    plt.yscale('log')
    plt.legend(loc=1)

    plt.tight_layout()

    if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
        filefig = 'PLOTS/ULTRASAT/EBL/noise_' + detector + ',' + gal_survey + '_thetamax.png'
    else:
        filefig = 'PLOTS/ULTRASAT/EBL/noise_' + detector + ',' + gal_survey + '.png'
    
    #plt.savefig(filefig,bbox_inches='tight')
    plt.show()

    return 

def plot_noise():

    sn = np.zeros(len(z_small))
    sf = np.zeros(len(z_small))
    sU = np.zeros(len(z_small))
    sCUV = np.zeros(len(z_small_castor))
    sCU = np.zeros(len(z_small_castor))
    sCG = np.zeros(len(z_small_castor))
    nG = np.zeros(len(z_small))
    nU = np.zeros(len(z_small))
    nCUVU = np.zeros(len(z_small_castor))
    nCG = np.zeros(len(z_small_castor))


    if scale_physical_max == 300*u.Mpc:
        filename_ULT = 'results/map_EBL/wJg_ULTRASAT,SPHEREx_thetamax.dat'
    else:
        filename_ULT = 'results/map_EBL/wJg_ULTRASAT,SPHEREx.dat'
    
    for i in tqdm(range(len(z_small))):
        sn[i] = abs(wJgz(z_small[i],'GALEX_NUV','SDSS',False,filename='results/map_EBL/wJg_GALEX_NUV,SDSS.dat'))
        sf[i] = abs(wJgz(z_small[i],'GALEX_FUV','SDSS',False,filename='results/map_EBL/wJg_GALEX_FUV,SDSS.dat'))
        sU[i] = abs(wJgz(z_small[i],'ULTRASAT','SPHEREx',False,filename=filename_ULT))
#
        nG[i] = sigma_wz(z_small[i],'GALEX_NUV','SDSS',True).value
        nU[i] = sigma_wz(z_small[i],'ULTRASAT','SPHEREx',True).value

    for i in tqdm(range(len(z_small_castor))):
        sCUV[i] = abs(wJgz(z_small_castor[i],'CASTOR_UV','SDSS',False,filename='results/map_EBL/wJg_CASTOR_UV,SDSS.dat'))
        sCU[i] = abs(wJgz(z_small_castor[i],'CASTOR_U','SDSS',False,filename='results/map_EBL/wJg_CASTOR_U,SDSS.dat'))
        sCG[i] = abs(wJgz(z_small_castor[i],'CASTOR_G','SDSS',False,filename='results/map_EBL/wJg_CASTOR_G,SDSS.dat'))

        nCUVU[i] = sigma_wz(z_small_castor[i],'CASTOR_UV','SDSS',True).value
        nCG[i] = sigma_wz(z_small_castor[i],'CASTOR_G','SDSS',True).value

        
    #plt.plot(z_small,sU,color_ULTRASAT,label=r'$ULTRASAT$')
    #plt.plot(z_small,sn,label=r'$GALEX\,\,NUV$')
    #plt.plot(z_small,sf,label=r'$GALEX\,\,FUV$')
    #plt.plot(z_small_castor,sCUV,label=r'$CASTOR\,\,uv$')
    #plt.plot(z_small_castor,sCU,label=r'$CASTOR\,\,u$')
    #plt.plot(z_small_castor,sCG,label=r'$CASTOR\,\,g$')

    plt.plot(z_small_castor,np.ones(len(z_small_castor)),'k',linewidth=.9)

    plt.plot(z_small,nU/sU,color_ULTRASAT,label=r'$\rm ULTRASAT\times SPHEREx$')
    plt.plot(z_small,nG/sf,color_FUV,label=r'$\rm \{NUV, FUV\}\times SDSS$')
    plt.plot(z_small_castor,nCUVU/sCUV,color_CASTOR,linestyle=':',label=r'$\rm \{{\it uv,u}\}\times SDSS$')
    plt.plot(z_small_castor,nCG/sCG,color_CASTOR,linestyle='-.',label=r'$\rm {\it g}\times SDSS$')
#
    plt.yscale('log')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\sigma_{w_{J_\nu {\rm g}}(z)}/w_{J_\nu {\rm g}}$')
    plt.legend(loc=1,ncol=1,fontsize=fontsize)

#    plt.ylim(1e-7,5e2)
    plt.xlim(z_small_castor[0],z_small_castor[-1])
    plt.tight_layout()


    if scale_physical_max == 300*u.Mpc: 
        filename = 'PLOTS/ULTRASAT/EBL/wJg_noise_thetamax.png'
    else:
        filename = 'PLOTS/ULTRASAT/EBL/wJg_noise.png'
    
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return 



pars_fid = ['log10_eps_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110','EW_0.3','EW_1','log_fLyC_1','log_fLyC_2',
            'bias_150_0','gamma_bv','gamma_bz']
#
fid_vals = [np.log10(fiducials['eps150'][0]),
      fiducials['eps150'][1],
      fiducials['alpha150'][0],
      fiducials['alpha150'][1],
      fiducials['alpha110'][0],
      fiducials['alpha110'][1],
      fiducials['EW'][0],
      fiducials['EW'][1],
      fiducials['fescape'][0],
      fiducials['fescape'][1],
      fiducials['bias'][0],
      fiducials['bias'][1],
      fiducials['bias'][2]
      ]


delta_der = 0.1

def ders(parameter,detector,gal_survey):

    use_z = z_small if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else z_small_castor if detector == 'CASTOR_UV' or detector == 'CASTOR_U' or detector == 'CASTOR_G' else -1

    print('Doing der: ' + str(parameter))

    if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
        filename = 'results/map_EBL/der/d' + parameter + '_' + detector + ',' + gal_survey + '_thetamax.dat'
        filename_wm = 'results/map_EBL/wmz_thetamax.dat'
    else:
        filename = 'results/map_EBL/der/d' + parameter + '_' + detector + ',' + gal_survey + '.dat'
    
        if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT':
            filename_wm = 'results/map_EBL/wmz.dat'
        elif detector == 'CASTOR_UV' or detector == 'CASTOR_U' or detector == 'CASTOR_G':
            filename_wm = 'results/map_EBL/wmz_castor.dat'
        else:
            filename_wm = -1

    filename_dJ = 'results/map_EBL/der/dJ_d' + parameter + '_' + detector + ',' + gal_survey + '.dat'
    filename_bJ = 'results/map_EBL/der/bJ_d' + parameter + '_' + detector + ',' + gal_survey + '.dat'

    filename_dJ_fid = 'results/map_EBL/dJdz_' + detector + '.dat'
    filename_bJ_fid = 'results/map_EBL/bJ_' + detector + '.dat' 

    if os.path.exists(filename):
        z, der_wz = np.genfromtxt(filename)
        if list(z) != list(use_z):
            print('Redshift bins in file are different from z_arr!')

    else:
        use_pars_up = []
        for i in range(len(pars_fid)):
            
            if pars_fid[i] == parameter:
                if pars_fid[i] == 'log10_eps_150_0' or pars_fid[i] == 'log_fLyC_1' or pars_fid[i] == 'log_fLyC_2':
                    par_up = pow(10,fid_vals[i]*(1+delta_der))
                    step = pow(10,fid_vals[i]*(1+delta_der)) - pow(10,fid_vals[i])

                else:
                    par_up = fid_vals[i]*(1+delta_der)
                    step = fid_vals[i]*delta_der
            else:
                if pars_fid[i] == 'log10_eps_150_0' or pars_fid[i] == 'log_fLyC_1' or pars_fid[i] == 'log_fLyC_2':
                    par_up = pow(10,fid_vals[i])
                else:
                    par_up = fid_vals[i]

            use_pars_up.append(par_up)

        with Pool(6) as pool:

            dJ_up_f = partial(dJdz,detector=detector,run=True,\
                    vals_eps150=[use_pars_up[pars_fid.index('log10_eps_150_0')],use_pars_up[pars_fid.index('gamma_150')]],\
                    vals_alpha150=[use_pars_up[pars_fid.index('alpha_150_0')],use_pars_up[pars_fid.index('C_alpha_150')]],\
                    vals_alpha110=[use_pars_up[pars_fid.index('alpha_110_0')],use_pars_up[pars_fid.index('C_alpha_110')]],\
                    val_EW=[use_pars_up[pars_fid.index('EW_0.3')],use_pars_up[pars_fid.index('EW_1')]],\
                    val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],\
                    val_alpha90=False,\
                    filename=filename_dJ)
            
            dJ_up = pool.map(dJ_up_f, use_z)

            bJ_up_f = partial(bJ_z,detector=detector,run=True,\
                    vals_eps150=[use_pars_up[pars_fid.index('log10_eps_150_0')],use_pars_up[pars_fid.index('gamma_150')]],\
                    vals_alpha150=[use_pars_up[pars_fid.index('alpha_150_0')],use_pars_up[pars_fid.index('C_alpha_150')]],\
                    vals_alpha110=[use_pars_up[pars_fid.index('alpha_110_0')],use_pars_up[pars_fid.index('C_alpha_110')]],\
                    val_EW=[use_pars_up[pars_fid.index('EW_0.3')],use_pars_up[pars_fid.index('EW_1')]],\
                    val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],\
                    val_alpha90=False,\
                    val_bias=[use_pars_up[pars_fid.index('bias_150_0')],use_pars_up[pars_fid.index('gamma_bv')],use_pars_up[pars_fid.index('gamma_bz')]],\
                    filename=filename_bJ)
            
            bJ_up = pool.map(bJ_up_f, use_z)

        np.savetxt(filename_dJ,(use_z,dJ_up))
        np.savetxt(filename_bJ,(use_z,bJ_up))


        if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
            filename_fid = 'results/map_EBL/wJg_' + detector + ',' + gal_survey + '_thetamax.dat'
        else:
            filename_fid = 'results/map_EBL/wJg_' + detector + ',' + gal_survey + '.dat'
       
        if os.path.exists(filename_fid):
                run_fid = False 
        else:
                run_fid = True       
       
        der_wz = np.zeros(len(use_z))
        for z in tqdm(range(len(use_z))):
            wz_up = wJgz(use_z[z],detector, gal_survey,
                True,
                vals_eps150=[use_pars_up[pars_fid.index('log10_eps_150_0')],use_pars_up[pars_fid.index('gamma_150')]],
                vals_alpha150=[use_pars_up[pars_fid.index('alpha_150_0')],use_pars_up[pars_fid.index('C_alpha_150')]],
                vals_alpha110=[use_pars_up[pars_fid.index('alpha_110_0')],use_pars_up[pars_fid.index('C_alpha_110')]],
                val_EW=[use_pars_up[pars_fid.index('EW_0.3')],use_pars_up[pars_fid.index('EW_1')]],
                val_flyc=[use_pars_up[pars_fid.index('log_fLyC_1')],use_pars_up[pars_fid.index('log_fLyC_2')]],
                val_alpha90=False,
                val_bias=[use_pars_up[pars_fid.index('bias_150_0')],use_pars_up[pars_fid.index('gamma_bv')],use_pars_up[pars_fid.index('gamma_bz')]],\
                #filename = filename, 
                filename_wm = filename_wm,
                filename_dJ = filename_dJ,
                filename_bJ =filename_bJ)

            wz = wJgz(use_z[z],detector,gal_survey,
                run=run_fid,
                vals_eps150=False,
                vals_alpha150=False,
                vals_alpha110=False,
                val_EW=False,
                val_flyc=False,
                val_alpha90=False,
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


def fisher_matrix(pars,detector,gal_survey,group_vox,run = False):

    if group_vox: 
        if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
            filename = 'results/map_EBL/FISHER_' + detector + ',' + gal_survey + '_thetamax.dat'
        else:
            filename = 'results/map_EBL/FISHER_' + detector + ',' + gal_survey + '.dat'
    else:
        if scale_physical_max == 300*u.Mpc and detector == 'ULTRASAT': 
            filename = 'results/map_EBL/FISHER_' + detector + ',' + gal_survey + '_singleVOX_thetamax.dat'
        else:
            filename = 'results/map_EBL/FISHER_' + detector + ',' + gal_survey + '_singleVOX.dat'

    if not run and os.path.exists(filename):
        Fisher_start = np.genfromtxt(filename)

        Fisher = np.zeros((len(pars),len(pars)))
        for i in range(len(pars)):
            for j in range(len(pars)):
                Fisher[i,j] = Fisher_start[pars_fid.index(pars[i])][pars_fid.index(pars[j])]
        return Fisher

    Fisher = np.zeros((len(pars),len(pars)))
    
    z = z_SDSS if gal_survey == 'SDSS' else z_SPHEREx if gal_survey == 'SPHEREx' else -1

    use_z = z_small if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else z_small_castor if detector == 'CASTOR_UV' or detector == 'CASTOR_U' or detector == 'CASTOR_G' else -1

    der_pars = {}
    for p in range(len(pars)):
        temp = ders(pars[p],detector,gal_survey)
        all_zbin_der = interp1d(use_z,temp)

        vals = np.zeros(len(z))
        for i in range(len(vals)):
            try:
                vals[i] = all_zbin_der(z[i])
            except:
                if z[i] < use_z[0]:
                    vals[i] = all_zbin_der(use_z[0])
                elif z[i] > use_z[-1]:
                    vals[i] = all_zbin_der(use_z[-1])
        der_pars[pars[p]] = vals 

    sigma2 = np.zeros(len(z))
    for zv in range(len(z)):

        sigma2[zv] = sigma_wz(z[zv],detector,gal_survey,group_vox).value**2

    Fisher = np.zeros((len(pars),len(pars)))
    for p1 in range(len(pars)):
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

    dlog10eps_dP = 1
    db_dP = np.log(10) * fid_vals[pars_fid.index('bias_150_0')]

    Fisher_prime = np.zeros((len(pars),len(pars)))
    
    for i in range(len(pars)):

        if pars[i] == 'log10_epsbias_150_0':

            Fisher_prime[i,i] = F[pars_fid.index('log10_eps_150_0')][pars_fid.index('log10_eps_150_0')] * dlog10eps_dP**2 + F[pars_fid.index('log10_eps_150_0')][pars_fid.index('bias_150_0')] * dlog10eps_dP*db_dP + F[pars_fid.index('bias_150_0')][pars_fid.index('bias_150_0')] * db_dP**2

            for j in range(len(pars)):
                
                if i < j:
                    Fisher_prime[i,j] = F[pars_fid.index('log10_eps_150_0')][pars_fid.index(pars[j])] * dlog10eps_dP + F[pars_fid.index('bias_150_0')][pars_fid.index(pars[j])] * db_dP

        else:
            for j in range(len(pars)):
                if i <= j:
                    Fisher_prime[i,j] = F[pars_fid.index(pars[i])][pars_fid.index(pars[j])]

    for i in range(len(pars)):
        for j in range(len(pars)):
            if i > j:
                Fisher_prime[i,j] = Fisher_prime[j,i]

    return Fisher_prime




def sigma_epsbias(group_vox):

    print('Doing contour plots')

    pars = ['log10_eps_150_0','bias_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110','EW_0.3','EW_1','log_fLyC_1','log_fLyC_2', 'gamma_bv','gamma_bz']
#
    F_N = fisher_matrix(pars,'GALEX_NUV','SDSS',group_vox,False)
    F_F = fisher_matrix(pars,'GALEX_FUV','SDSS',group_vox,False)
    F_GALEX = F_N + F_F

    F_ULTRASAT = fisher_matrix(pars,'ULTRASAT','SPHEREx',group_vox,False)

    F_both = F_ULTRASAT + F_GALEX

    F_CASTOR_UV = fisher_matrix(pars,'CASTOR_UV','SDSS',group_vox,False)
    F_CASTOR_U = fisher_matrix(pars,'CASTOR_U','SDSS',group_vox,False)
    F_CASTOR_G = fisher_matrix(pars,'CASTOR_G','SDSS',group_vox,False)
    F_CASTOR = F_CASTOR_UV + F_CASTOR_U + F_CASTOR_G    

    for j in range(len(F_CASTOR)):
        if pars[j] == 'gamma_150':
            F_GALEX[j,j] += 1/.3**2
            F_ULTRASAT[j,j] += 1/.3**2
            F_both[j,j] += 1/.3**2
            F_CASTOR[j,j] += 1/.3**2
        if pars[j] == 'C_alpha_150':
            F_GALEX[j,j] += 1/1.5**2
            F_ULTRASAT[j,j] += 1/1.5**2
            F_both[j,j] += 1/1.5**2
            F_CASTOR[j,j] += 1/1.5**2
        if pars[j] == 'C_alpha_110':
            F_GALEX[j,j] += 1/1.5**2
            F_ULTRASAT[j,j] += 1/1.5**2
            F_both[j,j] += 1/1.5**2
            F_CASTOR[j,j] += 1/1.5**2
        if pars[j] == 'EW_0.3':
            F_GALEX[j,j] += 1/10**2
            F_ULTRASAT[j,j] += 1/10**2
            F_both[j,j] += 1/10**2
            F_CASTOR[j,j] += 1/10**2
        if pars[j] == 'EW_1':
            F_GALEX[j,j] += 1/10**2
            F_ULTRASAT[j,j] += 1/10**2
            F_both[j,j] += 1/10**2
            F_CASTOR[j,j] += 1/10**2


    inv_F_GALEX = np.linalg.inv(F_GALEX)
    inv_F_ULTRASAT = np.linalg.inv(F_ULTRASAT)
    inv_F_both = np.linalg.inv(F_both)
    inv_F_CASTOR = np.linalg.inv(F_CASTOR)

    print('EPS_150 GALEX = ' + str(np.sqrt(inv_F_GALEX[0][0])))
    print('EPS_150 ULTRASAT = ' + str(np.sqrt(inv_F_ULTRASAT[0][0])))
    print('EPS_150 both = ' + str(np.sqrt(inv_F_both[0][0])))
    print('EPS_150 CASTOR = ' + str(np.sqrt(inv_F_CASTOR[0][0])))

    print('b_150 GALEX = ' + str(np.sqrt(inv_F_GALEX[1][1])))
    print('b_150 ULTRASAT = ' + str(np.sqrt(inv_F_ULTRASAT[1][1])))
    print('b_150 both = ' + str(np.sqrt(inv_F_both[1][1])))
    print('b_150 CASTOR = ' + str(np.sqrt(inv_F_CASTOR[1][1])))

    return


def compare_surveys(pars = ['log10_epsbias_150_0', 'gamma_150', 'alpha_150_0', 'C_alpha_150', 'alpha_110_0', 'C_alpha_110', 'EW_0.3', 'EW_1', 'gamma_bv', 'gamma_bz'], detector_list = ['GALEX','ULTRASAT','GALEX_ULTRASAT'], plot_flag = True, group_vox = True, prior = False):

    print('Doing contour plots')

    F_list = []
    inv_F_list = []

    for i in detector_list:

        if i == 'GALEX':

            F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
            F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
            temp = F_N + F_F

        elif i == 'GALEX_NUV':

            F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
            temp = F_N

        elif i == 'GALEX_FUV':

            F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
            temp = F_F


        elif i == 'ULTRASAT':

            F = Fisher_change_var(pars,'ULTRASAT','SPHEREx',group_vox,False)
            temp = F

        elif i == 'GALEX_ULTRASAT':

            F_N = Fisher_change_var(pars,'GALEX_NUV','SDSS',group_vox,False)
            F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
            F = Fisher_change_var(pars,'ULTRASAT','SPHEREx',group_vox,False)
            temp = F_F + F + F_N  

        elif i == 'CASTOR_UV':

            F_CUV = Fisher_change_var(pars,'CASTOR_UV','SDSS',group_vox,False)
            temp = F_CUV 

        elif i == 'CASTOR_U':

            F_CU = Fisher_change_var(pars,'CASTOR_U','SDSS',group_vox,False)
            temp = F_CU

        elif i == 'CASTOR_G':

            F_G = Fisher_change_var(pars,'CASTOR_G','SDSS',group_vox,False)
            temp = F_G

        elif i == 'CASTOR':

            F_CUV = Fisher_change_var(pars,'CASTOR_UV','SDSS',group_vox,False)
            F_CU = Fisher_change_var(pars,'CASTOR_U','SDSS',group_vox,False)
            F_G = Fisher_change_var(pars,'CASTOR_G','SDSS',group_vox,False)
            temp = F_CUV + F_CU + F_G

        elif i == 'combined':

            F_G = Fisher_change_var(pars,'CASTOR_G','SDSS',group_vox,False)
            F_F = Fisher_change_var(pars,'GALEX_FUV','SDSS',group_vox,False)
            F = Fisher_change_var(pars,'ULTRASAT','SPHEREx',group_vox,False)
            temp = F_F + F + F_G

        if prior: 
            for j in range(len(temp)):
                if pars[j] == 'gamma_150':
                    temp[j,j] += 1/.3**2
                if pars[j] == 'C_alpha_150':
                    temp[j,j] += 1/1.5**2
                if pars[j] == 'C_alpha_110':
                    temp[j,j] += 1/1.5**2
                if pars[j] == 'EW_0.3':
                    temp[j,j] += 1/10**2
                if pars[j] == 'EW_1':
                    temp[j,j] += 1/10**2

        F_list.append(temp)
        inv_F_list.append(np.linalg.inv(temp))


    names = []
    
    print('\n')
    fiducials_pars = []
    for i in pars:
        if i == 'log10_epsbias_150_0':
            fiducials_pars.append(np.log10(10**fid_vals[pars_fid.index('log10_eps_150_0')]*fid_vals[pars_fid.index('bias_150_0')]))
        else:
            try:
                fiducials_pars.append(fid_vals[pars_fid.index(i)].value)
            except:
                fiducials_pars.append(fid_vals[pars_fid.index(i)])

        name = r'$\log_{10}(\epsilon_{150}^{z=0}b_{150}^{z=0})$' if i == 'log10_epsbias_150_0' else r'$\log_{10}\epsilon_{150}^{z=0}$' if i == 'log10_eps_150_0' else r'$\gamma_{150}$' if i == 'gamma_150' else r'$\alpha_{150}^{z=0}$' if i == 'alpha_150_0' else r'$C_{\alpha_{150}}$' if i == 'C_alpha_150' else r'$\alpha_{110}^{z=0}$' if i == 'alpha_110_0' else r'$C_{\alpha_{110}}$' if i == 'C_alpha_110' else r'$EW^{z=0.3}$' if i == 'EW_0.3' else r'$EW^{z=1}$' if i == 'EW_1' else r'$\log_{10}f_{\rm LyC}^{z=1}$' if i == 'log_fLyC_1' else r'$\log_{10}f_{\rm LyC}^{z=2}$' if i == 'log_fLyC_2' else r'$\gamma_{b_v}$' if i == 'gamma_bv' else r'$\gamma_{b_z}$' if i == 'gamma_bz' else -1

        names.append(name)

    if plot_flag:
        plot_dist_list = []
        label_names = []
        line_args_list = []
        colors_list = []
        filled_list = []
        for d in range(len(detector_list)):
            print('DETECTOR: ' + detector_list[d])
            for i in pars:
                print('--- ' + i + ': ' + str(fiducials_pars[pars.index(i)]) + ' +- ' + str(round(np.sqrt(inv_F_list[d][pars.index(i)][pars.index(i)]),4)))

            label_det = r'$\rm GALEX\times SDSS$' if detector_list[d] == 'GALEX' else r'$\rm NUV\times SDSS$' if detector_list[d] == 'GALEX_NUV' else r'$\rm FUV\times SDSS$' if detector_list[d] == 'GALEX_FUV' else  r'$\rm ULTRASAT\times SPHEREx$' if detector_list[d] == 'ULTRASAT' else r'$\rm CASTOR\times SDSS$'  if detector_list[d] == 'CASTOR' else r'$\rm uv {\times \rm SDSS}$'  if detector_list[d] == 'CASTOR_UV' else r'$\rm u {\times \rm SDSS}$'  if detector_list[d] == 'CASTOR_U' else r'$g {\times \rm SDSS}$'  if detector_list[d] == 'CASTOR_G' else r'$\rm Combined$' if detector_list[d] == 'GALEX_ULTRASAT' else -1

            label_names.append(label_det)

            color = color_FUV if detector_list[d] == 'GALEX' else color_ULTRASAT if detector_list[d] == 'ULTRASAT' else color_ULTRASAT if detector_list[d] == 'GALEX_ULTRASAT' else color_CASTOR  if detector_list[d] == 'CASTOR'else color_CASTOR  if detector_list[d] == 'CASTOR_UV' else color_CASTOR  if detector_list[d] == 'CASTOR_U' else color_CASTOR  if detector_list[d] == 'CASTOR_G' else color_FUV if detector_list[d] == 'GALEX_FUV' else color_NUV if detector_list[d] == 'GALEX_NUV' else -1
            colors_list.append(color)

            filled = True if detector_list[d] == 'GALEX_ULTRASAT' else False
            filled_list.append(filled)

            line_arg = {'ls':'-', 'color':color} 
            line_args_list.append(line_arg) 

            plot_dist_list.append(GaussianND(fiducials_pars, inv_F_list[d], names=names))
        
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

        g.triangle_plot(plot_dist_list, names, filled=filled_list,          
            colors=colors_list,line_args=line_args_list,\
            legend_labels=label_names,legend_loc='upper right',\
            markers={'x2':0})

        for i in range(len(pars)-1):
            ax = g.subplots[i+1, 0]
            ax.yaxis.set_label_coords(-0.7, 0.5)

        plt.tight_layout()
        name_plot = 'ellipse_combined' if 'GALEX_ULTRASAT' in detector_list else 'ellipse_surveys' 

        if scale_physical_max == 300*u.Mpc: 
            filename = 'PLOTS/ULTRASAT/EBL/' + name_plot + '_thetamax.png'
        else:
            filename = 'PLOTS/ULTRASAT/EBL/' + name_plot + '.png'

        plt.savefig(filename,bbox_inches='tight')
        plt.show()

    return





###################################################




def plot_parameters_errors(group_vox,run=False):

    z = z_small 

    use_pars_fid = ['log10_eps_150_0','bias_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110',#'EW_0.3','EW_1',
                    #'log_fLyC_1','log_fLyC_2',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps150'][0]),
        fiducials['bias'][0],
        fiducials['eps150'][1],
        fiducials['alpha150'][0],
        fiducials['alpha150'][1],
        fiducials['alpha110'][0],
        fiducials['alpha110'][1],
        fiducials['EW'][0],
        fiducials['EW'][1],
        #fiducials['fescape'][0],
        #fiducials['fescape'][1],
        fiducials['bias'][1],
        fiducials['bias'][2]
        ]

    F_N = fisher_matrix(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_F = fisher_matrix(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_U = fisher_matrix(use_pars_fid,'ULTRASAT','SPHEREx',group_vox,run)
    F = F_F#N + F_F

    F_b = F + F_U

    for j in range(len(F)):
        if use_pars_fid[j] == 'gamma_150':
            F[j,j] += 1/.3**2
            F_U[j,j] += 1/.3**2
            F_b[j,j] += 1/.3**2
        if use_pars_fid[j] == 'C_alpha_150':
            F[j,j] += 1/1.5**2
            F_U[j,j] += 1/1.5**2
            F_b[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_110':
            F[j,j] += 1/1.5**2
            F_U[j,j] += 1/1.5**2
            F_b[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'EW_0.3':
            F[j,j] += 1/10**2
            F_U[j,j] += 1/10**2
            F_b[j,j] += 1/10**2
        if use_pars_fid[j] == 'EW_1':
            F[j,j] += 1/10**2
            F_U[j,j] += 1/10**2
            F_b[j,j] += 1/10**2


    #inv_F = np.linalg.inv(F)
    inv_F_G = np.linalg.inv(F)
    inv_F_U = np.linalg.inv(F_U)
    inv_F_b = np.linalg.inv(F_b)

    eps_150 = eps150(z).value

    deps_deps0 = (1+z)**use_fid_vals[use_pars_fid.index('gamma_150')]
    deps_dgamma = pow(10,use_fid_vals[use_pars_fid.index('log10_eps_150_0')]) * ((1+z)**use_fid_vals[use_pars_fid.index('gamma_150')])*np.log(1+z)

    sigma2_eps0_G = inv_F_G[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('log10_eps_150_0')]
    sigma2_gamma_G = inv_F_G[use_pars_fid.index('gamma_150')][use_pars_fid.index('gamma_150')]
    sigma_epsgamma_G =  inv_F_G[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('gamma_150')]

    sigma2_eps0_U = inv_F_U[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('log10_eps_150_0')]
    sigma2_gamma_U = inv_F_U[use_pars_fid.index('gamma_150')][use_pars_fid.index('gamma_150')]
    sigma_epsgamma_U =  inv_F_U[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('gamma_150')]

    sigma2_eps0_b = inv_F_b[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('log10_eps_150_0')]
    sigma2_gamma_b = inv_F_b[use_pars_fid.index('gamma_150')][use_pars_fid.index('gamma_150')]
    sigma_epsgamma_b =  inv_F_b[use_pars_fid.index('log10_eps_150_0')][use_pars_fid.index('gamma_150')]

    sigma_eps150_G = np.sqrt( deps_deps0**2 * sigma2_eps0_G + deps_dgamma**2 * sigma2_gamma_G + 2*deps_deps0 * deps_dgamma * sigma_epsgamma_G )
    sigma_eps150_U = np.sqrt( deps_deps0**2 * sigma2_eps0_U + deps_dgamma**2 * sigma2_gamma_U + 2*deps_deps0 * deps_dgamma * sigma_epsgamma_U )
    sigma_eps150_b = np.sqrt( deps_deps0**2 * sigma2_eps0_b + deps_dgamma**2 * sigma2_gamma_b + 2*deps_deps0 * deps_dgamma * sigma_epsgamma_b )

    # print('\nsigma eps150 GALEX = ' + str(sigma_eps150_G))
    # print('\nsigma eps150 ULTRASAT = ' + str(sigma_eps150_U))
    # print('\nsigma eps150 both = ' + str(sigma_eps150_b))

    #plt.subplot(511)
    plt.plot(z, eps_150,color=color_ULTRASAT)
    plt.fill_between(z, eps_150 - sigma_eps150_G, eps_150 + sigma_eps150_G, color=color_FUV, alpha = 0.2, label=r'$\rm GALEX\,FUV$')

    plt.fill_between(z, eps_150 - sigma_eps150_U, eps_150 + sigma_eps150_U, color=color_ULTRASAT, alpha = 0.2, label=r'$\rm ULTRASAT$')

    plt.plot(z, eps_150 - sigma_eps150_b,  color=color_ULTRASAT, linestyle = '--')
    plt.plot(z, eps_150 + sigma_eps150_b, color=color_ULTRASAT, linestyle = '--', label=r'$\rm FUV+ULTRASAT$')

    #plt.ylim(1e25,1e27)
    plt.ylabel(r'$\epsilon_{150}$')
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$z$')
    plt.ylim(1e25,1e28)
    plt.xlim(z[0],z[-1])
    if scale_physical_max == 300*u.Mpc: 
        filename = 'PLOTS/ULTRASAT/EBL/eps_150_err_combined_thetamax.png'
    else:
        filename = 'PLOTS/ULTRASAT/EBL/eps_150_err_combined.png'

    plt.savefig(filename)
    plt.show()
    return 
    alpha_150 = alpha150(z)

    dalpha150_dalpha150_0 = 1.
    dalpha150_dC150 = np.log10(1+z)
    
    sigma2_alpha150_0 = inv_F[use_pars_fid.index('alpha_150_0')][use_pars_fid.index('alpha_150_0')]
    sigma2_C150 = inv_F[use_pars_fid.index('C_alpha_150')][use_pars_fid.index('C_alpha_150')]
    sigma_alphaC150 =  inv_F[use_pars_fid.index('alpha_150_0')][use_pars_fid.index('C_alpha_150')]

    sigma_alpha150 = np.sqrt( dalpha150_dalpha150_0**2 * sigma2_alpha150_0 + dalpha150_dC150**2 * sigma2_C150 + 2*dalpha150_dalpha150_0 * dalpha150_dC150 * sigma_alphaC150 )


    #plt.subplot(512)
    plt.plot(z, alpha_150, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.fill_between(z, alpha_150 - sigma_alpha150, alpha_150 + sigma_alpha150, color='r', alpha = 0.2)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{150}$')
    plt.xlabel(r'$z$')
    plt.xlim(z[0],z[-1])
    plt.savefig('PLOTS/ULTRASAT/EBL/alpha_150_err.png')

    plt.show()

    alpha_110 = alpha110(z)

    dalpha110_dalpha110_0 = 1.
    dalpha110_dC110 = np.log10(1+z)
    
    sigma2_alpha110_0 = inv_F[use_pars_fid.index('alpha_110_0')][use_pars_fid.index('alpha_110_0')]
    sigma2_C110 = inv_F[use_pars_fid.index('C_alpha_110')][use_pars_fid.index('C_alpha_110')]
    sigma_alphaC110 =  inv_F[use_pars_fid.index('alpha_110_0')][use_pars_fid.index('C_alpha_110')]

    sigma_alpha110 = np.sqrt( dalpha110_dalpha110_0**2 * sigma2_alpha110_0 + dalpha110_dC110**2 * sigma2_C110 + 2*dalpha110_dalpha110_0 * dalpha110_dC110 * sigma_alphaC110 )

    #plt.subplot(513)
    plt.plot(z, alpha_110, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.fill_between(z, alpha_110 - sigma_alpha110, alpha_110 + sigma_alpha110, color='r', alpha = 0.2)
    plt.ylim(-6,6)
    plt.ylabel(r'$\alpha_{110}$')
    plt.xlabel(r'$z$')
    plt.xlim(z[0],z[-1])
    plt.savefig('PLOTS/ULTRASAT/EBL/alpha_110_err.png')

    plt.show()

    EW = EW_val(z).value
    dEW_dEW03 = 1 - np.log10((1+z)/(1+0.3))/np.log10((1+1)/(1+0.3))
    dEW_dEW1 = np.log10((1+z)/(1+0.3))/np.log10((1+1)/(1+0.3))
    
    sigma2_EW03 = inv_F[use_pars_fid.index('EW_0.3')][use_pars_fid.index('EW_0.3')]
    sigma2_EW1 = inv_F[use_pars_fid.index('EW_1')][use_pars_fid.index('EW_1')]
    sigma_EW03EW1 =  inv_F[use_pars_fid.index('EW_0.3')][use_pars_fid.index('EW_1')]

    sigma_EW = np.sqrt( dEW_dEW03**2 * sigma2_EW03 + dEW_dEW1**2 * sigma2_EW1 + 2*dEW_dEW03 * dEW_dEW1 * sigma_EW03EW1 ) 

    print(EW, sigma_EW)
    #plt.subplot(514)
    plt.plot(z, EW, label=r'$\lambda = %g\,{\rm nm}$'%150)
    plt.fill_between(z, EW - sigma_EW, EW + sigma_EW, color='r', alpha = 0.2)
    plt.ylim(-10,30)
    plt.ylabel(r'$\rm EW_{Ly\alpha} [A]$')

    plt.xlabel(r'$z$')
    plt.xlim(z[0],z[-1])
    plt.savefig('PLOTS/ULTRASAT/EBL/EW_err.png')

    plt.show()

    return


def error_emissivity_biasfix(group_vox,run = False):

    # !!! WE FIX THE BIAS 

    z_val = [0.5]#.3,0.6,1,1.5] # 

    rest_wave = np.concatenate((np.linspace(91.2,121.6,200), np.linspace(121.6,300,200)))*u.nm

    use_pars_fid = ['log10_eps_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110','EW_0.3','EW_1',
                    #'log_fLyC_1','log_fLyC_2',
                    # 'bias_150_0',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps150'][0]),
        fiducials['eps150'][1],
        fiducials['alpha150'][0],
        fiducials['alpha150'][1],
        fiducials['alpha110'][0],
        fiducials['alpha110'][1],
        fiducials['EW'][0],
        fiducials['EW'][1],
        #fiducials['fescape'][0],
        #fiducials['fescape'][1],
        #fiducials['bias'][0],
        fiducials['bias'][1],
        fiducials['bias'][2]
        ]
    
    #if detector == 'ULTRASAT':
    F_U = fisher_matrix(use_pars_fid,'ULTRASAT','SPHEREx',group_vox,run)
    F_F = fisher_matrix(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    F_C = fisher_matrix(use_pars_fid,'CASTOR_G','SDSS',group_vox,run)
    F_N = fisher_matrix(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    F_all = F_F + F_N + F_U #+ F_C

    for j in range(len(F_U)):
        if use_pars_fid[j] == 'gamma_150':
            F_U[j,j] += 1/.3**2
            F_F[j,j] += 1/.3**2
            F_N[j,j] += 1/.3**2
            F_all[j,j] += 1/.3**2
        if use_pars_fid[j] == 'C_alpha_150':
            F_F[j,j] += 1/1.5**2
            F_N[j,j] += 1/1.5**2
            F_U[j,j] += 1/1.5**2
            F_all[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_110':
            F_F[j,j] += 1/1.5**2
            F_N[j,j] += 1/1.5**2
            F_U[j,j] += 1/1.5**2
            F_all[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'EW_0.3':
            F_F[j,j] += 1/1**2
            F_N[j,j] += 1/1**2
            F_U[j,j] += 1/1**2
            F_all[j,j] += 1/1**2
        if use_pars_fid[j] == 'EW_1':
            F_F[j,j] += 1/1**2
            F_N[j,j] += 1/1**2
            F_U[j,j] += 1/1**2
            F_all[j,j] += 1/1**2

    #else:
    #    F_N = fisher_matrix(use_pars_fid,'GALEX_NUV','SDSS',group_vox,run)
    #    F_F = fisher_matrix(use_pars_fid,'GALEX_FUV','SDSS',group_vox,run)
    #    F_all = F_F #+ F_F


    inv_F_FUV_all = np.linalg.inv(F_F+F_N)
    inv_F_UL_all = np.linalg.inv(F_U)
    inv_F_all = np.linalg.inv(F_all)
 
    # remove the 2 bias parameters
    inv_F_FUV = np.zeros((len(inv_F_FUV_all)-2,len(inv_F_FUV_all)-2))
    inv_F_UL = np.zeros((len(inv_F_UL_all)-2,len(inv_F_UL_all)-2))
    inv_F = np.zeros((len(inv_F_all)-2,len(inv_F_all)-2))

    for i in range(len(inv_F)):
        for j in range(len(inv_F)):
            inv_F[i,j] = inv_F_all[i,j]

    for i in range(len(inv_F_FUV)):
        for j in range(len(inv_F_FUV)):
            inv_F_FUV[i,j] = inv_F_FUV_all[i,j]

    for i in range(len(inv_F_UL)):
        for j in range(len(inv_F_UL)):
            inv_F_UL[i,j] = inv_F_UL_all[i,j]

    #if detector == 'GALEX_NUV':
    #    nu_min = nu_min_gNUV
    #    nu_max = nu_max_gNUV
#
    #elif detector == 'GALEX_FUV':
    nu_min_FUV = nu_min_gFUV
    nu_max_FUV = nu_max_gFUV

    #elif detector == 'ULTRASAT':
    nu_min = nu_min_US
    nu_max = nu_max_US
    #else:
    ##    print('Check detector!')
    #    return 

    for use_z in z_val:

        use_signal = lambda w: signal(w,use_z,False,False,False,False,False,False)

        depsnu_deps150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) / pow(10,use_fid_vals[use_pars_fid.index('log10_eps_150_0')])

        depsnu_dgamma150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * np.log10(1+use_z)

        coef_alpha150 = lambda w: np.log(nu_from_lambda(w) / nu_from_lambda(150*u.nm)) if w.value > 121.6 else np.log(nu_from_lambda(121.6*u.nm) / nu_from_lambda(150*u.nm))

        depsnu_dalpha150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w)

        depsnu_dC150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w) * np.log10(use_z)

        depsnu_dalpha110_0 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm))).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))).value if w.value <= 121.6 else 0.


        depsnu_dC110 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 121.6 else 0.


        dEW_dEW03 = 1 - np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        dEW_dEW1 = np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        depsnu_dEW03 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW03).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm) else 0. 

        depsnu_dEW1 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW1).value  if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm) else 0. 
        
        Jacobian = lambda w: np.asarray((depsnu_deps150_0(w).value, depsnu_dgamma150(w).value, depsnu_dalpha150_0(w).value,depsnu_dC150(w).value,depsnu_dalpha110_0(w),depsnu_dC110(w),depsnu_dEW03(w),depsnu_dEW1(w)))

        signal_val = np.zeros(len(rest_wave))
        sigma_eps_nu = np.zeros(len(rest_wave))
        sigma_eps_nu_FUV = np.zeros(len(rest_wave))
        sigma_eps_nu_UL = np.zeros(len(rest_wave))
        for ww in tqdm(range(len(rest_wave))):
            
            signal_val[ww] = use_signal(rest_wave[ww]).value

            sigma_eps_nu[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F,Jacobian(rest_wave[ww])]))
            sigma_eps_nu_FUV[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F_FUV,Jacobian(rest_wave[ww])]))
            sigma_eps_nu_UL[ww] = np.sqrt(np.linalg.multi_dot([Jacobian(rest_wave[ww]),inv_F_UL,Jacobian(rest_wave[ww])]))
        
        plt.figure()
        plt.plot(rest_wave,signal_val,label=r'$z= %g$'%round(use_z,3),color=color_ULTRASAT)
        plt.fill_between(rest_wave.value, signal_val - sigma_eps_nu, signal_val - sigma_eps_nu_FUV , color = color_FUV, alpha = 0.3,label=r'$\rm NUV+FUV$')
        plt.fill_between(rest_wave.value, signal_val + sigma_eps_nu, signal_val + sigma_eps_nu_FUV , color = color_FUV, alpha = 0.3)
        #plt.fill_between(rest_wave.value, signal_val - sigma_eps_nu_UL, signal_val + sigma_eps_nu_UL , color = color_ULTRASAT, alpha = 0.3,label=r'$\rm ULTRASAT$')
        #plt.loglog(rest_wave.value, signal_val - sigma_eps_nu, color = color_ULTRASAT, linestyle='--', label = r'$\rm FUV+ULTRASAT$')
        #plt.loglog(rest_wave.value, signal_val + sigma_eps_nu, color = color_ULTRASAT, linestyle='--')
        print(sigma_eps_nu_FUV, sigma_eps_nu)
        plt.fill_between(rest_wave.value, signal_val - sigma_eps_nu, signal_val + sigma_eps_nu , color = color_ULTRASAT, alpha = 0.3,label=r'$\rm GALEX+ULTRASAT$')#+ {\it g}$')
        plt.legend(loc=4)
        plt.xlim(rest_wave[0].value,rest_wave[-1].value)
        plt.ylim(1e25,1e27)
        plt.xlabel(r'$\lambda\,[{\rm nm}]$',fontsize=fontsize)
        plt.ylabel(r'$\epsilon(\nu,z)$',fontsize = fontsize)

        plt.axvspan(rest_wave[0].value,wavelenght_min('GALEX_FUV').value/(1+use_z),color='k',alpha=0.1)
        plt.axvspan(wavelenght_max('GALEX_NUV').value/(1+use_z),rest_wave[-1].value,color='k',alpha=0.1)
        plt.axvline(wavelenght_max('ULTRASAT').value/(1+use_z),3e25,3e27,color='k',linestyle='--')
        plt.axvline(wavelenght_min('ULTRASAT').value/(1+use_z),3e25,3e27,color='k',linestyle='--')

        plt.yscale('log')
        plt.tight_layout()

        if scale_physical_max == 300*u.Mpc: 
            filename = 'PLOTS/ULTRASAT/workshop/err_emissivity_z' + str(use_z) + '_thetamax.png'
        else:
            filename = 'PLOTS/ULTRASAT/workshop/err_emissivity_z' + str(use_z) + '.png'

        plt.savefig(filename,bbox_inches='tight')
        plt.show()

    return 



def error_dJ_biasfix(group_vox,run = False):

    # !!! WE FIX THE BIAS 

    z_val = list(z_small) 

    rest_wave = np.concatenate((np.linspace(91.2,121.6,50), np.linspace(121.6,300,50)))*u.nm

    use_pars_fid = ['log10_eps_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110','EW_0.3','EW_1',
                    #'log_fLyC_1','log_fLyC_2','bias_150_0',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps150'][0]),
        fiducials['eps150'][1],
        fiducials['alpha150'][0],
        fiducials['alpha150'][1],
        fiducials['alpha110'][0],
        fiducials['alpha110'][1],
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
    F_C = fisher_matrix(use_pars_fid,'CASTOR_G','SDSS',group_vox,run)
    F_all = F_N + F_F
    F_UF = F_U + F_F + F_N 

    for j in range(len(F_UF)):
        if use_pars_fid[j] == 'gamma_150':
            F_UF[j,j] += 1/.3**2
            F_all[j,j] += 1/.3**2
        if use_pars_fid[j] == 'C_alpha_150':
            F_UF[j,j] += 1/1.5**2
            F_all[j,j] += 1/1.5**2
        if use_pars_fid[j] == 'C_alpha_110':
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
    for use_z in z_val:

        use_signal = lambda w: signal(w,use_z,False,False,False,False,False,False)

        depsnu_deps150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) / pow(10,use_fid_vals[use_pars_fid.index('log10_eps_150_0')])

        depsnu_dgamma150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * np.log10(1+use_z)

        coef_alpha150 = lambda w: np.log(nu_from_lambda(w) / nu_from_lambda(150*u.nm)) if w.value > 121.6 else np.log(nu_from_lambda(121.6*u.nm) / nu_from_lambda(150*u.nm))

        depsnu_dalpha150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w)

        depsnu_dC150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w) * np.log10(use_z)

        depsnu_dalpha110_0 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm))).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))).value if w.value <= 121.6 else 0.


        depsnu_dC110 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 121.6 else 0.


        dEW_dEW03 = 1 - np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        dEW_dEW1 = np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
        depsnu_dEW03 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW03).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 

        depsnu_dEW1 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW1).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 
        
        Jacobian = lambda w: np.asarray((depsnu_deps150_0(w).value, depsnu_dgamma150(w).value, depsnu_dalpha150_0(w).value,depsnu_dC150(w).value,depsnu_dalpha110_0(w),depsnu_dC110(w),depsnu_dEW03(w),depsnu_dEW1(w)))

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

        sigma_dJ_FN[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_FN,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + 'GALEX_FUV' + '.dat')

        sigma_dJ_FU[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_FU,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + 'GALEX_FUV' + '.dat')

        dJ_F[z_val.index(use_z)] = dJdz(use_z,'GALEX_FUV',False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + 'GALEX_FUV' + '.dat') * bJ_z(use_z,'GALEX_FUV',run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + 'GALEX_FUV' + '.dat')
                
        sigma_dJ_U[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_U,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,'ULTRASAT',run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + 'ULTRASAT' + '.dat')

        dJ_U[z_val.index(use_z)] = dJdz(use_z,'ULTRASAT',False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + 'ULTRASAT' + '.dat') * bJ_z(use_z,'ULTRASAT',run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + 'ULTRASAT' + '.dat')
        
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

    if scale_physical_max == 300*u.Mpc: 
        filename = 'PLOTS/ULTRASAT/workshop/err_bJdJ_thetamax.png'
    else:
        filename = 'PLOTS/ULTRASAT/workshop/err_bJdJ.png'

    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return sigma_dJ_FU

# !!!!
# !!!! TO FIX 
def error_dJ_ULTRASAT(group_vox,run = False):

    detector_all = ['GALEX_NUV','GALEX_FUV','ULTRASAT']
    z_val = list(z_small)

    rest_wave = np.concatenate((np.linspace(91.2,121.6,200), np.linspace(121.6,350,200)))*u.nm

    use_pars_fid = ['log10_eps_150_0','gamma_150','alpha_150_0','C_alpha_150','alpha_110_0','C_alpha_110','EW_0.3','EW_1',#'log_fLyC_1','log_fLyC_2',
                    'bias_150_0',
                    'gamma_bv','gamma_bz']
    #
    use_fid_vals = [np.log10(fiducials['eps150'][0]),
        fiducials['eps150'][1],
        fiducials['alpha150'][0],
        fiducials['alpha150'][1],
        fiducials['alpha110'][0],
        fiducials['alpha110'][1],
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
    #    if use_pars_fid[i] == 'gamma_150':
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

            use_signal = lambda w: signal(w,use_z,False,False,False,False,False,False)

            depsnu_deps150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) / pow(10,use_fid_vals[use_pars_fid.index('log10_eps_150_0')])

            depsnu_dgamma150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * np.log10(1+use_z)

            coef_alpha150 = lambda w: np.log(nu_from_lambda(w) / nu_from_lambda(150*u.nm)) if w.value > 121.6 else np.log(nu_from_lambda(121.6*u.nm) / nu_from_lambda(150*u.nm))

            depsnu_dalpha150_0 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w)

            depsnu_dC150 = lambda w: signal(w,use_z,False,False,False,False,False,False) * coef_alpha150(w) * np.log10(use_z)

            depsnu_dalpha110_0 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm))).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))).value if w.value <= 121.6 else 0.


            depsnu_dC110 = lambda w: (signal(w,use_z,False,False,False,False,False,False) * np.log(nu_from_lambda(91.2*u.nm) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 91.2 else (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * ( nu_from_lambda(w) / nu_from_lambda(121.6*u.nm))**alpha110(use_z,False) * np.log(nu_from_lambda(w) / nu_from_lambda(121.6*u.nm)) * np.log10(1+use_z)).value if w.value <= 121.6 else 0.


            dEW_dEW03 = 1 - np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
            dEW_dEW1 = np.log10((1+use_z)/(1+0.3))/np.log10((1+1)/(1+0.3))
            depsnu_dEW03 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW03).value if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 

            depsnu_dEW1 = lambda w: (non_ionizing_continuum(nu_from_lambda(121.6*u.nm),use_z,False,False) * nu_from_lambda(w)**2 / (cu.c.to(u.nm/u.s)) * dEW_dEW1).value  if nu_from_lambda(w) == nu_from_lambda(121.6*u.nm)  else 0. 
            
            Jacobian = lambda w: np.asarray((depsnu_deps150_0(w).value, depsnu_dgamma150(w).value, depsnu_dalpha150_0(w).value,depsnu_dC150(w).value,depsnu_dalpha110_0(w),depsnu_dC110(w),depsnu_dEW03(w),depsnu_dEW1(w)))

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
            
                sigma_dJ[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')

                sigma_dJ_det[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_det,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')

                dJ[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')
            
                dJ_eps[z_val.index(use_z)] = dJdz(use_z,detector,True,eps_large,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '_eps150_0_large.dat') * bJ_z(use_z,detector,run=True,vals_eps150=eps_large,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '_eps150_0_large.dat')
            
                dJ_bias[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=True,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = bias_large,filename='results/map_EBL/bJ_' + detector + '_b150_0_large.dat')
            
            else:
                sigma_dJ[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')

                sigma_dJ_det[z_val.index(use_z)] = (((cu.c.to(u.km/u.s) / (4*np.pi*H(use_z)*(1+use_z)) * np.trapz(intg_det,nu_obs_val))*unit*u.Hz/u.steradian).to(u.Jy/u.steradian)).value * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')

                dJ[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=False,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '.dat')

                dJ_eps[z_val.index(use_z)] = dJdz(use_z,detector,True,eps_large,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '_eps150_0_large.dat') * bJ_z(use_z,detector,run=True,vals_eps150=eps_large,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = False,filename='results/map_EBL/bJ_' + detector + '_eps150_0_large.dat')
            
                dJ_bias[z_val.index(use_z)] = dJdz(use_z,detector,False,False,False,False,False,False,False,filename='results/map_EBL/dJdz_' + detector + '.dat') * bJ_z(use_z,detector,run=True,vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,val_bias = bias_large,filename='results/map_EBL/bJ_' + detector + '_b150_0_large.dat')

        plt.subplot(subplot[detector_all.index(detector)])

        detector_name = r'$GALEX\, NUV$' if detector == 'GALEX_NUV' else r'$GALEX\, FUV$' if detector == 'GALEX_FUV' else r'$ULTRASAT$' if detector == 'ULTRASAT' else -1
        plt.plot(z_val,dJ,'k',linestyle='-',label=detector_name)
        plt.plot(z_val,dJ_eps,'k',linestyle='--',label=r'$3\epsilon_{150}^0$')
        plt.plot(z_val,dJ_bias,'k',linestyle='-.',label=r'$3b_{150}^0$')
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
        filename = 'PLOTS/ULTRASAT/EBL/forecast_dJ_all_beps_thetamax.png'
    else:
        filename = 'PLOTS/ULTRASAT/EBL/forecast_dJ_all_beps.png'
    plt.savefig(filename,bbox_inches='tight')
    plt.show()

    return sigma_dJ
