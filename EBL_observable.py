from EBL_signal import *
from scipy.special import jv
from scipy.interpolate import CubicSpline
from scipy import signal as scipy_signal
from matplotlib import patches

scale_physical_min = 1.*u.Mpc #0.5*u.Mpc
scale_physical_max = 5*u.Mpc
theta_min = lambda angular_distance: np.arctan(scale_physical_min.value / angular_distance) 
theta_max = lambda angular_distance, detector: ((4*u.deg).to(u.rad)).value if detector == 'ULTRASAT' else np.arctan(scale_physical_max.value / angular_distance)


def Jnu_monopole(detector,reduced_z = False):

    #if reduced_z:
    #    use_z = z_gals_interp if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1
    #else:
    gal_survey = 'SPHEREx' if detector == 'ULTRASAT' else 'SDSS'
    use_z = z_gals(gal_survey) if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else -1

    reduced_label = '_reduced' #if reduced_z else ''

    filename = use_results_dir + 'EBL/dJdz_' + detector + reduced_label + '.dat'
    
    J = np.zeros(len(use_z))
    for i in range(len(use_z)):

        J[i] = dJdz(use_z[i], detector,run = False,\
         vals_eps1500=False,vals_alpha1500=False,vals_alpha1100=False,val_EW=False,val_flyc=False,val_alpha900=False,filename=filename)

    Jnu = np.trapz(J,use_z)

    return Jnu


def Pkz(k,z):

    Pkz = PK.P(z,k.value)

    return Pkz

# comoving distance
def chi(z):

    comoving_distance = cosmo.comoving_radial_distance(z) #* u.Mpc

    return comoving_distance


def dz_dchi(z):

    E = H(z) / (light * 1e-3) # km/s/Mpc / (km/s)

    return E 


def bessel(arg):

    J0 = jv(0,arg) 

    return J0


def plot_bessel():

    z = [0.5,1,2]
    colors = ['b','teal','springgreen','orange','tomato','r']

    theta_deg = np.linspace(0.,180,100)
    theta = np.radians(theta_deg)

    for i in z:

        b1=  np.zeros(len(theta))
        b2=  np.zeros(len(theta))
        b3=  np.zeros(len(theta))
        b4=  np.zeros(len(theta))

        print('Doing z = ' + str(i))
        for t in tqdm(range(len(theta))):

            b1[t] = bessel(k_arr[0].value*theta[t]*chi(i))
            b2[t] = bessel(k_arr[50].value*theta[t]*chi(i))
            b3[t] = bessel(k_arr[500].value*theta[t]*chi(i))
            b4[t] = bessel(k_arr[-1].value*theta[t]*chi(i))

        plt.plot(theta,b1,linestyle='-',color=colors[z.index(i)],label=r'$k = %g$'%round(k_arr[0].value,4) + r'$,\,z =%g$'%i)
        plt.plot(theta,b2,linestyle='--',color=colors[z.index(i)],label=r'$k = %g$'%round(k_arr[50].value,4) + r'$,\,z =%g$'%i)
        plt.plot(theta,b3,linestyle='-.',color=colors[z.index(i)],label=r'$k = %g$'%round(k_arr[500].value,4) + r'$,\,z =%g$'%i)
        plt.plot(theta,b4,linestyle=':',color=colors[z.index(i)],label=r'$k = %g$'%round(k_arr[-1].value,4) + r'$,\,z =%g$'%i)

    plt.xlabel(r'$\theta$',fontsize=fontsize)
    plt.ylabel(r'$J_0(k\theta\chi(z))$',fontsize=fontsize)
    plt.legend()

    plt.savefig(use_results_dir + 'PLOTS/EBL/bessel.png',bbox_inches='tight')
    plt.show()

    return 


def kernel(ktheta):

    zm = 0.0741
    dNdz = lambda z: 3*z**2 / (2*(zm/1.412)**3) * np.exp(-(1.412*z/zm)**(3/2))

    intg = lambda z: bessel(ktheta*chi(z)) *dz_dchi(z).value /(2*np.pi) *(dNdz(z)/quad(dNdz,0,np.inf)[0])**2 

    ker = quad(intg,0,np.inf)[0]

    return ker

def plot_kernel():

    ktheta = np.linspace(0.001,0.1)
    ker = np.zeros(len(ktheta))
    for i in range(len(ker)):
        ker[i] = np.abs(kernel(ktheta[i]))
    
    plt.figure(figsize=(15,8))
    plt.subplot(121)
    plt.loglog(ktheta,ker)
    plt.xlabel(r'$k\theta$',fontsize=fontsize)
    plt.ylabel(r'$|ker_{\rm REF}|$',fontsize=fontsize)


    zv = [z_gals('DESI')[0],z_gals('DESI')[50],z_gals('DESI')[100],z_gals('DESI')[500],z_gals('DESI')[-1]] 
    ktheta = np.logspace(-4,2)
    plt.subplot(122)
    for zi in zv:
        ker = np.zeros(len(ktheta))
        for i in range(len(ker)):
            intg = lambda ktheta_v: bessel(ktheta_v*chi(zi)) * dz_dchi(zi).value /(2*np.pi) * dNgdz(zi,'SDSS')
            ker[i] = np.abs(intg(ktheta[i]))

        plt.loglog(ktheta,ker,label=r'$z=%g$'%zi)
    plt.xlabel(r'$k\theta$',fontsize=fontsize)
    plt.ylabel(r'$|ker|$',fontsize=fontsize)

    plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig(use_results_dir + 'PLOTS/EBL/kernel.png',bbox_inches='tight')
    plt.show()

    return 


def w_thetaz_m(theta,z):

    kernel_z = lambda k: (bessel(k*theta*chi(z)))  /(2*np.pi) * dz_dchi(z) 

    intg = lambda k: (k * Pkz(k/u.Mpc,z) * kernel_z(k)).value

    #unit = k_arr[0].unit*k_arr[0].unit*Pkz(k_arr[0],0).unit * dz_dchi(0.).unit ---> adimensional, we don't put it to save time
    
    w = quad(intg, 1e-3,k_arr[-1].value)[0]


    return w #* unit



def plot_wtz():

    z = [0.5,1]

    theta_deg = np.linspace(0,3.14,200) 
    theta = np.radians(theta_deg)
    for i in z:
        w = np.zeros(len(theta))
        print('Doing z = ' + str(i))
        
        for t in tqdm(range(len(theta))):
            w[t] = w_thetaz_m(theta[t],i)
        
        plt.loglog(theta/0.000290888,np.abs(w),label=r'$z =%g$'%i)
        plt.axvline(theta_min(cosmo.angular_diameter_distance(0.5))/0.000290888,ymin=0,ymax=1,color='k',linewidth=1.5)
        plt.axvline(theta_max(cosmo.angular_diameter_distance(0.5),'GALEX_NUV')/0.000290888,ymin=0,ymax=1,color='k',linewidth=1.5)
        plt.axvline(theta_min(cosmo.angular_diameter_distance(1))/0.000290888,0,1,color='k',linestyle='--',linewidth=1.5)
        plt.axvline(theta_max(cosmo.angular_diameter_distance(1),'ULTRASAT')/0.000290888,0,1,color='k',linestyle='--',linewidth=1.5)

    plt.xlabel(r'$\theta\,[{\rm arcmin}]$',fontsize=fontsize)
    plt.ylabel(r'$|w(\theta,z)|$',fontsize=fontsize)
    plt.legend()
    plt.ylim(1e-6,3e-3)

    if scale_physical_max == 300*u.Mpc:
        filefig = use_results_dir + 'PLOTS/EBL/wthetaz_thetamax.png'
    else:
        filefig = use_results_dir + 'PLOTS/EBL/wthetaz.png'

    plt.savefig(filefig,bbox_inches='tight')
    plt.show()

    return 


def dNgdz(z,gal_survey):

    if gal_survey == 'SDSS':
       
       if z < 0.1:
            zval, dNdzval = np.genfromtxt('dat/SDSS/MAIN_dNdz.dat').T          

       elif 0.1 <= z <= 0.2:
            zval_main, dNdzval_main = np.genfromtxt('dat/SDSS/MAIN_dNdz.dat').T
            zval_lowz, dNdzval_lowz = np.genfromtxt('dat/SDSS/LOWZ_dNdz.dat').T
            dNdzval_main = interp1d(zval_main,dNdzval_main)
            dNdzval_lowz = interp1d(zval_lowz,dNdzval_lowz)

            zval = np.linspace(max(min(zval_lowz),min(zval_main),0.1),min(max(zval_lowz),max(zval_main),0.2))
            dNdzval = dNdzval_main(zval) + dNdzval_lowz(zval)

       elif 0.2 < z < 0.4:
            zval, dNdzval = np.genfromtxt('dat/SDSS/LOWZ_dNdz.dat').T

       elif 0.4 <= z < 0.65:
            zval, dNdzval = np.genfromtxt('dat/SDSS/CMASS_dNdz.dat').T
       else:
            zval, dNdzval = np.genfromtxt('dat/SDSS/QSO_dNdz.dat').T

       dNdz = interp1d(zval,dNdzval/0.1)(z) if z >= zval[0] and z <= zval[-1] else dNdzval[0] if z < zval[0] else 0.

    elif gal_survey == 'SPHEREx':
        zval, nval, bval = np.genfromtxt('dat/SPHEREx.dat')
        Nval = np.zeros(len(zval))

        for i in range(len(zval)):
            if i < len(zval)-1:
                if zval[i] != 0.9:
                    halfbin = (zval[i+1]-zval[i])/2
                else:
                    halfbin = (zval[i]-zval[i-1])/2
            else:
                halfbin = (zval[i]-zval[i-1])/2

            V = 4*np.pi/3 * (cosmo.comoving_radial_distance(zval[i]+halfbin)**3 - cosmo.comoving_radial_distance(zval[i]-halfbin)**3)

            Nval[i] = (H(0.)/100.).value**3 *nval[i]*V/(2*halfbin)

        dNdz = interp1d(zval,Nval)(z) if z >= zval[0] and z <= zval[-1] else Nval[0] if z < zval[0] else 0.

    elif gal_survey == 'DESI':

        zval, nval = np.genfromtxt('dat/DESI.dat')
        Nval = np.zeros(len(zval))

        Adesi = 14e3 # deg^2 since table 7 in 2306.06307 gives density deg^-2
        for i in range(len(zval)):
            if i < len(zval)-1:
                if zval[i] != 0.9:
                    halfbin = (zval[i+1]-zval[i])/2
                else:
                    halfbin = (zval[i]-zval[i-1])/2
            else:
                halfbin = (zval[i]-zval[i-1])/2

            Nval[i] = nval[i]*Adesi /(2*halfbin)

        dNdz = interp1d(zval,Nval)(z) if z >= zval[0] and z <= zval[-1] else Nval[0] if z < zval[0] else 0.


    else: 
       print('Galaxy survey not recognized!')
       return -1 

    return dNdz




def b_gal(z, gal_survey):

    if gal_survey == 'SDSS':
       
       if z < 0.1:
            zval, bval = np.genfromtxt('dat/SDSS/MAIN_bias.dat').T           
            m = (bval[-1]-bval[0])/(zval[-1]-zval[0])

       elif 0.1 <= z <= 0.2:
            zval_main, bval_main = np.genfromtxt('dat/SDSS/MAIN_bias.dat').T
            zval_lowz, bval_lowz = np.genfromtxt('dat/SDSS/LOWZ_bias.dat').T
            b_main = interp1d(zval_main,bval_main)
            b_lowz = interp1d(zval_lowz,bval_lowz)

            zval = np.linspace(max(min(zval_lowz),min(zval_main),0.1),min(max(zval_lowz),max(zval_main),0.2))
            bval = (b_main(zval) + b_lowz(zval))/2.

       elif 0.2 < z < 0.4:
            zval, bval = np.genfromtxt('dat/SDSS/LOWZ_bias.dat').T

       elif 0.4 <= z < 0.65:
            zval, bval = np.genfromtxt('dat/SDSS/CMASS_bias.dat').T

       else:
           zval, bval = np.genfromtxt('dat/SDSS/QSO_bias.dat').T
           m = (bval[-1]-bval[0])/(zval[-1]-zval[0])

       if 0.1 <= z < 0.65:   
         b = CubicSpline(zval,bval,extrapolate=True)(z)           
       else:
         b = m*(z-zval[0]) + bval[0]

    elif gal_survey == 'SPHEREx': 
        zval, nval, bval = np.genfromtxt('dat/SPHEREx.dat')
        b = interp1d(zval,bval)(z) if z >= zval[0] else bval[0]

    elif gal_survey == 'DESI': # sec 6.1 in 2306.06307
        if z < 0.4:
            b = 1.34 / Dgrowth(z)
        elif 0.4 <= z < 1.1:
            b = 1.7 / Dgrowth(z)
        elif 1.1 <= z < 1.6:
            b = 0.84 / Dgrowth(z)
        elif 1.6 <= z < 2.1:
            b = 1.2 / Dgrowth(z)
        else:
            b = 1.1 / Dgrowth(z)
    else: 
       print('Galaxy survey not recognized!')
       return -1 

    return b


def plot_bgal():


    z_sdss = np.linspace(zmin_gal,zmax_gal,80) #z_gals('SDSS')
    z_spherex = z_gals('SPHEREx')
    z_desi = np.linspace(zmin_gal,zmax_gal,80) #z_gals('DESI')
    b = np.zeros(len(z_sdss))
    bsp = np.zeros(len(z_spherex))
    bdesi = np.zeros(len(z_desi))
    for i in range(len(z_sdss)):
        b[i] = b_gal(z_sdss[i],'SDSS')
    for i in range(len(z_spherex)):
        bsp[i] = b_gal(z_spherex[i],'SPHEREx')
    for i in range(len(z_desi)):
        bdesi[i] = b_gal(z_desi[i],'DESI')

    #window_size = 5  # Adjust this value based on your preference
    bsdss_smoothed = b #moving_average(b, window_size)
    
    plt.plot(z_desi,bdesi,color=color_ULTRASAT,label=r'$\rm DESI$')
    plt.plot(z_sdss[:len(bsdss_smoothed)],bsdss_smoothed,label=r'$\rm SDSS$',color=color_FUV)
    #plt.plot(z_spherex,bsp,color=color_ULTRASAT,label=r'$\rm SPHEREx$')

    plt.xlim(zmin_gal,zmax_gal)

    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_g$',fontsize=fontsize)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig(use_results_dir + 'PLOTS/EBL/bg.png',bbox_inches='tight')
    plt.show()

    return 




Window = lambda theta: theta**-0.8
Window_normed = lambda theta: Window(theta) / (quad(Window,0,np.pi)[0])

def intg_f(tv,z):
    return Window_normed(tv) *  w_thetaz_m(tv,z)


def wz_m(z, detector, run = False, filename = ''):
    
    if run:
        DA_z = cosmo.angular_diameter_distance(z)
        theta_min_use = theta_min(DA_z) 
        theta_max_use = theta_max(DA_z,detector)

        theta = np.linspace(theta_min_use,theta_max_use,200)
        # theta_min_standard =  np.arctan(0.5 / DA_z)
        # theta = np.arange(theta_min_use,theta_max_use,(theta_max_use-theta_min_standard)/199)
        
        with Pool(6) as pool:
            intg_ff = partial(intg_f, z=z)
            intg = pool.map(intg_ff, theta)

        #intg = np.zeros(len(theta))
        #for i in tqdm(range(len(theta))):
        #    intg[i] = Window_normed(theta[i]) *  w_thetaz_m(theta[i],z)

        bar_wm = np.trapz(intg,theta) 

    else:

        zv, bv = np.genfromtxt(filename)
        bar_wm = interp1d(zv,bv)(z)

    return bar_wm


def wJgz(z,detector,gal_survey,run=False,
            vals_eps1500=False,
            vals_alpha1500=False,
            vals_alpha1100=False,
            val_EW=False,
            val_flyc=False,
            val_alpha900=False,
            val_bias=False,
            filename = '',
            filename_wm = '',
            filename_dJ = '',
            filename_bJ = ''):

    if not run and os.path.exists(filename):

        z_arr, bar_wJg_all = np.genfromtxt(filename)
        bar_wJg = interp1d(z_arr,bar_wJg_all,bounds_error=False,fill_value=bar_wJg_all[-1])(z)

    else:

        bJdJdz = lambda zv: bJ_z(zv, detector,False,vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,val_bias,filename_bJ) * dJdz(zv, detector, False, vals_eps1500,vals_alpha1500,vals_alpha1100,val_EW,val_flyc,val_alpha900,filename_dJ)#.value

        bar_wm = lambda zv: wz_m(zv,detector,False,filename_wm)

        #dN = lambda zv: dNgdz(zv,gal_survey) 
        #A = quad(dN,0.,5.)[0]
        #dN_norm = lambda zv: dN(zv)/A

        bg = lambda zv: b_gal(zv,gal_survey)  

        intg = lambda zv: bar_wm(zv) * bg(zv) * bJdJdz(zv)

        bar_wJg = intg(z) 

    return bar_wJg


def run_wm(detector,reduced_z = False):

    if reduced_z:
        use_z = z_gals_interp
        reduced_label = '_reduced'
    else:
        use_z = z_gals('SPHEREx') if detector == 'ULTRASAT' else z_gals('SDSS')
        reduced_label = ''

    wm = np.zeros(len(use_z))
    filename = use_results_dir + 'EBL/wmz_' + detector + reduced_label + '.dat'

    for i in (range(len(use_z))):
        print('\nDoing z = ' + str(use_z[i]))
        print('DM')
        wm[i] = wz_m(use_z[i],detector,run = True,filename =filename)

    np.savetxt(filename,(use_z, wm))

    return wm


def run_wJg(reduced_z = False):

    if reduced_z:
        use_z = z_gals_interp
        use_z_ultrasat = z_gals_interp
        reduced_label = '_reduced'
    else:
        use_z = z_gals('SDSS')
        reduced_label = ''
        use_z_ultrasat = z_gals('DESI')
        reduced_label = ''

    wn = np.zeros(len(use_z))
    wf = np.zeros(len(use_z))
    wu = np.zeros(len(use_z_ultrasat))


    for i in (range(len(use_z))):
       print('\nDoing z = ' + str(use_z[i]))
       print('NUV')
       filename_wm = use_results_dir + 'EBL/wmz_GALEX_NUV' + reduced_label + '.dat'
       wn[i] = wJgz(use_z[i],'GALEX_NUV','SDSS',True,
           filename_wm = filename_wm,
           filename_dJ = use_results_dir + 'EBL/dJdz_GALEX_NUV' + reduced_label + '.dat',
           filename_bJ = use_results_dir + 'EBL/bJ_GALEX_NUV' + reduced_label + '.dat')
       print('FUV')
       filename_wm = use_results_dir + 'EBL/wmz_GALEX_FUV' + reduced_label + '.dat'
       wf[i] = wJgz(use_z[i],'GALEX_FUV','SDSS',True,
           filename_wm = filename_wm,
           filename_dJ = use_results_dir + 'EBL/dJdz_GALEX_FUV' + reduced_label + '.dat',
           filename_bJ = use_results_dir + 'EBL/bJ_GALEX_FUV' + reduced_label + '.dat')
#
    for i in (range(len(use_z_ultrasat))):
        print('ULTRASAT')
        filename_wm = use_results_dir + 'EBL/wmz_ULTRASAT' + reduced_label + '.dat'
        wu[i] = wJgz(use_z_ultrasat[i],'ULTRASAT','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = use_results_dir + 'EBL/dJdz_ULTRASAT' + reduced_label + '.dat',
           filename_bJ = use_results_dir + 'EBL/bJ_ULTRASAT' + reduced_label + '.dat')
        
    np.savetxt(use_results_dir + 'EBL/wJg_GALEX_NUV,SDSS' + reduced_label + '.dat',(use_z, wn))
    np.savetxt(use_results_dir + 'EBL/wJg_GALEX_FUV,SDSS' + reduced_label + '.dat',(use_z, wf))
    filename_ULT = use_results_dir + 'EBL/wJg_ULTRASAT,DESI' + reduced_label + '.dat'

    np.savetxt(filename_ULT,(use_z_ultrasat, wu))
    
    return 


def run_wJg_galexdesi(reduced_z = False):

    if reduced_z:
        use_z = z_gals_interp
        reduced_label = '_reduced'
    else:
        use_z = z_gals('DESI')
        reduced_label = ''

    wn = np.zeros(len(use_z))
    wf = np.zeros(len(use_z))

    # use the DM file with name ultrasat because it has the same z bins as desi, ma there is no info on the lim part
    for i in (range(len(use_z))):
        filename_wm = use_results_dir + 'EBL/wmz_ULTRASAT' + reduced_label + '.dat'
        wn[i] = wJgz(use_z[i],'GALEX_NUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = use_results_dir + 'EBL/dJdz_GALEX_NUV' + reduced_label + '.dat',
           filename_bJ = use_results_dir + 'EBL/bJ_GALEX_NUV' + reduced_label + '.dat')
        wf[i] = wJgz(use_z[i],'GALEX_FUV','DESI',True,
           filename_wm = filename_wm,
           filename_dJ = use_results_dir + 'EBL/dJdz_GALEX_FUV' + reduced_label + '.dat',
           filename_bJ = use_results_dir + 'EBL/bJ_GALEX_FUV' + reduced_label + '.dat')
        
    np.savetxt(use_results_dir + 'EBL/wJg_GALEX_NUV,DESI' + reduced_label + '.dat',(use_z, wn))
    np.savetxt(use_results_dir + 'EBL/wJg_GALEX_FUV,DESI' + reduced_label + '.dat',(use_z, wf))
    
    return 



def plot_wz(reduced_z = False):

    wm = np.zeros(len(z_gals('SDSS')))
    wn = np.zeros(len(z_gals('SDSS')))
    wf = np.zeros(len(z_gals('SDSS')))
    wnD = np.zeros(len(z_gals('DESI')))
    wfD = np.zeros(len(z_gals('DESI')))
    wu = np.zeros(len(z_gals('SPHEREx')))
    wuD = np.zeros(len(z_gals('DESI')))

    reduced_label = '_reduced' #if reduced_z else ''
    use_filename_wm = use_results_dir + 'EBL/wmz_GALEX_NUV' + reduced_label + '.dat'
    use_filename_nuv = use_results_dir + 'EBL/wJg_GALEX_NUV,SDSS'+ reduced_label + '.dat'
    use_filename_fuv = use_results_dir + 'EBL/wJg_GALEX_FUV,SDSS'+ reduced_label + '.dat'
    use_filename_nuv_desi = use_results_dir + 'EBL/wJg_GALEX_NUV,DESI'+ reduced_label + '.dat'
    use_filename_fuv_desi = use_results_dir + 'EBL/wJg_GALEX_FUV,DESI'+ reduced_label + '.dat'
    for i in (range(len(z_gals('SDSS')))):
        print('\nDoing z = ' + str(z_gals('SDSS')[i]))
        print('DM')
        wm[i] = wz_m(z_gals('SDSS')[i],detector='GALEX_NUV',run = False,filename = use_filename_wm)

    #filename_ULT = use_results_dir + 'EBL/wJg_ULTRASAT,SPHEREx' + reduced_label + '.dat'
#
    filename_ULT_DESI = use_results_dir + 'EBL/wJg_ULTRASAT,DESI' + reduced_label + '.dat'

    for i in (range(len(z_gals('SDSS')))):
        print('\nDoing z = ' + str(z_gals('SDSS')[i]))
        print('NUV')
        wn[i] = wJgz(z_gals('SDSS')[i],'GALEX_NUV','SDSS',False,filename=use_filename_nuv)
        print('FUV')
        wf[i] = wJgz(z_gals('SDSS')[i],'GALEX_FUV','SDSS',False,
            filename=use_filename_fuv)
        print('ULTRASAT')
    #for i in (range(len(z_gals('SPHEREx')))):
    #    wu[i] = wJgz(z_gals('SPHEREx')[i],'ULTRASAT','SPHEREx',False,
    #        filename=filename_ULT)
    for i in (range(len(z_gals('DESI')))):
        wuD[i] = wJgz(z_gals('DESI')[i],'ULTRASAT','DESI',False,
            filename=filename_ULT_DESI)
        wnD[i] = wJgz(z_gals('DESI')[i],'GALEX_NUV','DESI',False,filename=use_filename_nuv_desi)
        wfD[i] = wJgz(z_gals('DESI')[i],'GALEX_FUV','DESI',False,
            filename=use_filename_fuv_desi)
        
    window_size = 70  # Adjust this value based on your preference
    wGN_smoothed = moving_average(wn, window_size)
    wGF_smoothed = moving_average(wf, window_size)
    wU_smoothed = moving_average(wuD, window_size)

    wGND_smoothed = moving_average(wnD, window_size)
    wGFD_smoothed = moving_average(wfD, window_size)

    plt.figure(figsize=(20,14))

    #plt.plot(z_gals('SPHEREx'),(wu),color_ULTRASAT,label=r'$\rm ULTRASAT\times SPHEREx$')
    plt.plot(z_gals('SDSS')[:len(wGF_smoothed)],wGF_smoothed,color_FUV,label=r'$\rm FUV\times SDSS$')
    plt.plot(z_gals('DESI')[:len(wGFD_smoothed)],wGFD_smoothed,color_FUV,label=r'$\rm FUV\times DESI$',linestyle='--')

    plt.plot(z_gals('SDSS')[:len(wGN_smoothed)],wGN_smoothed,color_NUV,label=r'$\rm NUV\times SDSS$')
    plt.plot(z_gals('DESI')[:len(wGND_smoothed)],wGND_smoothed,color_NUV,label=r'$\rm NUV\times DESI$',linestyle='--')
    plt.plot(z_gals('DESI')[:len(wU_smoothed)],wU_smoothed,color_ULTRASAT,label=r'$\rm ULTRASAT\times DESI$')

    plt.xlabel(r'$z$',fontsize=fontsize*1.2)
    plt.legend(bbox_to_anchor = (1.005,-0.1), ncol=3,fontsize=fontsize*1.2)
    plt.ylabel(r'$\bar{\omega}_{J_\nu{\rm g}}(z)$',fontsize=fontsize*1.2)

    plt.tight_layout()
    filefig = use_results_dir + 'PLOTS/EBL/wJg' + reduced_label + '.png'

    plt.savefig(filefig,bbox_inches='tight')
    plt.show()

    return 







