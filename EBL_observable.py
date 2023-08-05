from EBL_signal import *
from scipy.special import jv
from scipy.interpolate import CubicSpline
from scipy import signal as scipy_signal

scale_physical_min = 0.5*u.Mpc
scale_physical_max = 300*u.Mpc #5*u.Mpc
theta_min = lambda angular_distance: np.arctan(scale_physical_min.value / angular_distance) 
theta_max = lambda angular_distance: np.arctan(scale_physical_max.value / angular_distance)


def Jnu_monopole(detector):

    use_z = z_small if detector == 'GALEX_NUV' or detector == 'GALEX_FUV' or detector == 'ULTRASAT' else z_small_castor if detector == 'CASTOR_UV' or detector == 'CASTOR_U' or detector == 'CASTOR_G' else -1

    filename = 'dat/dJdz_' + detector + '.dat'
    
    J = np.zeros(len(use_z))
    for i in range(len(use_z)):

        J[i] = dJdz(use_z[i], detector,run = False,\
         vals_eps150=False,vals_alpha150=False,vals_alpha110=False,val_EW=False,val_flyc=False,val_alpha90=False,filename=filename)

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

    plt.savefig('PLOTS/ULTRASAT/EBL/bessel.png',bbox_inches='tight')
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


    zv = [z_SDSS[0],z_SDSS[50],z_SDSS[100],z_SDSS[500],z_SDSS[-1]] 
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
    plt.savefig('PLOTS/ULTRASAT/EBL/kernel.png',bbox_inches='tight')
    plt.show()

    return 


def w_thetaz_m(theta,z):
    #,k_low,ktheta_stop):

    kernel_z = lambda k: (bessel(k*theta*chi(z)))  /(2*np.pi) * dz_dchi(z) 

    intg = lambda k: (k * Pkz(k/u.Mpc,z) * kernel_z(k)).value

    #unit = k_arr[0].unit*k_arr[0].unit*Pkz(k_arr[0],0).unit * dz_dchi(0.).unit ---> adimensional, we don't put it to save time
    
    w = quad(intg, k_arr[0].value,k_arr[-1].value)[0]


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
        plt.axvline(theta_max(cosmo.angular_diameter_distance(0.5))/0.000290888,ymin=0,ymax=1,color='k',linewidth=1.5)
        plt.axvline(theta_min(cosmo.angular_diameter_distance(1))/0.000290888,0,1,color='k',linestyle='--',linewidth=1.5)
        plt.axvline(theta_max(cosmo.angular_diameter_distance(1))/0.000290888,0,1,color='k',linestyle='--',linewidth=1.5)

    plt.xlabel(r'$\theta\,[{\rm arcmin}]$',fontsize=fontsize)
    plt.ylabel(r'$|w(\theta,z)|$',fontsize=fontsize)
    plt.legend()
    plt.ylim(1e-6,3e-3)

    if scale_physical_max == 300*u.Mpc:
        filefig = 'PLOTS/ULTRASAT/EBL/wthetaz_thetamax.png'
    else:
        filefig = 'PLOTS/ULTRASAT/EBL/wthetaz.png'

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

            Nval[i] = nval[i]*V/(2*halfbin)

        dNdz = interp1d(zval,Nval)(z) if z >= zval[0] and z <= zval[-1] else Nval[0] if z < zval[0] else 0.

    else: 
       print('Galaxy survey not recognized!')
       return -1 

    return dNdz

def plot_dNgdz():

    z = z_small 
    N = np.zeros(len(z))
    Nsp =  np.zeros(len(z))

    for i in range(len(z)):
        N[i] = dNgdz(z[i],'SDSS')#/norm_SDSS
        Nsp[i] = dNgdz(z[i],'SPHEREx')#/norm_SPHEREx

    print(np.trapz(N,z))
    print(np.trapz(Nsp,z))
    #zval_main, Nval_main = np.genfromtxt('dat/SDSS/MAIN_dNdz.dat').T
    #zval_lowz, Nval_lowz = np.genfromtxt('dat/SDSS/LOWZ_dNdz.dat').T
    #zval_cmass, Nval_cmass = np.genfromtxt('dat/SDSS/CMASS_dNdz.dat').T
    #zval_qso, Nval_qso = np.genfromtxt('dat/SDSS/QSO_dNdz.dat').T

    #plt.plot(zval_main,Nval_main/0.1,'--',label=r'$MAIN$')
    #plt.plot(zval_lowz,Nval_lowz/0.1,'--',label=r'$LOWZ$')
    #plt.plot(zval_cmass,Nval_cmass/0.1,'--',label=r'$CMASS$')
    #plt.plot(zval_qso,Nval_qso/0.1,'--',label=r'$QSO$')
    plt.plot(z,N,label=r'$\rm SDSS$',color=color_FUV)
    plt.plot(z,Nsp,color_ULTRASAT,label=r'$\rm SPHEREx$')
    plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$dN_g/d z$',fontsize=fontsize)
    plt.legend(loc=1)
    plt.tight_layout()
    
    plt.savefig('PLOTS/ULTRASAT/EBL/dNg.png',bbox_inches='tight')
    plt.show()

    return 

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

    else: 
       print('Galaxy survey not recognized!')
       return -1 

    return b


def plot_bgal():

    z = z_SDSS
    b = np.zeros(len(z))
    bsp = np.zeros(len(z))
    for i in range(len(z)):
        b[i] = b_gal(z[i],'SDSS')
        bsp[i] = b_gal(z[i],'SPHEREx')

    #zval, dNval, bval = np.genfromtxt('dat/SPHEREx.dat')
#
    #zval_main, bval_main = np.genfromtxt('dat/SDSS/MAIN_bias.dat').T
    #zval_lowz, bval_lowz = np.genfromtxt('dat/SDSS/LOWZ_bias.dat').T
    #zval_cmass, bval_cmass = np.genfromtxt('dat/SDSS/CMASS_bias.dat').T
    #zval_qso, bval_qso = np.genfromtxt('dat/SDSS/QSO_bias.dat').T

    #plt.plot(zval_main,bval_main,'--',label=r'$MAIN$')
    #plt.plot(zval_lowz,bval_lowz,'--',label=r'$LOWZ$')
    #plt.plot(zval_cmass,bval_cmass,'--',label=r'$CMASS$')
    #plt.plot(zval_qso,bval_qso,'--',label=r'$QSO$')
    plt.plot(z,b,label=r'$\rm SDSS$',color=color_FUV)
    plt.plot(z,bsp,color=color_ULTRASAT,label=r'$\rm SPHEREx$')
#    plt.plot(zval,bval,'D--',color=color_ULTRASAT)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$b_g$',fontsize=fontsize)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig('PLOTS/ULTRASAT/EBL/bg.png',bbox_inches='tight')
    plt.show()

    return 




Window = lambda theta: theta**-0.8
Window_normed = lambda theta: Window(theta) / (quad(Window,0,np.pi)[0])

def intg_f(tv,z):
    return Window_normed(tv) *  w_thetaz_m(tv,z)


def wz_m(z, run = False, filename = ''):
    

    if run:
        theta_min_use = theta_min(cosmo.angular_diameter_distance(z)) 
        theta_max_use = theta_max(cosmo.angular_diameter_distance(z))

        theta = np.linspace(theta_min_use,theta_max_use,200)

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
            vals_eps150=False,
            vals_alpha150=False,
            vals_alpha110=False,
            val_EW=False,
            val_flyc=False,
            val_alpha90=False,
            val_bias=False,
            filename = '',
            filename_wm = '',
            filename_dJ = '',
            filename_bJ = ''):

    if not run and os.path.exists(filename):

        z_arr, bar_wJg_all = np.genfromtxt(filename)
        bar_wJg = interp1d(z_arr,bar_wJg_all)(z)

    else:

        bJdJdz = lambda zv: bJ_z(zv, detector,False,vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90,val_bias,filename_bJ) * dJdz(zv, detector, False, vals_eps150,vals_alpha150,vals_alpha110,val_EW,val_flyc,val_alpha90,filename_dJ)#.value

        bar_wm = lambda zv: wz_m(zv,False,filename_wm)

        #dN = lambda zv: dNgdz(zv,gal_survey) 
        #A = quad(dN,0.,5.)[0]
        #dN_norm = lambda zv: dN(zv)/A

        bdNg = lambda zv: b_gal(zv,gal_survey) #*  dN_norm(zv) 

        intg = lambda zv: bar_wm(zv) * bdNg(zv) * bJdJdz(zv)

        bar_wJg = intg(z) 

    return bar_wJg


def run_wm():

    wm = np.zeros(len(z_small))
    if scale_physical_max == 300*u.Mpc:
        filename = 'dat/wmz_thetamax.dat'
    else:
        filename = 'dat/wmz.dat'

    for i in (range(len(z_small))):
        print('\nDoing z = ' + str(z_small[i]))
        print('DM')
        wm[i] = wz_m(z_small[i],run = True,filename =filename)

    np.savetxt(filename,(z_small, wm))

#    wm = np.zeros(len(z_small_castor))
#    for i in (range(len(z_small_castor))):
#        print('\nDoing z = ' + str(z_small_castor[i]))
#        print('DM')
#        wm[i] = wz_m(z_small_castor[i],run = True,filename = 'dat/wmz_castor.dat')
#
#    np.savetxt('dat/wmz_castor.dat',(z_small_castor, wm))

    return 

def run_wJg():

    wn = np.zeros(len(z_small))
    wf = np.zeros(len(z_small))
    wu = np.zeros(len(z_small))
    wcuv = np.zeros(len(z_small_castor))
    wcu = np.zeros(len(z_small_castor))
    wcg = np.zeros(len(z_small_castor))

    if scale_physical_max == 300*u.Mpc:
        filename_dm = 'dat/wmz_thetamax.dat'
    else:
        filename_dm = 'dat/wmz.dat'

    for i in (range(len(z_small))):
        print('\nDoing z = ' + str(z_small[i]))
        print('NUV')
        wn[i] = wJgz(z_small[i],'GALEX_NUV','SDSS',True,
            filename_wm = 'dat/wmz.dat',
            filename_dJ = 'dat/dJdz_GALEX_NUV.dat',
            filename_bJ = 'dat/bJ_GALEX_NUV.dat')
        print('FUV')
        wf[i] = wJgz(z_small[i],'GALEX_FUV','SDSS',True,
            filename_wm = 'dat/wmz.dat',
            filename_dJ = 'dat/dJdz_GALEX_FUV.dat',
            filename_bJ = 'dat/bJ_GALEX_FUV.dat')
        print('ULTRASAT')
        wu[i] = wJgz(z_small[i],'ULTRASAT','SPHEREx',True,
            filename_wm = filename_dm,
            filename_dJ = 'dat/dJdz_ULTRASAT.dat',
            filename_bJ = 'dat/bJ_ULTRASAT.dat')
        
    np.savetxt('dat/wJg_GALEX_NUV,SDSS.dat',(z_small, wn))
    np.savetxt('dat/wJg_GALEX_FUV,SDSS.dat',(z_small, wf))
    if scale_physical_max == 300*u.Mpc:
        filename_ULT = 'dat/wJg_ULTRASAT,SPHEREx_thetamax.dat'
    else:
        filename_ULT = 'dat/wJg_ULTRASAT,SPHEREx.dat'

    np.savetxt(filename_ULT,(z_small, wu))
    
    for i in (range(len(z_small_castor))):
        print('\nDoing z = ' + str(z_small_castor[i]))
        print('CASTOR UV')
        wcuv[i] = wJgz(z_small_castor[i],'CASTOR_UV','SDSS',True,
            filename_wm = 'dat/wmz_castor.dat',
            filename_dJ = 'dat/dJdz_CASTOR_UV.dat',
            filename_bJ = 'dat/bJ_CASTOR_UV.dat')
        print('CASTOR U')
        wcu[i] = wJgz(z_small_castor[i],'CASTOR_U','SDSS',True,
            filename_wm = 'dat/wmz_castor.dat',
            filename_dJ = 'dat/dJdz_CASTOR_U.dat',
            filename_bJ = 'dat/bJ_CASTOR_U.dat')
        print('CASTOR G')
        wcg[i] = wJgz(z_small_castor[i],'CASTOR_G','SDSS',True,
            filename_wm = 'dat/wmz_castor.dat',
            filename_dJ = 'dat/dJdz_CASTOR_G.dat',
            filename_bJ = 'dat/bJ_CASTOR_G.dat')
        
    np.savetxt('dat/wJg_CASTOR_UV,SDSS.dat',(z_small_castor, wcuv))
    np.savetxt('dat/wJg_CASTOR_U,SDSS.dat',(z_small_castor, wcu))
    np.savetxt('dat/wJg_CASTOR_G,SDSS.dat',(z_small_castor, wcg))

    return 



def plot_wz():

    wm = np.zeros(len(z_small))
    wn = np.zeros(len(z_small))
    wf = np.zeros(len(z_small))
    wu = np.zeros(len(z_small))
    wcuv = np.zeros(len(z_small_castor))
    wcu = np.zeros(len(z_small_castor))
    wcg = np.zeros(len(z_small_castor))

    for i in (range(len(z_small))):
        print('\nDoing z = ' + str(z_small[i]))
        print('DM')
        wm[i] = wz_m(z_small[i],run = False,filename = 'dat/wmz.dat')

    if scale_physical_max == 300*u.Mpc:
        filename_ULT = 'dat/wJg_ULTRASAT,SPHEREx_thetamax.dat'
    else:
        filename_ULT = 'dat/wJg_ULTRASAT,SPHEREx.dat'

    for i in (range(len(z_small))):
        print('\nDoing z = ' + str(z_small[i]))
        print('NUV')
        wn[i] = wJgz(z_small[i],'GALEX_NUV','SDSS',False,filename='dat/wJg_GALEX_NUV,SDSS.dat')
        print('FUV')
        wf[i] = wJgz(z_small[i],'GALEX_FUV','SDSS',False,
            filename='dat/wJg_GALEX_FUV,SDSS.dat')
        print('ULTRASAT')
        wu[i] = wJgz(z_small[i],'ULTRASAT','SPHEREx',False,
            filename=filename_ULT)
        
    for i in (range(len(z_small_castor))):
        print('\nDoing z = ' + str(z_small_castor[i]))
        print('CASTOR UV')
        wcuv[i] = wJgz(z_small_castor[i],'CASTOR_UV','SDSS',False,filename='dat/wJg_CASTOR_UV,SDSS.dat')
        print('CASTOR U')
        wcu[i] = wJgz(z_small_castor[i],'CASTOR_U','SDSS',False,
            filename='dat/wJg_CASTOR_U,SDSS.dat')
        print('CASTOR G')
        wcg[i] = wJgz(z_small_castor[i],'CASTOR_G','SDSS',False,
            filename='dat/wJg_CASTOR_G,SDSS.dat')
        
    plt.figure()
    plt.plot(z_small,(wu),color_ULTRASAT,label=r'$\rm ULTRASAT\times SPHEREx$')
    plt.plot(z_small,(wn),color_NUV,label=r'$\rm NUV\times SDSS$')
    plt.plot(z_small,(wf),color_FUV,label=r'$\rm FUV\times SDSS$')
    plt.plot(z_small_castor,(wcuv),color_CASTOR,linestyle=':',label=r'$uv\times {\rm SDSS}$')
    plt.plot(z_small_castor,(wcu),color_CASTOR,linestyle='--',label=r'$u\times {\rm SDSS}$')
    plt.plot(z_small_castor,(wcg),color_CASTOR,linestyle='-.',label=r'$g\times {\rm SDSS}$')
    #plt.plot(z_small,abs(wm),'k',label=r'$\bar{\omega}_m(z)$')

    #plt.yscale('log')
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$\bar{\omega}_{J_\nu{\rm g}}(z)$',fontsize=fontsize)
    plt.legend(loc=1,ncol=2)

    plt.tight_layout()
    if scale_physical_max == 300*u.Mpc:
        filefig = 'PLOTS/ULTRASAT/EBL/wJg_thetamax.png'
    else:
        filefig = 'PLOTS/ULTRASAT/EBL/wJg.png'

    plt.savefig(filefig,bbox_inches='tight')
    plt.show()

    return 







