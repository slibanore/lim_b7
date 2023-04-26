# SL: last update 01/17/2023

from .overall_setting import *


###########################################
# CO LINE
###########################################

# Li et al. (2016) 
def TonyLi(self,Mvec, MLpar, z):

# SFR --> LIR --> LCO
# Uses Behroozi et al. SFR(Mh) results

    # slope of logLIR/logLCO relation
    alpha = MLpar['alpha'] 
    # intercept of logLIR/logLCO relation
    beta = MLpar['beta']
    # 10^10 times SFR/LIR normalization
    dMF = MLpar['dMF']
    # read and interpolate SFR data in file
    SFR_file = MLpar['SFR_file']
    SFR = SFR_Mz_2dinterp(Mvec,z,SFR_file)

    # includes quenching from Universe Machine
    try:
        do_quench = MLpar['do_quench']
        if do_quench:
            fQint = process_fq()
            SFR *= (1-fQint(np.log10(Mvec.value),1+z))
    except:
        pass
        
    # LIR in Lsun
    LIR = SFR/(dMF*1e-10)
    
    # L'_CO in K km/s pc^2
    Lprime = (10.**-beta * LIR)**(1./alpha)
    
    # LCO in Lsun
    L = (4.9e-5*u.Lsun)*Lprime*(self.nu/(115.27*u.GHz))**3

    return L


###########################################
# Lyalpha LINE
###########################################

# Libanore, Kovetz (2023) - ULTRASAT white paper
def Lyalpha(self, Mvec, MLpar, z):

    # read and interpolate SFR data in file
    SFR_file = MLpar['SFR_file']
    SFR = SFR_Mz_2dinterp(Mvec,z,SFR_file)

    C = 1.6e42 # erg s-1
    L_Lsun = lambda L_ergsec: L_ergsec / 3.826e33 *u.Lsun

    if 'csi' in MLpar:
        csi = MLpar['csi']
        z0 = MLpar['z0']
        zeta = MLpar['zeta']
        f0 = MLpar['f0']
        SFR0 = MLpar['SFR0'] # Msun yr-1
        t = MLpar['t']

        f_esc = ((1 + np.exp(-csi*(z-z0)))**(-zeta)*(f0 + (1-f0)/(1+(SFR/SFR0)**t)))**2 # erg s-1

        L = L_Lsun(C*f_esc*SFR)

    elif 'A_lya' in MLpar:
        A_lya = MLpar['A_lya']
        B_lya = MLpar['B_lya']
        D_lya = MLpar['D_lya']

        f_esc = B_lya*(A_lya+(1+3*A_lya)/(1.5*A_lya+(D_lya*SFR)**2))

        L = L_Lsun(C*f_esc*D_lya*SFR)

    return L 


###########################################
# UV CONTINUUM EMISSION
###########################################

# Libanore, Kovetz (2023) - ULTRASAT white paper
#!!! to be checked
def UV_continuum_SFR(self, Mvec, MLpar, z):

	SFR_file = MLpar['SFR_file']
	rescale = MLpar['rescale']
	
	# Read and interpolate Behroozi SFR(M) data
	SFR_h = SFR_Mz_2dinterp(Mvec,z,SFR_file) * u.Msun * u.yr**-1

	lambda_obs = (cu.c / self.nuObs).to(u.nm) 
	use_DL = interp1d(self.zcosmo[::-1],self.cosmo.comoving_radial_distance(self.zcosmo),kind='cubic',
                            bounds_error=False,fill_value='extrapolate')
	DL = (use_DL(z)/(1.+z)) * u.Mpc
	z_160 = -1. + (cu.c / self.nuObs).to(u.nm) / (160*u.nm) 
	DL_160 = (use_DL(z_160)/(1.+z_160)) * u.Mpc
	K_UV = 0.63 * 2.5 * 1e-10 * u.Msun * u.yr**-1 * u.Lsun**-1
	K_IR = 0.63 * 1.73 * 1e-10 * u.Msun * u.yr**-1 * u.Lsun**-1
	Ms =  10**9.15 * u.Msun
	alpha = 0.97
	IRX = lambda M: (M/Ms)**alpha 
	beta = lambda M: -2.3 + 2.5/1.42 * np.log10(((M/Ms)**alpha)/1.7+1)
	N = 0.0351 - 0.0248*(z/(1+z))
	M1 = 10**(11.59 + 1.195 * (z/(1+z))) * u.Msun
	delta = - (1.376 - 0.826*(z/(1+z)))
	gamma = 0.608 + 0.329 * (z/(z+1))
	
	M = (2*N*Mvec*((Mvec/M1)**delta+(Mvec/M1)**gamma)**-1)

	L = rescale * SFR_h / (K_UV + K_IR * IRX(M)) * \
		(lambda_obs / (160*u.nm) )**(beta(M)) * \
		(DL/DL_160)**2

	return L



###########################################
# OTHER FUNCTIONS REQUIRED
###########################################

# interpolate SFR(Mh, z)
def SFR_Mz_2dinterp(M,z,SFR_file):

    # columns in the table are:
    # (1+z), log10(Mhalo/Msun), log10(SFR / (Msun/yr))

    SFR_folder = os.path.dirname(os.path.realpath(__file__)).split("source")[0] 
    try:
        x = np.loadtxt(SFR_folder+SFR_file)
    except:
        x = np.loadtxt(SFR_file)
        
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(len(zb),len(logMb),order='F')

    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=-40.)
    
    logM = np.log10((M.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)
    
    return SFR

# interpolate quenching 
def process_fq():
 
    # columns in the table are:
    # (1+z), log10(Mhalo/Msun), fQ
 
    SFR_folder = os.path.dirname(os.path.realpath(__file__)).split("source")[0] 
    files = os.listdir(SFR_folder+'qf/')
    a = []
    mat = np.zeros((32,len(files)))
    
    counter=0
    for name in files:
        a.append(float((name.split('qf_hm_a')[1]).split('.dat')[0]))
    
    inds = np.argsort(np.array(a))
    a = np.sort(np.array(a))
    
    names = []
    for i in range(len(inds)):
        names.append(files[inds[i]])
    
    for ia in range(len(a)):  
        data = np.genfromtxt('tables/qf/'+names[ia])
        mat[:,ia] = data[:,1]
        
        
    logMh = data[:,0]
    zp1 = 1/a
        
    fQ_interp = interp2d(logMh,zp1,mat.T,bounds_error=True)
    
    return fQ_interp
    
