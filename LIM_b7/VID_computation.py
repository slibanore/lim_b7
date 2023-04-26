# SL: last update 01/18/2023
from .modelling import *

class VID(LineObs):

    def __init__(self,
        
            # VID - RELATED QUANTITES
            nT=2**22,
            Ngal_max=500,
            subtract_VID_mean=True,
            linear_VID_bin=False,
            do_sigma_G = True,
            sigma_G_input = 1.6,

            **line_kwargs):
                    
        # Initiate LineObs() parameters
        LineObs.__init__(self,**line_kwargs)

        self._update_cosmo_list = self._update_cosmo_list
        
        self._vid_params = locals()
        self._vid_params.pop('self')
        self._vid_params.pop('line_kwargs')
        self._default_vid_params = get_default_params(VID.__init__)
        check_params(self._vid_params,self._default_vid_params)
        
        # Set instrument parameters
        for key in self._vid_params:
            setattr(self,key,self._vid_params[key])
        
        # Combine lim_params with vid_params
        self._input_params.update(self._vid_params)
        self._default_params.update(self._default_vid_params)

    ###########################################
    # VID QUANTITIES
    ###########################################

# 1) Intensity bins

    # edges of intensity bins
    @cached_vid_property
    def Tedge(self):

        if self.do_fast_VID:
            Te = ulinspace(self.Tmin_VID,self.Tmax_VID,self.nT+1)
        else:
            Te = ulogspace(self.Tmin_VID,self.Tmax_VID,self.nT+1)
            
        if self.subtract_VID_mean:
            return Te-self.Tmean
        else:
            return Te

    # outputs centers of histogram bins given their edges
    def binedge_to_binctr(self,binedge):

        Nedge = binedge.size        
        binctr = (binedge[0:Nedge-1]+binedge[1:Nedge])/2.
        
        return binctr

    def binctr_to_binedge(self,binctr):

        unit = binctr.unit
        binctr = binctr.value

        binedge = np.concatenate((np.array((0,(binctr[0]+ binctr[1])/2.-2*abs((binctr[0]+ binctr[1])/2.-binctr[0]))), (binctr[0:len(binctr)-1]+ binctr[1:len(binctr)])/2., np.array(((binctr[-1]+ binctr[-2])/2. +2*abs(-(binctr[-2]+ binctr[-1])/2.+binctr[-1]),0))))

        binedge = binedge[1:-1]*unit
        
        return binedge

    # centers of intensity bins
    @cached_vid_property
    def T(self):

        return self.binedge_to_binctr(self.Tedge)
        
    # widths of intensity bins
    @cached_vid_property
    def dT(self):

        return np.diff(self.Tedge)
        

# 2) number count probability distribution 

    # mean number of galaxies per voxel
    @cached_vid_property
    def Nbar(self):

        return self.nbar*self.Vvox

    # vector of galaxy number counts, from 0 to self.Ngal_max
    @cached_vid_property
    def Ngal(self):

        return np.array(range(0,self.Ngal_max+1))
        
    # rms of fluctuations in a voxel (Gaussian window)
    @cached_vid_property
    def sigma_G(self):

        if self.do_sigma_G:
            #kvalues, and kpar, kperp. Power spectrum in observed redshift
            k = np.logspace(np.log10(self.kmin.value),np.log10(self.kmax.value),self.nk)*u.Mpc**-1
            ki,mui = np.meshgrid(k,self.mu)
            Pk = self.PKint(self.z,ki.value)*u.Mpc**3
            
            kpar = ki*mui
            kperp = ki*np.sqrt(1.-mui**2.)
            
            #Gaussian window for voxel -> FT
            Wkpar2 = np.exp(-((kpar*self.sigma_par)**2).decompose())
            Wkperp2 = np.exp(-((kperp*self.sigma_perp)**2).decompose())
            Wk2 = Wkpar2*Wkperp2
            
            #Compute sigma_G
            if self.f_NL != 0. and self.hmf_pars['fNL'] != 0.:
                bias = interp1d(self.k,self.bavg,kind='linear',bounds_error=False,fill_value=(self.bavg[0],self.bavg[-1]))
                integrnd = bias(self.k)**2*Pk*Wk2*ki**2/(4.*np.pi**2)
            else:
                bias = self.bavg[-1]
                integrnd = bias**2*Pk*Wk2*ki**2/(4.*np.pi**2)
            integrnd_mu = np.trapz(integrnd,self.mu,axis=0)
            sigma = np.sqrt(np.trapz(integrnd_mu,ki[0,:]))
            
            return sigma
        else:
            return self.sigma_G_input
    

    # probability of a voxel containing N galaxies
    @cached_vid_property
    def PofN(self):
       
        #lognormal + Poisson model from Breysse et al. 2017
        
        # PDF of galaxy density field mu
        Nbar = np.trapz(self.dndL,self.L)*self.Vvox

        logMuMin = np.log10(Nbar)-20*self.sigma_G
        logMuMax = np.log10(Nbar)+5*self.sigma_G
        mu = np.logspace(logMuMin.value,logMuMax.value,10**4)
        mu2,Ngal2 = np.meshgrid(mu,self.Ngal) # Keep arrays for fast integrals
        Pln = (np.exp(-(np.log(mu2/Nbar)+self.sigma_G**2/2.)**2/(2*self.sigma_G**2)) / (mu2*np.sqrt(2.*np.pi*self.sigma_G**2.)))

        P_poiss = poisson.pmf(Ngal2,mu2)

        return np.trapz(P_poiss*Pln,mu)


# 3) Intensity PDF 

    # relate total luminosity in a voxel to observed intensity
    @cached_vid_property
    def XLT(self):
        
        if self.developer == 'decayDM':
            if self.do_Jysr:

                x = cu.c/(4.*np.pi*(self.nu)*self.H) 
            else:
                x = cu.c**3*(1+self.z)**2/(8*np.pi*cu.k_B*self.nu**3*self.H)
            return x.to(u.uK*u.Mpc**3/u.Lsun)/self.Vvox
        
        else:
            return self.CLT/self.Vvox

    # probability of observing a given intensity in a voxel which contains exactly one emitter
    @cached_vid_property
    def P1(self):

        dndL_T = interp1d(self.L,self.dndL,bounds_error=False,fill_value='extrapolate')

        if self.subtract_VID_mean:
            # ! test
            # Te = ulogspace(self.Tmin_VID,self.Tmax_VID,self.nT+1)
            # Tedge = Te - self.Tmean
            # T = self.binedge_to_binctr(Tedge)
            # LL = ((T+self.Tmean)/self.XLT).to(u.Lsun)
            LL = ((self.T+self.Tmean)/self.XLT).to(u.Lsun)
        else:
            LL = (self.T/self.XLT).to(u.Lsun)

        dndL = dndL_T(LL.value)*self.dndL.unit
        PT1 = dndL/(np.trapz(self.dndL,self.L)*self.XLT)

        # ! test
        # dndL = dndL_T(LL.value)
        # PT1_log = dndL/(quad(dndL_T,LL[0].value,LL[-1].value)[0]*self.XLT.value)*(self.T.unit)**-1

        # PT1_interp = interp1d(T+self.Tmean,PT1_log,bounds_error=False,fill_value='extrapolate')

        # PT1 = PT1_interp(self.T+self.Tmean)
        # plt.loglog(T+self.Tmean,PT1_log)
        # plt.loglog(self.T+self.Tmean,PT1)
        # plt.show()

        return PT1

    # probability of intensity between T and T+dT in a given voxel
    @cached_vid_property
    def PT(self):
        
        # does NOT include the delta function at T=0 from voxels containing zero sources (use self.PT_zero)

        if self.do_fast_VID:
            if self.developer == 'decayDM':
                
                print('Using PT for DM decays')
                use_sigmaM = np.trapz(self.sigmaM,self.M)/np.trapz(np.ones(len(self.M)),self.M)
                gz = 0.075 + 0.25/(1+self.z)**5
                xi = -2 -0.05/gz * (1-2.4*use_sigmaM**0.05*np.exp(-(4.7/(use_sigmaM*gz))**2))
                rho_0 = 0.048 + 0.77/use_sigmaM
                rho_1 = 4.7*use_sigmaM**1.9*np.exp(-2/use_sigmaM)

                rho_crit = 2.77536627e11*(self.Msunh*self.Mpch**-3)

                rhoM_mean = rho_crit*(self.camb_pars.omegam-self.camb_pars.omeganu)
                      
                logRhoMin = np.log10(rho_crit.value)-5*use_sigmaM
                logRhoMax = np.log10(rho_crit.value)+5*use_sigmaM

                tilde_rho_m = np.logspace(logRhoMin.value,logRhoMax.value,10**4)*rho_crit.unit / rhoM_mean

                Prho_un = tilde_rho_m**xi * np.exp(-(rho_0/tilde_rho_m)**1.1-(tilde_rho_m/rho_1)**0.55) 
                A = np.trapz(Prho_un,tilde_rho_m)**-1
                Prho = A * Prho_un

                PDM = Prho / self.CLT

                PDM_all = interp1d(tilde_rho_m,PDM,bounds_error=False,fill_value='extrapolate')

                rhoL_DM_val = (self.Theta_DM*(u.s**-1)*(rhoM_mean)*(1+self.z)**3*(cu.c.to(u.Mpc/u.s))**2)
               
                if self.subtract_VID_mean:
                    dd = (self.T+self.Tmean).to(u.Msun*u.GHz**-1*u.s**-3*u.sr**-1)*(4.*np.pi*(self.nu*self.H*(1.*u.sr)))/(cu.c.to(u.km/u.s)*rhoL_DM_val)
                else:
                    dd = self.T.to(u.Msun*u.GHz**-1*u.s**-3*u.sr**-1)*(4.*np.pi*(self.nu*self.H*(1.*u.sr)))/(cu.c.to(u.km/u.s)*rhoL_DM_val)

                PDM_vals = PDM_all(dd.value)*PDM.unit / np.trapz(PDM_all(dd.value),dd.value) 
                
                return PDM_vals

            else:
                # use fft's and linearly spaced T points 
                fP1 = fft(self.P1)*self.dT

                fP1 = ((fP1*self.P1.unit).decompose()).value 
                fPT = np.zeros(self.T.size,dtype=complex)

                for ii in range(1,self.Ngal_max+1):
                    fPT += fP1**(ii)*self.PofN[ii].value
                            
                # Errors in fft's leave a small imaginary part, remove for output
                return (ifft(fPT)/self.dT).real
            
        else:
            #  uses brute-force convolutions and logarithmically spaced T points 
            
            P_N = np.zeros([self.Ngal_max,self.T.size])*self.P1.unit
            P_N[0,:] = self.P1
            
            for ii in range(1,self.Ngal_max):
                P_N[ii,:] = conv_parallel(self.T,P_N[ii-1,:], self.T,self.P1,self.T)
            
            PT = np.zeros(self.T.size)

            for ii in range(0,self.Ngal_max):
                PT = PT+P_N[ii,:]*self.PofN[ii+1]
                
            return PT
            
    # PT delta function at T = 0
    @cached_vid_property
    def PT_zero(self):

        return self.PofN[0]
                

    # Converts continuous probability distribution PT(T) to predicted histogram values for bin edges Tedge in Nvov
    def pdf_to_histogram(self,T,PT,Tedge,Nvox,Tmean_sub,PT_zero):

        # P(T) is assumed to be zero outside of the given range of T

        # Bins are assumed to be small compared to the width of P(T)
        
        Tbin = self.binedge_to_binctr(Tedge)
        PTi = interp1d(T,PT,bounds_error=False,fill_value=0)
        
        h = np.zeros(Tbin.size)
        
        for ii in range(0,Tbin.size):
            h[ii] = quad(PTi,Tedge[ii].value,Tedge[ii+1].value)[0]*Nvox
            if Tedge[ii]<=-Tmean_sub<=Tedge[ii+1]:
                h[ii] = h[ii]+PT_zero
        
        return h


    def Add_noise(self,T,Tout,PT,Pzero,sig2):

        conv = np.zeros(len(Tout))*PT.unit

        # convolute signal PDF and gaussian noise

        for i in range(len(Tout)):
            conv[i] = np.trapz(Gauss(Tout[i]-T,0.,sig2)*PT,T)

        conv += Pzero*Gauss(Tout,0.,sig2)
        
        return conv
    

    # VID computation
    def get_VID(self,Tbin,T,Tmean,PT,PT_zero,Nbin,minBi):

        # noise dispersion
        sig2 = self.sigma_N**2
        Nsigma = 5

        # create T array between -(noise dispersion) and max T 
        TT2=np.concatenate((np.linspace(-Nsigma*sig2.value**0.5,0,400)*T.unit, T + Tmean))

        # subtract average T to temperature arrays
        Tout_nomean = np.linspace(-Nsigma*sig2.value**0.5-Tmean.value,T[-1].value,2**8+1)*T.unit

        Tout = self.binedge_to_binctr(Tout_nomean)
        TTnomean = TT2-Tmean

        # convolute noise to PT 
        use_PT=np.concatenate((np.zeros(400)*PT.unit,PT))
        PT_total = self.Add_noise(TTnomean,Tout,use_PT,PT_zero,sig2)

        if not Tbin:

            # define the temperature bins if not input

            # maximum temperature between maximum T and dispersion noise
            Tmax = max(Nsigma*sig2.value**0.5,T[-1].value)

            # final temperature bins boundaries
            Tbin = ulinspace(-Nsigma*sig2**0.5, Tmax*self.T.unit, Nbin)
            
            # compute Ti binning and Bi 
            Ti_all = self.binedge_to_binctr(Tbin)

        else:
            Ti_all = deepcopy(Tbin)
            Tbin = self.binctr_to_binedge(Ti_all)
     
        Bi_all = self.pdf_to_histogram(Tout,PT_total,Tbin,self.Nvox,0.*Tbin.unit,0.)

        if not self.do_Jysr:
            if not Tbin:

                index = 0
                for i in range(len(Ti_all)):
                    if Ti_all[i].value >= 0.:
                        index = i 
                        break

                Ti_vals = Ti_all[index:]
                # return Bi in a single slice of the survey
                Bi_slice = Bi_all[index:] / self.Nvox * self.Npix

                if minBi:

                    # consider only Bi above a certain temperature
                    stop = -1
                    for i in range(len(Bi_slice)):
                        if Ti_vals[i].value > 1. and Bi_slice[i] <= minBi / self.Nvox * self.Npix:
                            stop = i
                            break
                    
                    Ti_vals = Ti_vals[:stop]

                    # return Bi in a single slice of the survey
                    Bi_slice = Bi_slice[:stop] 

            else:
                    Ti_vals = Ti_all
                    Bi_slice = Bi_all / self.Nvox * self.Npix
        else:
            Ti_vals = Ti_all
            Bi_slice = Bi_all / self.Nvox * self.Npix
            for i in range(len(Bi_slice)):
                Bi_slice[i] = max(Bi_slice[i],0.)

        return Ti_vals, Bi_slice