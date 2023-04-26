# SL: last update 01/18/2023
from .modelling import *

class Pkline(LineObs):

    def __init__(self,

            do_onehalo=True,
            do_RSD=True,
            sigma_NL=7*u.Mpc,
            FoG_damp='Lorentzian',
            smooth= True,
            do_conv_Wkmin = False,

            **line_kwargs):
                    
        # Initiate LineObs() parameters
        LineObs.__init__(self,**line_kwargs)

        self._update_cosmo_list = self._update_cosmo_list
        
        self._pk_params = locals()
        self._pk_params.pop('self')
        self._pk_params.pop('line_kwargs')
        self._default_pk_params = get_default_params(Pkline.__init__)
        check_params(self._pk_params,self._default_pk_params)
        
        # Set instrument parameters
        for key in self._pk_params:
            setattr(self,key,self._pk_params[key])
        
        # Combine lim_params with pk_params
        self._input_params.update(self._pk_params)
        self._default_params.update(self._default_pk_params)

###########################################
# POWER SPECTRUM QUANTITIES
###########################################

    # Kaiser factor and FoG for RSD
    @cached_property
    def RSD(self):

        if self.do_RSD == True:
            kaiser = (1.+self.f_eff/self.bavg*self.mui_grid**2.)**2.
            
            if self.FoG_damp == 'Lorentzian':
                FoG = (1.+0.5*(self.k_par*self.sigma_NL).decompose()**2.)**-2.
            elif self.FoG_damp == 'Gaussian':
                FoG = np.exp(-((self.k_par*self.sigma_NL)**2.).decompose()) 
            else:
                raise ValueError('Only Lorentzian or Gaussian damping terms for FoG')
                
            return FoG*kaiser
        else:
            return np.ones(self.Pm.shape)

            
    # shot noise amplitude for target line at frequency nuObs
    @cached_property
    def Pshot(self):

        if self.developer == 'decayDM':
            return np.zeros(self.Pm.shape)*(self.Pk_twohalo).unit
        else:
            return self.CLT**2*self.L2mean

    # two-halo term in power spectrum      
    @cached_property
    def Pk_twohalo(self):
    
        # if not one halo term equals Tmean^2*bavg^2*Pm 
        if self.do_onehalo:
            if self.developer == 'decayDM':
                wt = self.CLT
            else:
                Mass_Dep = self.LofM*self.dndM
                itgrnd = np.tile(Mass_Dep,(self.k.size,1)).T*self.ft_NFW*self.bofM
                # Special case for SFR(M) scatter in Tony Li model
                wt = self.CLT*np.trapz(itgrnd,self.M,axis=0)*self.fduty

            if self.model_name=='TonyLi':
                alpha = self.model_par['alpha']
                sig_SFR = self.model_par['sig_SFR']
                wt = wt*np.exp((alpha**-2-alpha**-1)
                    *sig_SFR**2*np.log(10)**2/2.)
        else:
            if self.developer == 'decayDM':
                wt = self.CLT
            else:
                wt = self.Tmean*self.bavg
        
        return wt**2*self.Pm
        

    # one halo term in power spectrum 
    @cached_property
    def Pk_onehalo(self):

        if self.do_onehalo:

            if self.developer == 'decayDM':
                
                return np.zeros(self.Pm.shape)*self.Pk_twohalo.unit
            else:
                Mass_Dep = self.LofM**2.*self.dndM
                itgrnd = np.tile(Mass_Dep,(self.nk,1)).T*self.ft_NFW**2.
                #add effect for the scatter in LCO
                itgrnd = itgrnd*np.exp(self.sigma_scatter**2*np.log(10)**2)
                            
                # Tony Li model- scatter does not preserve LCO
                if self.model_name=='TonyLi':
                    alpha = self.model_par['alpha']
                    sig_SFR = self.model_par['sig_SFR']
                    itgrnd = itgrnd*np.exp((2.*alpha**-2-alpha**-1)
                                        *sig_SFR**2*np.log(10)**2)
                wt = np.trapz(itgrnd,self.M,axis=0)*self.fduty
                
                return np.tile(self.CLT**2.*wt,(self.nmu,1))

        else:
            return np.zeros(self.Pm.shape)*self.Pshot.unit
    
    
    # clustering power spectrum of target line without shot noise
    @cached_property    
    def Pk_clust(self):

        return (self.Pk_twohalo+self.Pk_onehalo)*self.RSD
        

    # shot noise power spectrum
    @cached_property    
    def Pk_shot(self):

        return self.Pshot*np.ones(self.Pm.shape)
        

    # Full line power spectrum including both clustering and shot noise as function of k and mu
    @cached_property    
    def Pk(self):
         
        if self.smooth:
            if self.do_conv_Wkmin:
            # convolve with the survey mask window assuming a cylindrical volume

                Pkres = self.Wk*(self.Pk_clust+self.Pk_shot)
                
                #Get the vector to integrate over
                qe = np.logspace(-4,2,self.nk+1)
                q = 0.5*(qe[:-1]+qe[1:])*u.Mpc**-1
                muq = self.mu
                
                qi_grid,muqi_grid = np.meshgrid(q,muq)
                q_par = qi_grid*muqi_grid
                q_perp = qi_grid*np.sqrt(1-muqi_grid**2)
                
                #get the window to convolve with
                L_perp=np.sqrt(self.Sfield/np.pi)
                Wpar = 2*np.sin((q_par*self.Lfield/2).value)/q_par
                Wperp = 2*np.pi*L_perp*j1(q_perp*L_perp)/q_perp
                Wconv = Wpar*Wperp
                
                #Do the convolution
                Pconv = np.zeros(Pkres.shape)*Pkres.unit*self.Vfield.unit
                Pkres_interp = RegularGridInterpolator((self.k.value,self.mu),Pkres.T.value, bounds_error=False, fill_value=0)
                for ik in range(self.nk):
                    for imu in range(self.nmu):
                        #Get the unconvolved power spectrum in the sum of vectors
                        qsum_grid,musum_grid = add_vector(self.k[ik],self.mu[imu],qi_grid,-muqi_grid)
                        Pconv[imu,ik] = np.trapz(np.trapz(qi_grid**2*Pkres_interp((qsum_grid.value,musum_grid.value))*Pkres.unit*np.abs(Wconv**2)/(2*np.pi)**2,muq,axis=0),q)

                return Pconv/self.Vfield
            else:
                return self.Wk*(self.Pk_clust+self.Pk_shot)
        else:
            return self.Pk_clust+self.Pk_shot
          
            
    # power spectrum monopole   
    @cached_property
    def Pk_0(self):

        return 0.5*np.trapz(self.Pk,self.mu,axis=0)
        
    # power spectrum quadrupole
    @cached_property
    def Pk_2(self):

        L2 = legendre(2)

        return 2.5*np.trapz(self.Pk*L2(self.mui_grid),self.mu,axis=0)
        
    # power spectrum hexadecapole        
    @cached_property
    def Pk_4(self):

        L4 = legendre(4)

        return 4.5*np.trapz(self.Pk*L4(self.mui_grid),self.mu,axis=0)
        
    # power spectrum multipole L
    def Pk_l(self,l):

        if l == 0:
            return self.Pk_0
        elif l == 2:
            return self.Pk_2
        elif l == 4:
            return self.Pk_4
        else:
            Ll = legendre(l)
            return (2.*l+1.)/2.*np.trapz(self.Pk*Ll(self.mui_grid),self.mu,axis=0)


    # instrumental noise power spectrum amplitude
    @cached_obs_property    
    def Pnoise(self):

        if self.do_Jysr:
            return (self.sigma_N / (self.Nvox/self.combined_pixels**2))**2*self.Vvox 
        else:
            return self.sigma_N**2*self.Vvox

###########################################
# COMPUTE POWER SPECTRUM S and N COVARIANCE
###########################################
                
    # Number of modes between k and k+dk.        
    @cached_obs_property
    def Nmodes(self):

        # Multiply by dmu/2 to get the number of modes between k and k+dk and mu and mu+dmu

        return self.ki_grid**2*self.dk*self.Vfield*self.Nfield/4./np.pi**2.
                
    
    # 00 term of the covariance matrix from CV
    @cached_obs_property
    def covmat_CV_00(self):

        return 0.5*np.trapz(self.Pk**2/self.Nmodes,self.mu,axis=0)


    # 02 term of the covariance matrix from CV       (equal to 20)
    @cached_obs_property
    def covmat_CV_02(self):

        L2 = legendre(2)(self.mui_grid)

        return 5./2.*np.trapz(self.Pk**2*L2/self.Nmodes,self.mu,axis=0)
        
        
    # 04 term of the covariance matrix from CV         (equal to 20)
    @cached_obs_property
    def covmat_CV_04(self):

        L4 = legendre(4)(self.mui_grid)

        return 9./2.*np.trapz(self.Pk**2*L4/self.Nmodes,self.mu,axis=0)
        
    # 22 term of the covariance matrix from CV
    @cached_obs_property
    def covmat_CV_22(self):

        L2 = legendre(2)(self.mui_grid)
        return 25./2.*np.trapz(self.Pk**2*L2*L2/self.Nmodes,self.mu,axis=0)


    # 24 term of the covariance matrix from CV
    #    (equal to 42)
    @cached_obs_property
    def covmat_CV_24(self):

        L2 = legendre(2)(self.mui_grid)
        L4 = legendre(4)(self.mui_grid)

        return 45./2.*np.trapz(self.Pk**2*L2*L4/self.Nmodes,self.mu,axis=0)


    # 44 term of the covariance matrix from CV     
    @cached_obs_property
    def covmat_CV_44(self):

        L4 = legendre(4)(self.mui_grid)

        return 81./2.*np.trapz(self.Pk**2*L4*L4/self.Nmodes,self.mu,axis=0)
        

    # l1l2 term of the covariance matrix from CV
    def covmat_CV_l1l2(self,l1,l2):

        if l1 == 0 and l2 == 0:
            return self.covmat_CV_00
        elif l1 == 0 and l2 == 2:
            return self.covmat_CV_02
        elif l1 == 0 and l2 == 4:
            return self.covmat_CV_04
        elif l1 == 2 and l2 == 2:
            return self.covmat_CV_22
        elif l1 == 2 and l2 == 4:
            return self.covmat_CV_24
        elif l1 == 4 and l2 == 4:
            return self.covmat_CV_44
        else:
            Ll1 = legendre(l1)(self.mui_grid)
            Ll2 = legendre(l2)(self.mui_grid)
            return (2.*l1+1.)*(2.*l2+1.)/2.*np.trapz(self.Pk**2*Ll1*Ll2/self.Nmodes,self.mu,axis=0)


    # 00 term of the covariance matrix from instrumental noise          
    @cached_obs_property
    def covmat_N_00(self):

        return 1./2.*np.trapz(self.Pnoise**2./(self.Nmodes),self.mu,axis=0)
        

    # 02 term of the covariance matrix from instrumental noise (equal to 20)        
    @cached_obs_property
    def covmat_N_02(self):

        L2 = legendre(2)(self.mui_grid)

        return 5./2.*np.trapz(self.Pnoise**2.*L2/(self.Nmodes),self.mu,axis=0)


    # 04 term of the covariance matrix from instrumental noise (equal to 40)
    @cached_obs_property
    def covmat_N_04(self):

        L4 = legendre(4)(self.mui_grid)

        return 9./2.*np.trapz(self.Pnoise**2.*L4/(self.Nmodes),self.mu,axis=0)
        

    # 22 term of the covariance matrix from instrumental noise        
    @cached_obs_property
    def covmat_N_22(self):

        L2 = legendre(2)(self.mui_grid)

        return 25./2.*np.trapz(self.Pnoise**2.*L2*L2/(self.
        Nmodes),self.mu,axis=0)


    # 24 term of the covariance matrix from instrumental noise (equal to 42)
    @cached_obs_property
    def covmat_N_24(self):

        L2 = legendre(2)(self.mui_grid)
        L4 = legendre(4)(self.mui_grid)

        return 45./2.*np.trapz(self.Pnoise**2.*L2*L4/(self.Nmodes),self.mu,axis=0)
        

    # 44 term of the covariance matrix from instrumental noise        
    @cached_obs_property
    def covmat_N_44(self):

        L4 = legendre(4)(self.mui_grid)

        return 81./2.*np.trapz(self.Pnoise**2.*L4*L4/(self.Nmodes),self.mu,axis=0)
        

    # l1l2 term of the covariance matrix from instrumental noise        
    def covmat_N_l1l2(self,l1,l2):

        if l1 == 0 and l2 == 0:
            return self.covmat_N_00
        elif l1 == 0 and l2 == 2:
            return self.covmat_N_02
        elif l1 == 0 and l2 == 4:
            return self.covmat_N_04
        elif l1 == 2 and l2 == 2:
            return self.covmat_N_22
        elif l1 == 2 and l2 == 4:
            return self.covmat_N_24
        elif l1 == 4 and l2 == 4:
            return self.covmat_N_44
        else:
            Ll1 = legendre(l1)(self.mui_grid)
            Ll2 = legendre(l2)(self.mui_grid)
            return (2.*l1+1.)*(2.*l2+1.)*np.trapz(self.Pnoise**2.*Ll1*Ll2/(self.Nmodes),self.mu,axis=0)
        
        
    # Total error at k and mu
    @cached_obs_property
    def sk(self):

        # sample variance errror
        sk_CV = self.Pk/np.sqrt(self.Nmodes*self.dmu[0])

        # Error at k and mu due to instrumental noise
        sk_N = self.Pnoise/(np.sqrt(self.Nmodes*self.dmu[0]/2.))

        return sk_CV+sk_N
        
        
    # 00 term of the total covariance matrix
    @cached_obs_property
    def covmat_00(self):

        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 0.5*np.trapz(integrand**2,self.mu,axis=0)
        

    # 02 term of the total covariance matrix         
    @cached_obs_property
    def covmat_02(self):

        L2 = legendre(2)(self.mui_grid)
        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 5./2.*np.trapz(integrand**2*L2,self.mu,axis=0)
        

    # 04 term of the total covariance matrix        
    @cached_obs_property
    def covmat_04(self):

        L4 = legendre(4)(self.mui_grid)
        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 9./2.*np.trapz(integrand**2*L4,self.mu,axis=0)
        
        
    # 22 term of the total covariance matrix
    @cached_obs_property
    def covmat_22(self):

        L2 = legendre(2)(self.mui_grid)
        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 25./2.*np.trapz(integrand**2*L2*L2,self.mu,axis=0)
        

    # 24 term of the total covariance matrix        
    @cached_obs_property
    def covmat_24(self):

        L2 = legendre(2)(self.mui_grid)
        L4 = legendre(4)(self.mui_grid)
        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 45./2.*np.trapz(integrand**2*L2*L4,self.mu,axis=0)
        
        
    # 44 term of the total covariance matrix
    @cached_obs_property
    def covmat_44(self):

        L4 = legendre(4)(self.mui_grid)
        integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

        return 81./2.*np.trapz(integrand**2*L4*L4,self.mu,axis=0)


    # l1l2 term of the total covariance matrix
    def covmat_l1l2(self,l1,l2):

        if l1 == 0 and l2 == 0:
            return self.covmat_00
        elif l1 == 0 and l2 == 2:
            return self.covmat_02
        elif l1 == 0 and l2 == 4:
            return self.covmat_04
        elif l1 == 2 and l2 == 2:
            return self.covmat_22
        elif l1 == 2 and l2 == 4:
            return self.covmat_24
        elif l1 == 4 and l2 == 4:
            return self.covmat_44
        else:
            l1 = legendre(l1)(self.mui_grid)
            l2 = legendre(l2)(self.mui_grid)
            integrand = (self.Pk+self.Pnoise)/self.Nmodes**0.5

            return (2.*l1+1.)*(2.*l2+1.)*np.trapz(integrand**2*l1*l2,self.mu,axis=0)


    # Signal to noise ratio for given model and experiment
    @cached_obs_property
    def SNR(self):

        SNR_k = (self.Pk**2/self.sk**2).decompose()
        ind = np.logical_and(self.k>=self.kmin_field,self.k<=self.kmax_field)

        return np.sqrt(SNR_k[ind].sum())


    # Signal to noise ratio in the monopole for given model and experiment
    @cached_obs_property
    def SNR_0(self):

        SNR_k = (self.Pk_0**2/self.covmat_00).decompose()
        ind = np.logical_and(self.k>=self.kmin_field,self.k<=self.kmax_field)

        return np.sqrt(SNR_k[ind].sum())
        

    # Signal to noise ratio in the quadrupole for given model and experiment        
    @cached_obs_property
    def SNR_2(self):

        SNR_k = (self.Pk_2**2/self.covmat_22).decompose()
        ind = np.logical_and(self.k>=self.kmin_field,self.k<=self.kmax_field)

        return np.sqrt(SNR_k[ind].sum())
        
        
    # Signal to noise ratio in the hexadecapole for given model and experiment
    @cached_obs_property
    def SNR_4(self):

        SNR_k = (self.Pk_4**2/self.covmat_44).decompose()
        ind = np.logical_and(self.k>=self.kmin_field,self.k<=self.kmax_field)

        return np.sqrt(SNR_k[ind].sum())
        

    # Signal to noise ratio in the monopole, quadrupole and hexadecapole for given model and experiment     
    @cached_obs_property
    def SNR_multipoles(self):

        ind = np.where(np.logical_and(self.k>=self.kmin_field,
                                      self.k<=self.kmax_field))[0]
        Nkseen = len(ind)
        Pkvec = np.zeros(Nkseen*3)
        covmat = np.zeros((Nkseen*3,Nkseen*3))
        
        Pkvec[:Nkseen] = self.Pk_0[ind]
        Pkvec[Nkseen:Nkseen*2] = self.Pk_2[ind]
        Pkvec[Nkseen*2:Nkseen*3] = self.Pk_4[ind]
        
        covmat[:Nkseen,:Nkseen] = np.diag(self.covmat_00[ind])
        covmat[:Nkseen,Nkseen:Nkseen*2] = np.diag(self.covmat_02[ind])
        covmat[:Nkseen,Nkseen*2:Nkseen*3] = np.diag(self.covmat_04[ind])
        covmat[Nkseen:Nkseen*2,:Nkseen] = np.diag(self.covmat_02[ind])
        covmat[Nkseen:Nkseen*2,Nkseen:Nkseen*2] = np.diag(self.covmat_22[ind])
        covmat[Nkseen:Nkseen*2,Nkseen*2:Nkseen*3] = np.diag(self.covmat_24[ind])
        covmat[Nkseen*2:Nkseen*3,:Nkseen] = np.diag(self.covmat_04[ind])
        covmat[Nkseen*2:Nkseen*3,Nkseen:Nkseen*2] = np.diag(self.covmat_24[ind])
        covmat[Nkseen*2:Nkseen*3,Nkseen*2:Nkseen*3] = np.diag(self.covmat_44[ind])
        
        return np.sqrt(np.dot(Pkvec,np.dot(np.linalg.inv(covmat),Pkvec)))
        
        
    # covariance matrix for a given number of multipoles
    def get_covmat(self,Nmul):

        # starting always from the monopole and without skipping any

        if Nmul > 3:
            raise ValueError('Not implemented yet!\
            Implement covmat_66 and expand this function')
            
        covmat = np.zeros((self.nk*Nmul,self.nk*Nmul))
        covmat[:self.nk,:self.nk] = np.diag(self.covmat_00)
        
        if Nmul > 1:
            covmat[:self.nk,self.nk:self.nk*2] = np.diag(self.covmat_02)
            covmat[self.nk:self.nk*2,:self.nk] = np.diag(self.covmat_02)
            covmat[self.nk:self.nk*2,self.nk:self.nk*2] = np.diag(self.covmat_22)
            covmat[:self.nk,self.nk:self.nk*2] = np.diag(self.covmat_02)
        if Nmul > 2:
            covmat[:self.nk,self.nk*2:self.nk*3] = np.diag(self.covmat_04)
            covmat[self.nk:self.nk*2,self.nk*2:self.nk*3] = np.diag(self.covmat_24)
            covmat[self.nk*2:self.nk*3,:self.nk] = np.diag(self.covmat_04)
            covmat[self.nk*2:self.nk*3,self.nk:self.nk*2] = np.diag(self.covmat_24)
            covmat[self.nk*2:self.nk*3,self.nk*2:self.nk*3] = np.diag(self.covmat_44)

        return covmat
        
