# SL: last update 01/17/2023

from .cosmo_astro_prescriptions import *

# initialize line model + experimental observation
class LineObs(LineModel):

    def __init__(self, 
                 Tsys_NEFD = 40*u.K,
                 Nfeeds = 19,
                 beam_FWHM = 4.1*u.arcmin,
                 Delta_nu = 8*u.GHz,
                 dnu = 15.6*u.MHz,
                 tobs = 6000*u.hr, 
                 Omega_field = 2.25*u.deg**2,
                 Nfield = 1,
                 N_FG_par = 1, #Foreground multiplicative number to kmin_par in window
                 N_FG_perp = 1, #Foreground multiplicative number to kmin_perp in window
                 do_FG_wedge = False, #Apply foreground wedge removal
                 a_FG = 0.*u.Mpc**-1, #Constat superhorizon buffer for foregrounds
                 b_FG = 0., #Foreground parameter accounting for antenna chromaticity
                 input_sigmaN = -1.*u.uK,
                 combined_pixels = -1,
                 mag_limit = 22.4,
                 nsigma = 5,
                 px_side = 5.45*u.arcsec,
                 **line_kwargs):
                    
        # Initiate LineModel() parameters
        LineModel.__init__(self,**line_kwargs)
        
        self._update_cosmo_list = self._update_cosmo_list
        
        self._obs_params = locals()
        self._obs_params.pop('self')
        self._obs_params.pop('line_kwargs')
        self._default_obs_params = get_default_params(LineObs.__init__)
        check_params(self._obs_params,self._default_obs_params)
        
        # Set instrument parameters
        for key in self._obs_params:
            setattr(self,key,self._obs_params[key])
        
        # Combine lim_params with obs_params
        self._input_params.update(self._obs_params)
        self._default_params.update(self._default_obs_params)

###########################################
# DETECTOR RELATED QUANTITIES
###########################################

    # frequency channels
    @cached_obs_property
    def Nch(self):

        return np.round((self.Delta_nu/self.dnu).decompose())


    # beam width (1sigma of Gaussian beam)        
    @cached_obs_property
    def beam_width(self):

        return self.beam_FWHM*0.4247
  
  
    # number of pixels 
    @cached_obs_property
    def Npix(self):

        # pixel size = beam FWHM on a side
        theta_side = np.sqrt(self.Omega_field)
        Nside = np.floor((theta_side/self.beam_FWHM).decompose())

        if self.combined_pixels > 0:
            return Nside**2/self.combined_pixels**2
        else:
            return Nside**2


    # number of voxels
    @cached_obs_property
    def Nvox(self):

        return self.Npix*self.Nch


    # covered sky fraction        
    @cached_obs_property
    def fsky(self):

        return (self.Omega_field/(4*np.pi*u.rad**2)).decompose()


###########################################
# DETECTOR NOISE
###########################################

    # time spent observing each pixel with a single detector
    @cached_obs_property
    def tpix(self):

        return self.tobs/self.Npix
    
    # instrumental noise per voxel  
    @cached_obs_property
    def sigma_N(self):

        if self.do_Jysr:
            # equivalent to sigma_pix in this case
            if self.input_sigmaN.value > 0.: 
                  return self.input_sigmaN
            else:
                #return ((self.Tsys_NEFD/self.beam_FWHM**2).to(u.Jy*u.s**(1./2)/u.sr))
                # SL - ULTRASAT
                fnu = 10**((self.mag_limit - 8.90)/(-2.5))*u.Jy
                if self.combined_pixels > 0:
                    fnu_px = fnu/(((self.px_side*self.combined_pixels)**2).to(u.sr))
                else:
                    fnu_px = fnu/(((self.px_side)**2).to(u.sr))
                sigma_N = (fnu_px/self.nsigma) 
                return (sigma_N).to(u.Jy/u.sr)
        else:
            if self.input_sigmaN.value > 0.: 
                  return self.input_sigmaN      
            else:
                  return ((self.Tsys_NEFD/np.sqrt(self.Nfeeds*self.dnu*self.tpix)).to(u.uK))
    
    # sky-averaged brightness temperature at nuObs from target line
    @cached_property
    def Tmean(self):

        if self.developer == 'decayDM':
            return self.CLT
        else:    
            return self.CLT*self.Lmean
        
    # sky-averaged luminosity density at nuObs from target line
    @cached_property
    def Lmean(self):

        itgrnd = self.LofM*self.dndM
        Lbar = np.trapz(itgrnd,self.M)*self.fduty   
        # Tony Li model- scatter does not preserve LCO
        if self.model_name=='TonyLi':
            alpha = self.model_par['alpha']
            sig_SFR = self.model_par['sig_SFR']
            Lbar = Lbar*np.exp((alpha**-2-alpha**-1)
                    *sig_SFR**2*np.log(10)**2/2.)

        return Lbar


    # sky-averaged squared luminosity density at nuObs from target line
    @cached_property
    def L2mean(self):

        itgrnd = self.LofM**2*self.dndM
        L2bar = np.trapz(itgrnd,self.M)*self.fduty
        # Add L vs. M scatter
        L2bar = L2bar*np.exp(self.sigma_scatter**2*np.log(10)**2)
        # Tony Li model- scatter does not preserve LCO
        if self.model_name=='TonyLi':
            alpha = self.model_par['alpha']
            sig_SFR = self.model_par['sig_SFR']
            L2bar = L2bar*np.exp((2.*alpha**-2-alpha**-1)
                                *sig_SFR**2*np.log(10)**2)

        return L2bar
        
   
###########################################
# PROJECTED QUANTITIES IN THE SKY 
###########################################

    # comoving distance to central redshift of field
    @cached_obs_property
    def r0(self):

        return self.cosmo.comoving_radial_distance(self.z)*u.Mpc
    

    # area of single field in the sky in Mpc**2
    @cached_obs_property
    def Sfield(self):

        return (self.r0**2*(self.Omega_field/(1.*u.rad**2))).to(u.Mpc**2)


    # length side of a pixel in Mpc        
    @cached_obs_property
    def Lpix_side(self):

        pix_area = self.Sfield/self.Npix
        return pix_area**0.5
        
    # depth of a single field
    @cached_obs_property
    def Lfield(self):

        z_min = (self.nu/(self.nuObs+self.Delta_nu/2.)-1).value
        z_max = (self.nu/(self.nuObs-self.Delta_nu/2.)-1).value

        dr_los = (self.cosmo.comoving_radial_distance(z_max)-
                    self.cosmo.comoving_radial_distance(z_min))

        return dr_los*u.Mpc


    # comoving volume of a single field
    @cached_obs_property    
    def Vfield(self):

        return self.Sfield*self.Lfield
    

    # comoving volume of a single voxel    
    @cached_obs_property            
    def Vvox(self):

        return self.Vfield/self.Nvox


###########################################
# MODES CUT OFFS 
###########################################
    
    #  High-resolution cutoff for line-of-sight modes
    @cached_obs_property
    def sigma_par(self):

        return (cu.c*self.dnu*(1+self.z)/(self.H*self.nuObs)).to(u.Mpc)
    
    
    #  High-resolution cutoff for transverse modes
    @cached_obs_property
    def sigma_perp(self):

        return (self.r0*(self.beam_width/(1*u.rad))).to(u.Mpc)
            
 
    # Minimum k accessible in a single field
    @cached_obs_property
    def kmin_field(self):
        
        # set by the maximum side length
        # Minimum k in transverse direction    
        kmin_sky = 2*np.pi/self.Sfield**0.5

        # Minimum k in line of sight direction              
        kmin_los = 2*np.pi/self.Lfield
    
        return min([kmin_los,kmin_sky])

    
    # Maximum k accessible in a single field
    @cached_obs_property
    def kmax_field(self):
        
        # set by the best resolution
        # Maximum k in transverse direction    
        kmax_sky = 2.*np.pi/self.sigma_perp

        # Maximum k in line of sight direction              
        kmax_los = 2.*np.pi/self.sigma_par

        return max([kmax_los,kmax_sky])
        

    # Resolution cutoff in power spectrum         
    @cached_obs_property
    def Wkmax(self):

        # Resolution cutoff in power spectrum in the los direction
        Wkmax_par = np.exp(-((self.k_par*self.sigma_par)**2).decompose())
        # Resolution cutoff in power spectrum in the transverse direction
        Wkmax_perp = np.exp(-((self.k_perp*self.sigma_perp)**2).decompose())

        return Wkmax_par*Wkmax_perp
        

    # Precision cutoff in power spectrum 
    # due to volume 
    @cached_obs_property
    def Wkmin(self):

        # Only relevant if do_conv_Wkmin = False
        # Precision cutoff in power spectrum 
        # due to volume observed in los direction 
        kmin_los = 2*np.pi/self.Lfield
        Wkmin_par = 1.-np.exp(-((self.k_par/(self.N_FG_par*kmin_los))**2).decompose())

        # Precision cutoff in power spectrum 
        # due to volume observed in transverse direction
        kmin_sky = 2*np.pi/self.Sfield**0.5
        Wkmin_perp = 1.-np.exp(-((self.k_perp/(self.N_FG_perp*kmin_sky))**2).decompose())

        return Wkmin_par*Wkmin_perp
        
    
    # Apply foreground wedge removal
    @cached_obs_property
    def Wk_FGwedge(self):

        W = np.ones(self.ki_grid.shape)
        if self.do_FG_wedge:
            kpar_min_wedge = self.a_FG.to(self.k.unit) + self.b_FG*np.abs(self.k_perp)
            ind = np.where(np.abs(self.k_par)<kpar_min_wedge)
            W[ind] = 0
            return W
        else:
            return W
        

    # Resolution cutoff in power spectrum   
    @cached_obs_property
    def Wk(self):

        if not self.do_conv_Wkmin:
            return self.Wkmin*self.Wkmax*self.Wk_FGwedge
        else:
            return self.Wkmax*self.Wk_FGwedge

        
