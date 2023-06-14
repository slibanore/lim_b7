# SL: last update 01/17/2023

from .utils_functions import *
# TalAdi: used for reading PMF .dat files
import csv
from random import randint
import time

# initialize the class containing all the relevant quantities to model the line emission
class LineModel(object):

    def __init__(self,    

    ###########################################
    # GENERAL QUANTITIES
    ###########################################

        # COSMOLOGICAL MODEL
        cosmo_input_camb = dict(
        f_NL=0, H0=67.67, cosmomc_theta=None,
        ombh2=0.0224, omch2=0.1193, omk=0.0, neutrino_hierarchy='degenerate', 
        num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
        YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
        TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
        theta_H0_range=[10, 100], w=-1.0, wa=0., cs2=1.0, 
        dark_energy_model='ppf',As=2.105e-09, ns=0.967, nrun=0., nrunrun=0.0, r=0.0, nt=None, ntrun=0.0, 
        pivot_scalar=0.05, pivot_tensor=0.05,
        parameterization=2,halofit_version='mead'),

        # HALO PROPERTIES
        Mmin = 1e9*u.Msun,
        Mmax = 1e15*u.Msun,
        nM = 1000,
        hmf_model = 'Tinker',
        hmf_pars = dict(\
        A_tinker = 0.186, a_tinker = 1.47,
        b_tinker = 2.57, c_tinker = 1.19),
        bias_model = 'Tinker10',
        bias_par = {},

        # ASTROPHYSICAL MODEL FOR LINE EMISSION
        nu = 115.271*u.GHz,
        nuObs = 30*u.GHz,
        model_name = 'TonyLi',
        model_par = {\
        'alpha':1.37, 'beta':-1.74, \
		'dMF':1.,'sig_SFR':0.3, \
		'SFR_file':os.getcwd() + '/LIM_b7/modelling/SFR_tables/sfr_release.dat'},
        sigma_scatter=0.3,
        Lmin = 1e-5*u.Lsun,
        Lmax = 1e10*u.Lsun,
        nL = 5000,
        dndL_Lcut=50.*u.Lsun,

        # POWER SPECTRUM - RELATED QUANTITIES
        kmin = 1e-2*u.Mpc**-1,
        kmax = 100.*u.Mpc**-1,
        nk = 500,
        fduty=1.,
        nonlinear=False,
        nmu=1000,

        # VID - RELATED QUANTITIES
        do_Jysr=False,
        do_fast_VID=True,
        Tmin_VID=1e-5*u.uK,
        Tmax_VID=120.*u.uK,

        # set which kind of analysis you want to use
        developer = 'CDS', 
        
    ###########################################
    # PROPERTIES RELATED WITH SPECIFIC ANALYSIS
    ###########################################

        # CDS : arXiv:2208.01658
        # LCDM deviations and related scales
        CDS_alpha_1=1., 
        CDS_alpha_2=1., 
        CDS_alpha_3=1., 
        CDS_alpha_4=1., 
        k_12 = 0.5, 
        k_23 = 2., 
        k_34 = 10.,

        # axions
        # mass and damping scale for fuzzy DM 
        m_axion = 1e-24, 
        f_axion = 0.1, 

        # PMF
        n_B = -2.9,         # PMF spectral index
        sigma_B_0 = 0.1,    # PMF variance in nGauss
        smooth_scale = 1.,   # PMF smooth scale in Mpc

    ):

###########################################
# FUNCTIONS NEEDED FOR COSMOLOGY
###########################################

# 1) Set up the required stuff

        # Check type and units wrt default names and values
        self._lim_params = locals()
        self._lim_params.pop('self')
        self._default_lim_params = get_default_params(LineModel.__init__)
        check_params(self._lim_params,self._default_lim_params)
        
        # Set parameters
        for key in self._lim_params:
            setattr(self,key,self._lim_params[key])

        # Create parameter lists 
        self._input_params = {} 
        self._default_params = {}
        self._input_params.update(self._lim_params)
        self._default_params.update(self._default_lim_params)
        
        # Create cached properties list
        self._update_list = []
        self._update_cosmo_list = []
        self._update_obs_list = []
        self._update_vid_list = []
        
        # Check model_name 
        check_halo_mass_function_model(self.hmf_model)
        check_bias_model(self.bias_model)
        check_astro_model(self.model_name)
        
# 2) Set and compute cosmology through CAMB

        # set cosmo parameters
        self.cosmo_input_camb = self._default_params['cosmo_input_camb']
        for key in cosmo_input_camb:
            self.cosmo_input_camb[key] = cosmo_input_camb[key]

        self.camb_pars = camb.set_params(\
            H0=self.cosmo_input_camb['H0'], cosmomc_theta=self.cosmo_input_camb['cosmomc_theta'],
            ombh2=self.cosmo_input_camb['ombh2'], omch2=self.cosmo_input_camb['omch2'], omk=self.cosmo_input_camb['omk'],
            neutrino_hierarchy=self.cosmo_input_camb['neutrino_hierarchy'], 
            num_massive_neutrinos=self.cosmo_input_camb['num_massive_neutrinos'],
            mnu=self.cosmo_input_camb['mnu'], nnu=self.cosmo_input_camb['nnu'], YHe=self.cosmo_input_camb['YHe'], 
            meffsterile=self.cosmo_input_camb['meffsterile'], 
            standard_neutrino_neff=self.cosmo_input_camb['standard_neutrino_neff'], 
            TCMB=self.cosmo_input_camb['TCMB'], tau=self.cosmo_input_camb['tau'], 
            deltazrei=self.cosmo_input_camb['deltazrei'], 
            bbn_predictor=self.cosmo_input_camb['bbn_predictor'], 
            theta_H0_range=self.cosmo_input_camb['theta_H0_range'],
            w=self.cosmo_input_camb['w'], cs2=self.cosmo_input_camb['cs2'], 
            dark_energy_model=self.cosmo_input_camb['dark_energy_model'],
            As=self.cosmo_input_camb['As'], ns=self.cosmo_input_camb['ns'], 
            nrun=self.cosmo_input_camb['nrun'], nrunrun=self.cosmo_input_camb['nrunrun'], 
            r=self.cosmo_input_camb['r'], nt=self.cosmo_input_camb['nt'], ntrun=self.cosmo_input_camb['ntrun'], 
            pivot_scalar=self.cosmo_input_camb['pivot_scalar'], 
            pivot_tensor=self.cosmo_input_camb['pivot_tensor'],
            parameterization=self.cosmo_input_camb['parameterization'],
            halofit_version=self.cosmo_input_camb['halofit_version'])
                
        self.camb_pars.WantTransfer=True    
        self.camb_pars.Transfer.accurate_massive_neutrinos = True
    
    #  z array to interpolate cosmological quantities
    @cached_cosmo_property
    def zcosmo(self):

        zmax = 15.  
        Nz = 150
        if self.z > zmax:
            raise ValueError('Required z_obs outside interpolation region. Increase zmax or change nuObs')

        return np.linspace(0.,zmax,Nz)


    # cosmological evolution
    @cached_cosmo_property
    def cosmo(self):

        self.camb_pars.set_matter_power(redshifts=list(self.zcosmo))

        return camb.get_results(self.camb_pars)


    # matter transfer function (k in Mpc-1)
    @cached_property
    def transfer_m(self):

        zz = self.zcosmo[::-1] #camb sortes earlier first

        # Find above and below values for z in zcosmo
        iz_down = np.where(zz - self.z < 0)[0][0]
        iz_up = iz_down - 1
        dz = zz[iz_up] - zz[iz_down]
        
        # Get transfer
        T = self.cosmo.get_matter_transfer_data()
        kvec = (T.transfer_z('k/h',-1)*self.Mpch**-1).to(u.Mpc**-1)
        Tk_up = T.transfer_z('delta_tot',iz_up)
        Tk_down = T.transfer_z('delta_tot',iz_down)

        #interpolate in z (linear)
        Tz = Tk_down*(1.-(self.z-zz[iz_down])/dz) + Tk_up*(self.z-zz[iz_down])/dz

        #interpolate in k (linear)
        return log_interp1d(kvec,Tz)


    # CDM + baryons transfer function
    # change wrt transfer_m because of massive nu
    @cached_property
    def transfer_cb(self):

        zz = self.zcosmo[::-1] #camb sortes earlier first

        # Find above and below values for z in zcosmo
        iz_down = np.where(zz - self.z < 0)[0][0]
        iz_up = iz_down - 1
        dz = zz[iz_up] - zz[iz_down]
        
        #Get transfer
        T = self.cosmo.get_matter_transfer_data()
        kvec = (T.transfer_z('k/h',-1)*self.Mpch**-1).to(u.Mpc**-1)
        Tk_up = T.transfer_z('delta_nonu',iz_up)
        Tk_down = T.transfer_z('delta_nonu',iz_down)

        #interpolate in z (linear)
        Tz = Tk_down*(1.-(self.z-zz[iz_down])/dz) + Tk_up*(self.z-zz[iz_down])/dz

        #interpolate in k (linear)
        return log_interp1d(kvec,Tz)


    # Non Gaussianity
    @cached_cosmo_property
    def f_NL(self):

        return self.cosmo_input_camb['f_NL']


    # Alcock - Packynsky: Da/rs, H*rs, DV/rs
    @cached_cosmo_property
    def Alcock_Packynski_params(self):

        BAO_pars = self.cosmo.get_BAO(self.zcosmo[1:],self.camb_pars)
        #This is rs/DV, H, DA, F_AP
        rs = self.cosmo.get_derived_params()['rdrag']
        DA = BAO_pars[:,2]
        DV = rs/BAO_pars[:,0]
        Hz = BAO_pars[:,1]
            
        DA_over_rs_int = interp1d(self.zcosmo[1:],DA/rs,kind='cubic', bounds_error=False, fill_value='extrapolate')

        DV_over_rs_int = interp1d(self.zcosmo[1:],DV/rs,kind='cubic', bounds_error=False,fill_value='extrapolate')

        H_times_rs_int = interp1d(self.zcosmo[1:],Hz*rs,kind='cubic', bounds_error=False, fill_value='extrapolate')
    
        return DA_over_rs_int, H_times_rs_int,DV_over_rs_int

# TalAdi: PMF induced matter power spectrum function that runs WolframScript
    @cached_property
    def pk_pmf(self):
        '''
        Description:
        Get interpolated PMF induced matter power spectrum P_mM(k,z) = Pi(k) * M2(z).
        Returns Pi(k), M2(z) as interpolated functions.

        DIMENSIONS:
        [Pi(k)] = Mpc^3 / s^4 ; [M2(z)] = s^4   ->  [P_mM(k,z)] = Mpc^3
        '''
        # Use LCDM, by adding 0 instead
        if self.n_B <= -3. or self.sigma_B_0 == 0.:
            zero_list = np.zeros(2)
            zero_interp = interp1d(zero_list, zero_list,
                                bounds_error=False, fill_value=0)
            return zero_interp, zero_interp
        
        # ==========================
        # Conversions and Constants
        # ==========================
        nGauss_to_kgOs = 5.290144928591178e-32 * u.kg / u.s   # [kg/s]
        # ---------------------------------------------------------------
        c = 2.99792458e8 * u.m / u.s  # Speed of light in [m/s]
        hbar = 6.62607015e-34 / 2 / np.pi * u.Joule * u.s   # Planck Constant /2pi in [J*sec]
        # ---------------------------------------------------------------
        sigma_T = 6.6524616e-29 * (u.m)**2    # Thomson cross-section in [m^2]
        m_p = 1.6735575e-27 * u.kg  # Proton mass in [kg]
        G_newton = 6.6743e-11 * u.Joule * u.m / (u.kg)**2  # Newton gravitational constant in [J*m/kg^2]

        # ========================
        # Cosmological parameters
        # ========================
        h = self.cosmo_input_camb['H0'] / 100
        Omega_b = self.cosmo_input_camb['ombh2'] / h**2
        Omega_m = self.cosmo_input_camb['omch2'] / h**2 + Omega_b
        YHe = 0.2454
        Omega_g = 5.4e-05
        Omega_Lambda = 1- Omega_m - Omega_g
        H0 = 1e5 * h * u.m / u.s / u.Mpc
        rho_crit = (3 * H0**2 / (8 * np.pi * G_newton)).decompose()  # [kg/m^3]

        # ===============
        # PMF parameters
        # ===============
        sigma_B_0 = self.sigma_B_0  # [nGauss]
        sigma_B_0_MKS = self.sigma_B_0 * nGauss_to_kgOs # [kg/s]
        smooth_scale = self.smooth_scale * u.Mpc    # [Mpc]
        n_B = self.n_B
        # ---------------
        k_min = 2e-4    # [1/Mpc]
        N_k = 50
        
        # =================================================
        # Function for calculating Alfven cutoff and z_rec
        # =================================================
        def get_kcut_zrec():
            # Get ionization fraction (x_e = n_e / n_H)
            bknd_data = camb.get_background(self.camb_pars)
            z_sample = 10**np.linspace(-1,8,3000)
            x_e = bknd_data.get_background_redshift_evolution(z_sample, ['x_e'])['x_e']
            xe_intrp = interp1d(z_sample, x_e, fill_value='extrapolate')
            # Find z_rec
            ne_nb = lambda z: xe_intrp(z)*(1 - YHe) / (1 - 0.74811 * YHe)
            z_rec = 1300
            while ne_nb(z_rec) > 0.1:
                z_rec -= 0.01
            # The integrand 
            Ih_integrand = lambda z: (1+z)**2/(
                        (Omega_b * (1+z)**3 + 4/3 * Omega_g * (1+z)**4) *
                        (Omega_m * (1+z)**3 + 4/3 * Omega_g * (1+z)**4 + Omega_Lambda)**0.5 * xe_intrp(z))
            # Integrating over redshift
            I_h, Ih_err = quad(Ih_integrand, z_rec, np.inf)
            # Calculate kutoff scale and inverting into cutoff k
            cut_scale = smooth_scale * (((2*np.pi)**2 * (0.59*m_p) * sigma_B_0_MKS**2 / smooth_scale**2 / sigma_T / (
                1-YHe) / H0 / rho_crit**2 / Omega_b * ((c.value**2*u.Joule/u.kg) / hbar / c**2)).decompose() * I_h)**(1/(5+n_B))    # [Mpc]
            
            k_cut = 2*np.pi / cut_scale.value   # [1/Mpc]

            return k_cut,z_rec
        
        # ============================
        # Generate PMF power spectrum
        # ============================
        wd = os.getcwd()
        cd = wd + '/PMF_cache/'
        
        input_params = [h, Omega_b, Omega_m, sigma_B_0, smooth_scale.value, n_B, k_min, N_k]
        pmf_file_name = 'PmM'
        for val in input_params:
            nval = float(f'{val:.{6}g}')
            sval = format(nval, '.10f').rstrip('0').rstrip('.')
            pmf_file_name += '_' + sval
        pmf_file_name += '.dat'
        
        run_count = 0
        while not os.path.isfile(cd+pmf_file_name):
            # Routine to get k_cut and z_rec
            k_cut, z_rec = get_kcut_zrec()
            # write input file for wolframscript
            id_num = str(randint(1111,9999))
            ini_name = 'input_' + id_num + '.ini'
            while os.path.isfile(cd+ini_name):
                id_num = str(randint(1111,9999))
                ini_name = 'input_' + id_num + '.ini'
            os.makedirs(os.path.dirname(cd+ini_name), exist_ok=True)
            f = open(cd+ini_name,'w')
            f.write(str([h, Omega_b, Omega_m, sigma_B_0, smooth_scale.value, n_B, k_min, N_k, k_cut, z_rec])[1:-1])
            f.close()
            # run wolframscript
            cmd = str(wd) + '/LIM_b7/generate_pmf.wls' + ' ' + id_num
            print('Calculating PMF...')
            os.system(cmd)
            # Remove input file
            os.remove(cd+ini_name)
            run_count += 1
            time.sleep(1)
            if run_count > 3:
                print('PMF script failed more than 3 times!')
                break
        if run_count == 0:
            print('Using cached PMF files...')

        # read PMF induced matter power spectrum
        k_list, pi_list, z_list, m2_list = [], [], [], []
        with open(cd+pmf_file_name, newline='') as f:
            spamreader = csv.reader(f)
            for row in spamreader:
                if row[0] != 'Nan':
                    k_list.append(float(row[0]))
                    pi_list.append(float(row[1]))
                if row[2] != 'Nan':
                    z_list.append(float(row[2]))
                    m2_list.append(float(row[3]))
        k_list = np.array(k_list)
        pi_list = np.array(pi_list)
        z_list = np.array(z_list)
        m2_list = np.array(m2_list)


        Pi_k = interp1d(k_list, pi_list, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
        M2 = interp1d(z_list, m2_list, kind='cubic',
                                bounds_error=False, fill_value='extrapolate')
        return Pi_k,M2

    # compute matter power spectrum depending on the model
    # note that the growth factor is always the LCDM one
    def PKint(self,z,k):

        if self.developer == 'CDS':
        # interpolate LCDM + plug deviations
            zmax = self.zcosmo[-1]
            nz_step=64
            if self.camb_pars.num_nu_massive != 0:
                var = 8
            else:
                var = 7
            PK = camb.get_matter_power_interpolator(self.camb_pars, zmin=0, zmax=zmax, nz_step=nz_step, 
            zs=None, kmax=100, nonlinear=self.nonlinear,
            var1=var, var2=var, hubble_units=False, 
            k_hunit=False, return_z_k=False, k_per_logint=None, log_interp=False,  extrap_kmax=True)

            Pk0 = PK.P(0.,k)
    
            alpha = lambda k: \
                self.CDS_alpha_1 if k < self.k_12 else \
                self.CDS_alpha_2 if self.k_12 <= k < self.k_23 \
                else self.CDS_alpha_3 if self.k_23 <= k < self.k_34 \
                else self.CDS_alpha_4

            #output =  np.genfromtxt('pk_fid.txt')[:2000,1]
            output = Pk0 * np.vectorize(alpha)(k) * self.Dgrowth(z)**2

        elif self.developer == 'axions':
        # fuzzy dark matter - compute local power spectrum and project back to z

            camb_axions = axions.run_axion_camb(m_axion = self.m_axion, f_axion = self.f_axion,\
            transfer_kmax = (self.kmax.value / \
            (self.cosmo_input_camb['H0'] / 100.)) )
            camb_axions.run_camb()
            
            filename = open('axion_CAMB_py/input_matterpower.dat','r')

            temp = np.loadtxt(filename)
            imported_k = temp[:,0]  
            imported_pk = temp[:,1] 
            Pk0 = interp1d(imported_k,imported_pk)
            
            pk_z = Pk0(k / (self.cosmo_input_camb['H0'] / 100.) ) * self.Dgrowth(z)**2 

            # convert k input to k / h used in interpol
            # and convert (Mpc/h)^3 to Mpc^3
            output = pk_z / (self.cosmo_input_camb['H0'] / 100.)**3 
                
        elif self.developer == 'PMF':
            # TalAdi: use kmax from dictionary
            # (should set the value in fiducial_pars accordingly)
            
            kmax = self.kmax.value
            zmax = self.zcosmo[-1]
            nz_step = 64
            if self.camb_pars.num_nu_massive != 0:
                var = 8
            else:
                var = 7
            # calculate LCDM matter power spectrum
            PK = camb.get_matter_power_interpolator(self.camb_pars, zmin=0, zmax=zmax, nz_step=nz_step, 
            zs=None, kmax=kmax, nonlinear=self.nonlinear,
            var1=var, var2=var, hubble_units=False, 
            k_hunit=False, return_z_k=False, k_per_logint=None, log_interp=False,  extrap_kmax=True)
            Pk0 = PK.P(0.,k)
            output = Pk0 *self.Dgrowth(z)**2
            # add PMF induced mattewr power spectrum (if condition is satisfied)
            try:
                Pi_k, M2 = self.pk_pmf(self)
            except TypeError:
                # print('Using cached PMF...')
                Pi_k, M2 = self.pk_pmf
            output += Pi_k(k) * M2(z)   # [Mpc^3]
            # if self.n_B > -3.:
            #     try:
            #         Pi_k, M2 = self.pk_pmf(self)
            #     except TypeError:
            #         # print('Using cached PMF...')
            #         Pi_k, M2 = self.pk_pmf
            #     output += Pi_k(k) * M2(z)   # [Mpc^3]
            # else:
            #     print('PMF induced matter power spectrum cannot be computed for n_B <= -3.')
            #     print('Resuming without PMF...')

        else:
                print('--- Developer not recognized ---')
                output = -1

        return output


    # matter power spectrum depending on the model used
    @cached_property
    def Pm(self):

        #return self.PKint(self.z,self.ki_grid.value)*u.Mpc**3
        return self.PKint(self.z,self.ki_grid)*u.Mpc**3


    # interpolate effective f depending on k,z + include tiling to multply by mu 
    @cached_property
    def f_eff(self):

        fs8lin = self.cosmo.get_fsigma8()
        s8lin = self.cosmo.get_sigma8()
        fz = interp1d(self.zcosmo[::-1],fs8lin/s8lin,kind='cubic')(self.z)
        #Apply correction if massive nu
        if self.camb_pars.num_nu_massive != 0:
            factor = self.transfer_m(self.k.value)/self.transfer_cb(self.k.value)
        else:
            factor = self.transfer_m(self.k.value)/self.transfer_m(self.k.value)

        return np.tile(fz*factor,(self.nmu,1))


    # LCDM growth factor 
    @cached_cosmo_property
    def Dgrowth(self):

        s8lin = self.cosmo.get_sigma8()

        return interp1d(self.zcosmo[::-1],s8lin/s8lin[-1],kind='cubic', bounds_error=False,fill_value='extrapolate')

    # 1/h units
    @cached_cosmo_property
    def hubble(self):

        return self.camb_pars.H0/100.
    
    @cached_cosmo_property
    def Mpch(self):

        return u.Mpc / self.hubble        
        
    @cached_cosmo_property
    def Msunh(self):

        return u.Msun / self.hubble


###########################################
# PROPERTIES AT TARGET REDSHIFT 
###########################################

    # emission redshift of target line
    @cached_property
    def z(self):

        return (self.nu/self.nuObs-1.).value
    
    # Hubble parameter at target redshift
    @cached_property
    def H(self):

        return self.cosmo.hubble_parameter(self.z)*(u.km/u.Mpc/u.s)

    # Convert luminosity density to brightness temperature
    @cached_property
    def CLT(self):

        if self.do_Jysr:
            x = cu.c/(4.*np.pi*self.nu*self.H*(1.*u.sr))
            return x.to(u.Jy*u.Mpc**3/(u.Lsun*u.sr))
        else:
            x = cu.c**3*(1+self.z)**2/(8*np.pi*cu.k_B*self.nu**3*self.H)
            return x.to(u.uK*u.Mpc**3/u.Lsun)
    

###########################################
# ARRAYS DEFINITION
###########################################

    # mass array 
    @cached_property
    def M(self):

        return ulogspace(self.Mmin,self.Mmax,self.nM)
    
    # luminosity array
    @cached_property
    def L(self):

        return ulogspace(self.Lmin,self.Lmax,self.nL)
        
    # wavenumber bin edges        
    @cached_property
    def k_edge(self):

        return ulogspace(self.kmin,self.kmax,self.nk+1)
    
    # k array    
    @cached_property
    def k(self):

        #k_val = np.genfromtxt('pk_fid.txt')[:2000,0]
        
        Nedge = self.k_edge.size
        
        #return k_val 
        return (self.k_edge[0:Nedge-1]+self.k_edge[1:Nedge])/2.
    
    # width of k bins
    @cached_property
    def dk(self):

        return np.diff(self.k_edge)
        
    # cos theta bin edges        
    @cached_property
    def mu_edge(self):

        return np.linspace(-1,1,self.nmu+1)
        
    # cos theta array        
    @cached_property
    def mu(self):

        Nedge = self.mu_edge.size
        return (self.mu_edge[0:Nedge-1]+self.mu_edge[1:Nedge])/2.
        
    # width of cos theta bins  
    @cached_property
    def dmu(self):

        return np.diff(self.mu_edge)
        
    # grid of k for anisotropic
    @cached_property
    def ki_grid(self):

        return np.meshgrid(self.k,self.mu)[0]
        
    # grid of cos theta for anisotropic        
    @cached_property
    def mui_grid(self):

        return np.meshgrid(self.k,self.mu)[1]
        
    # grid of k parallel        
    @cached_property
    def k_par(self):
        '''
        Grid of k_parallel
        '''
        return self.ki_grid*self.mui_grid
        
    # grid of k perpendicular
    @cached_property
    def k_perp(self):

        return self.ki_grid*np.sqrt(1.-self.mui_grid**2.)


        
###########################################
# HALO PROPERTIES
###########################################

    # Halo mass function
    @cached_property
    def dndM(self):

        Mvec = self.M.to(self.Msunh)
        rho_crit = 2.77536627e11*(self.Msunh*self.Mpch**-3).to(self.Msunh*self.Mpch**-3) #h^2 Msun/Mpc^3

        #Use Omega_m or Omega_cdm+Omega_b wheter mnu = 0 or > 0
        rhoM = rho_crit*(self.camb_pars.omegam-self.camb_pars.omeganu)
        
        mf = getattr(hf,self.hmf_model)(self,Mvec,rhoM,self.z)
        
        return mf.to(u.Mpc**-3*u.Msun**-1)
        

    # Mass (or cdm+b) variance at target redshift
    @cached_property
    def sigmaM(self):

        #Get R(M) and P(k)
        rho_crit = 2.77536627e11*(self.Msunh*self.Mpch**-3).to(u.Msun*u.Mpc**-3) #Msun/Mpc^3

        k = np.logspace(np.log10(self.kmin.value),np.log10(self.kmax.value),self.nk)*u.Mpc**-1

        #Use rho_m or rho_cb depending on mnu
        Pk = self.PKint(self.z,k.value)*u.Mpc**3
        rhoM = rho_crit*(self.camb_pars.omegam-self.camb_pars.omeganu)

        R = (3.0*self.M/(4.0*np.pi*rhoM))**(1.0/3.0)
                
        # Get the window of a configuration space tophat
        kvec = (np.tile(k,[R.size,1]).T)
        Pk = np.tile(Pk,[R.size,1]).T
        R = np.tile(R,[k.size,1])
        x = ((kvec*R).decompose()).value
        W = 3.0*(np.sin(x) - x*np.cos(x))/(x)**3 
        
        # When Pk is too sharp, use Gaussian in k space 
        # for now only enters for axions
        #!!! sharpness defibed through drop wrt to LCDM pk at k = 100)
        if self.developer == 'PMF': 
            # print('\n-------\n SHARP K CUT OFF WINDOW \n-------\n')
            W = np.exp(-(x/2.)**2. / 2.)
               
        #Compute sigma(M)
        integrnd = Pk*W**2*kvec**2/(2.*np.pi**2)
        sigma = np.sqrt(np.trapz(integrnd,kvec[:,0],axis=0)) 

        return sigma


    # derivative of sigma(M) with respect to M 
    @cached_property
    def dsigmaM_dM(self):

        sigmaint = log_interp1d(self.M,self.sigmaM,fill_value='extrapolate')
        Mminus = self.M/1.001
        Mplus =  self.M*1.001
        sigma_minus = sigmaint(Mminus.value)
        sigma_plus = sigmaint(Mplus.value)

        return (sigma_plus-sigma_minus)/(Mplus-Mminus)

    # Halo bias as a function of mass (and scale, if fNL != 0)
    @cached_property
    def bofM(self):

        # nonlinear overdensity
        dc = 1.686
        nu = dc/self.sigmaM
        
        bias = np.tile(getattr(hf,self.bias_model)(self,dc,nu),(self.k.size,1)).T
        
        Delta_b = 0.
        if self.f_NL != 0:
            #get the transfer function, depending on whether mnu = 0 or mnu > 0

            if self.camb_pars.num_nu_massive != 0:
                Tk = self.transfer_cb(self.k.value)
            else:
                Tk = self.transfer_m(self.k.value)
            Om0 = self.camb_pars.omegam

            #Compute non-Gaussian correction Delta_b
            factor = self.f_NL*dc* 3.*Om0*(100.*self.hubble*(u.km/u.s/u.Mpc))**2./   \
            (cu.c.to(u.km/u.s)**2.*self.k**2*(Tk/np.max(Tk))*self.Dgrowth(self.z))

            Delta_b = (bias-1.)*np.tile(factor,(self.nM,1))
            
        return bias + Delta_b
        

    # Fourier transform of NFW halo profile
    @cached_property
    def ft_NFW(self):

        #smaller sampling of M
        Mvec = ulogspace(self.Mmin,self.Mmax,256).value
        #fit parameters
        kappa = 0.42
        a0 = 2.37
        a1 = 1.74
        b0 = 3.39
        b1 = 1.82
        ca = 0.2
        #Compute the effective slope of the growth factor
        dz = self.z*0.001
        alpha_eff = -(np.log(self.Dgrowth(self.z+dz))-np.log(self.Dgrowth(self.z-dz)))/ \
                    (np.log(1.+self.z+dz)-np.log(1.+self.z-dz))
        #Compute the effective slope to the power spectrum (as function of M)
        fun_int = -2.*3.*self.M/self.sigmaM*self.dsigmaM_dM-3.
        neff = interp1d(np.log10(self.M.value),fun_int,fill_value='extrapolate',kind='linear')(np.log10(kappa*Mvec))
        #Quantities for c
        A = a0*(1.+a1*(neff+3))
        B = b0*(1.+b1*(neff+3))
        C = 1.-ca*(1.-alpha_eff)
        nu = 1.686/log_interp1d(self.M.value,self.sigmaM)(Mvec)
        arg = A/nu*(1.+nu**2/B)
        #Compute G(x), with x = r/r_s, and evaluate c
        x = np.logspace(-3,3,256)
        g = np.log(1+x)-x/(1.+x)

        c = np.zeros(len(Mvec))
        for iM in range(len(Mvec)):
            G = x/g**((5.+neff[iM])/6.)
            invG = log_interp1d(G,x,fill_value='extrapolate',kind='linear')
            c[iM] = C*invG(arg[iM])
            
        c_NFW = log_interp1d(Mvec,c,fill_value='extrapolate',kind='cubic')(self.M.value)

        #Radii of the SO collapsed (assuming 200*rho_crit)
        Delta = 200.
        rho_crit = 2.77536627e11*(self.Msunh*self.Mpch**-3).to(u.Msun*u.Mpc**-3) #Msun/Mpc^3
        R_NFW = (3.*self.M/(4.*np.pi*Delta*rho_crit))**(1./3.)
        #get characteristic radius
        r_s = np.tile(R_NFW/c_NFW,(self.nk,1)).T
        #concentration to multiply with ki
        c = np.tile(c_NFW,(self.nk,1)).T
        gc = np.log(1+c)-c/(1.+c)
        #argument: k*rs
        ki = np.tile(self.k,(self.nM,1))
        x = ((ki*r_s).decompose()).value        
        si_x, ci_x = sici(x)
        si_cx, ci_cx = sici((1.+c)*x)
        u_km = (np.cos(x)*(ci_cx - ci_x) +
                  np.sin(x)*(si_cx - si_x) - np.sin(c*x)/((1.+c)*x))

        return u_km/gc


###########################################
# LINE RELATED PROPERTIES
###########################################
       
    # Luminosity - halo mass relation
    @cached_property
    def LofM(self):
    
        L = getattr(ml,self.model_name)(self,self.M,self.model_par,self.z)
    
        return L

    # Luminosity function 
    @cached_property
    def dndL(self):

        #compute LF from the conditional LF
        if self.Lmin > self.LofM[self.LofM.value>0][0]:
            print(self.LofM[self.LofM.value>0][0],self.Lmin)
            print('Warning! reduce Lmin to cover all luminosities of the model')
        if self.Lmax < np.max(self.LofM):
            print('Warning! increase Lmax to cover all luminosities of the model')

        # Tony Li model- scatter does not preserve LCO
        sigma = max(self.sigma_scatter,0.05)

        if self.model_name=='TonyLi':
            alpha = self.model_par['alpha']
            sig_SFR = self.model_par['sig_SFR']
            # sigma and sig_SFR are totally uncorrelated
            sigma = (sigma**2 + sig_SFR**2/alpha**2)**0.5

        sigma_base_e = sigma*2.302585

        #assume a lognormal PDF for the CLF with minimum logscatter of 0.05
        CLF_of_M = np.zeros((self.nM,self.nL))*self.dndM.unit*self.L.unit**-1

        for iM in range(self.nM):
            CLF_of_M[iM,:] = lognormal(self.L,np.log(self.LofM[iM].value)-0.5*sigma_base_e**2.,sigma_base_e)*self.dndM[iM]

        LF = np.zeros(self.nL)*self.L.unit**-1*self.dndM.unit*self.M.unit
        for iL in range(self.nL):
            LF[iL] = np.trapz(CLF_of_M[:,iL],self.M)
        #Add a cut off at low luminosities to ease computations. Default 0*u.Lsun
        LF *= np.exp(-self.dndL_Lcut/self.L)

        # Luminosity function for UV
        # integrated in all the z range
        # Libanore, Kovetz (2023) - ULTRASAT white paper
        #!!! needs to be checked         
        if self.model_name == 'UV_continuum_SFR_integrate':
            use_z = np.linspace(0.,3.5,100)
            L_z = np.zeros((len(self.L),len(use_z)))
            for i in range(len(use_z)):
                #assume a lognormal PDF for the CLF with minimum logscatter of 0.05
                CLF_of_M = np.zeros((self.nM,self.nL))*self.dndM.unit*self.L.unit**-1
    
                for iM in range(self.nM):
                    CLF_of_M[iM,:] = lognormal(self.L,np.log(self.LofM[iM].value)-0.5*sigma_base_e**2.,sigma_base_e)*self.dndM[iM]

                LF_intg = np.zeros(self.nL)*self.L.unit**-1*self.dndM.unit*self.M.unit
                for iL in range(self.nL):
                    LF_intg[iL] = np.trapz(CLF_of_M[:,iL],self.M)

                #Add a cut off at low luminosities to ease computations
                LF_intg *= np.exp(-self.dndL_Lcut/self.L)
                for l in range(len(LF_intg)):
                    L_z[l,i] = LF_intg[l].value

            LF = np.trapz(L_z,use_z,axis=1)	
            LF *= u.Lsun**-1*u.Mpc**-3

        return LF


    # Average luminosity-weighted bias for the given line
    @cached_property
    def bavg(self):

        # Integrands for mass-averaging
        factor = np.tile(self.LofM*self.dndM,(self.nk,1)).T
        itgrnd1 = self.bofM*factor
        itgrnd2 = factor
        
        b_line = np.trapz(itgrnd1,self.M,axis=0) / np.trapz(itgrnd2,self.M,axis=0)
        
        return b_line 


    # Mean number density of galaxies
    @cached_property
    def nbar(self):

        nbar = np.trapz(self.dndM,self.M)

        return nbar
