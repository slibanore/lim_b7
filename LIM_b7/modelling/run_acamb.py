from .overall_setting import *

class run_axion_camb():
    
    # setup cosmological model 
    def __init__(self,cosmo_code = 'camb',\
        cosmo_input_camb=dict(f_NL=0,H0=67.36,cosmomc_theta=None,ombh2=0.02237, omch2=0.12, 
                omk=0.0, neutrino_hierarchy='degenerate', 
                num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
                YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
                TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
                theta_H0_range=[10, 100],w=-1.0, wa=0., cs2=1.0, 
                dark_energy_model='ppf',As=2.1e-09, ns=0.9649, nrun=0, 
                nrunrun=0.0, r=0.0, nt=None, ntrun=0.0, 
                pivot_scalar=0.05, pivot_tensor=0.05,
                parameterization=2,halofit_version='mead'\
                ),\
        transfer_kmax = 100.,\
        f_axion = 0.1,\
        m_axion = 1.e-30,\
        fname= 'axion_CAMB_py/input_params.ini',\
        script= 'axion_CAMB_py/run_script.sh'\
        ):
        
        self.cosmo_code = cosmo_code
        self.cosmo_input_camb = cosmo_input_camb
        
        self.f_axion = f_axion
        self.m_axion = m_axion
        self.transfer_kmax = transfer_kmax
        
        self.fname = fname
        self.script = script
                 
    
    # compute the power spectrum with fuzzy DM
    def run_camb(self):

        #write parameters in an initial file  
        f = open(self.fname, 'r')
        list_of_lines = f.readlines()
        #Replace the lines with new variables
        list_of_lines[3] = 'output_root = ../input\n'
        list_of_lines[34] = 'ombh2 = %f\n'%self.cosmo_input_camb['ombh2']
        list_of_lines[35] = 'omch2 = %f\n'%self.cosmo_input_camb['omch2']
        list_of_lines[37] = 'omk = %f\n'%self.cosmo_input_camb['omk']
        list_of_lines[38] = 'hubble = %f\n'%self.cosmo_input_camb['H0']
        list_of_lines[40] = 'w = %d\n'%self.cosmo_input_camb['w']
        list_of_lines[42] = 'cs2_lam = %d\n'%self.cosmo_input_camb['cs2']
        
        list_of_lines[55] = 'm_ax = %e\n'%self.m_axion
        
        list_of_lines[58] = 'omdah2 = %f\n'%self.cosmo_input_camb['omch2']
        
        list_of_lines[59] = 'axfrac = %f\n'%self.f_axion
        
        list_of_lines[69] = 'temp_cmb = %f\n'%self.cosmo_input_camb['TCMB']
        list_of_lines[89] = 'pivot_scalar = %f\n'%self.cosmo_input_camb['pivot_scalar']
        list_of_lines[90] = 'pivot_tensor = %f\n'%self.cosmo_input_camb['pivot_tensor']
        list_of_lines[91] = 'scalar_amp(1) = %e\n'%self.cosmo_input_camb['As']
        list_of_lines[92] = 'scalar_spectral_index(1) = %f\n'%self.cosmo_input_camb['ns']
        list_of_lines[93] = 'scalar_nrun(1) = %d\n'%self.cosmo_input_camb['nrun']

        list_of_lines[146] = 'transfer_kmax = %f\n'%self.transfer_kmax
        
        f.close()
        f = open(self.fname, 'w')
        f.writelines(list_of_lines)
        f.close()
        
        f = open(self.script,'w')
        f.write('#/bin/bash\n\n')
        f.write('cd axion_CAMB_py/axionCAMB-master\n')
        f.write('./camb ../input_params.ini\n')
        f.close()
        
        cmd = 'chmod +x {:s}'.format(self.script)
        os.system(cmd)
        cmd = './{:s}\t../{:s}'.format(self.script,self.fname)
        os.system(cmd)
        
        return 0
            

            
            
            
            
            
            
            
