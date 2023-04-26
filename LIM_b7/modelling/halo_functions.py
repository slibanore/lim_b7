# SL: last update 01/17/2023

from .overall_setting import *

###########################################
# HALO MASS FUNCTION
###########################################

def Tinker(self,Mvec,rhoM,z):
    '''
    Tinker et al 2008 halo mass function for delta=200
    '''

    A = self.hmf_pars['A_tinker']
    a = self.hmf_pars['a_tinker']
    b = self.hmf_pars['b_tinker']
    c = self.hmf_pars['c_tinker']

    fs = A*((self.sigmaM/b)**(-a)+1.0)*np.exp(-c/self.sigmaM**2)

    dndM = fs*(rhoM/Mvec)*(-self.dsigmaM_dM.to(self.Msunh**-1)/self.sigmaM)

    return dndM


def NG_Riotto(self,Mvec,rhoM,z):

    fNL = self.fNL 

    dndM_G = Tinker(self,Mvec,rhoM,z)

    delta_c = 1.686

    S = 1.8e-4 * fNL / (self.sigmaM**0.838 * self.Dgrowth(z)**0.162)

    dSdsigma = - 0.838 * S / self.sigmaM

    dc = 0.949 * delta_c
    T = np.sqrt(1-dc*S/3)

    C_Ng = (dc**2 /(6*T) * dSdsigma + T) * np.exp( S*dc**3 / (6*self.sigmaM**2))

    dndM = dndM_G * C_Ng

    return dndM

###########################################
# HALO BIAS
###########################################

def Tinker10(self,dc,nu):

    if len(self.bias_par.keys()) == 0:
        y = np.log10(200.)
        A = 1. + 0.24*y*np.exp(-(4./y)**4.)
        a = 0.44*y - 0.88
        B = 0.183
        b = 1.5
        C = 0.019 + 0.107*y + 0.19*np.exp(-(4./y)**4.)
        c = 2.4
    else:
        y = self.bias_par['y']
        B = self.bias_par['B']
        b = self.bias_par['b']
        c = self.bias_par['c']
        A = 1. + 0.24*y*np.exp(-(4./y)**4.)
        C = 0.019 + 0.107*y + 0.19*np.exp(-(4./y)**4.)
        a = 0.44*y - 0.88
    
    return 1.- A*nu**a/(nu**a+dc**a) + B*nu**b + C*nu**c
