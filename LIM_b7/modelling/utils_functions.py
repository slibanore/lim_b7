# SL: last update 01/17/2023

from .overall_setting import *
from . import halo_functions as hf
from . import Mh_to_Lline as ml

###########################################
# CHECK INPUTS 
###########################################

# Get default parameters
def get_default_params(func):
        
    args = inspect.getargspec(func)
    
    param_names = args.args
    if 'self' in param_names:
        param_names.remove('self')
    
    default_values = args.defaults
    
    default_params = dict(zip(param_names,default_values))

    return default_params

# Check parameter type and units wrt defaults
def check_params(input_params, default_params):
    
    # 1) Check if valid parameter
    for key in input_params.keys():
        if key not in default_params.keys():
            raise AttributeError(key+" is not a valid parameter")
        
        input_value = input_params[key]
        default_value = default_params[key]
        
    # 2) Check if correct type
        if type(input_value)!=type(default_value):
            # Some inputs can have multiple types
            if key=='scatter_seed':
                if type(input_value)==int or type(input_value)==float:
                    pass
                
            elif type(default_value)== u.Quantity:
                raise TypeError("Parameter "+key+
                        " must be an astropy quantity")
            else:
                raise TypeError("Parameter "+key+" must be a "+
                                    str(type(default_value)))
            
    # 2A) If quantity, check correct dimension
        elif (type(default_value)== u.Quantity and not
                 input_value.unit.is_equivalent(default_value.unit)):
            
            # Tmin/Tmax may be in either uK or Jy/sr depending on do_Jysr     
            if key=='Tmin' or key=='Tmax':
                if (input_params['do_Jysr'] and 
                   input_value.unit.is_equivalent(u.Jy/u.sr)):
                    pass
                else:
                    raise TypeError("Parameter "+key+
                                " must have units equivalent to "
                                +str(default_value.unit))
                                
    # 2B) Special requirements for certain parameters
        elif (key=='model_type' and not 
                (input_value=='ML' or input_value=='LF' or input_value=='TOY')):
            # model_type can only be ML or LF
            raise ValueError("model_type must be either 'ML' or 'LF' ot 'TOY' ")

# Check halo mass function existence                    
def check_halo_mass_function_model(hmf_model):

    if not hasattr(hf,hmf_model):
        raise ValueError(hmf_model+
        " mass function not found in halo_functions.py")

# Check bias model existence
def check_bias_model(bias_name):

    if not hasattr(hf,bias_name):
        raise ValueError(bias_name+
        " bias not found in halo_functions.py")

# Check astrophysical models existence
def check_astro_model(model_name):

    if not hasattr(ml,model_name):
            raise ValueError(model_name+
            " model not found in Mh_to_Lline.py")
            


###########################################
# CACHED PROPERTIES
###########################################

# Things are run once and then stored

class cached_property(object):

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        
        instance._update_list.append(self.func.__name__)
        
        res = instance.__dict__[self.func.__name__] = self.func(instance)
        return res


class cached_cosmo_property(object):

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        
        # ADDED THIS CODE TO LIST PROPERTY FOR UPDATING
        instance._update_cosmo_list.append(self.func.__name__)
        
        res = instance.__dict__[self.func.__name__] = self.func(instance)
        return res

class cached_vid_property(object):

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        
        # ADDED THIS CODE TO LIST PROPERTY FOR UPDATING
        instance._update_vid_list.append(self.func.__name__)
        
        res = instance.__dict__[self.func.__name__] = self.func(instance)
        return res

class cached_obs_property(object):

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        
        # ADDED THIS CODE TO LIST PROPERTY FOR UPDATING
        instance._update_obs_list.append(self.func.__name__)
        
        res = instance.__dict__[self.func.__name__] = self.func(instance)
        return res



###########################################
# OTHER FUNCTIONS
###########################################

# Computes logarithmically-spaced numpy array between xmin and xmax with nx points
def ulogspace(xmin,xmax,nx):

    return np.logspace(np.log10(xmin.value),np.log10(xmax.value),nx)*xmin.unit

# Computes linearly-spaced numpy array between xmin and xmax with nx points
def ulinspace(xmin,xmax,nx):

    return np.linspace(xmin.value,xmax.value,nx)*xmin.unit

# log-log interpolation
def log_interp1d(xx, yy, kind='linear',bounds_error=False,fill_value='extrapolate'):
 
    try:
        logx = np.log10(xx.value)
    except:
        logx = np.log10(xx)
    try:
        logy = np.log10(yy.value)
    except:
        logy = np.log10(yy)
    lin_interp = interp1d(logx, logy, kind=kind,bounds_error=bounds_error,fill_value=fill_value)
    
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))

    return log_interp


# lognormal PDF
def lognormal(x,mu,sigma):

    try: 
        return 1/x/sigma/(2.*np.pi)**0.5*np.exp(-(np.log(x.value) - mu)**2/2./sigma**2)
    except:
        return 1/x/sigma/(2.*np.pi)**0.5*np.exp(-(np.log(x) - mu)**2/2./sigma**2)


# gaussian
def Gauss(x,mu,sig2):

    exparg = -0.5*x**2/sig2
    norm = (2.*np.pi*sig2)**0.5

    return np.exp(exparg)/norm


# add two vectors, where k* is the module, mu* is the cosine of the polar angle 
def add_vector(k1, mu1, k2, mu2):

    k1par,k1perp = k1*mu1,k1*(1-mu1**2)
    k2par,k2perp = k2*mu2,k2*(1-mu2**2)
    
    ksum = np.sqrt((k1perp + k2perp)**2 + (k1par + k2par)**2)
    musum = (k1par + k2par)/ksum
    
    return ksum,musum

# convolution brute force between y1(x1) and y2(x2)
def conv_bruteforce(x1,y1,x2,y2,x):
    
    #y1 and y2 are assumed to be zero outside the given range of x1 and x2
    # Interpolate input functions for integration
    y1f = interp1d(x1,y1,bounds_error=False,fill_value=0)
    y2f = interp1d(x2,y2,bounds_error=False,fill_value=0)
    
    xmin = min(np.append(x1,x2))
    Imax = x-xmin
    if Imax<=xmin:
        y = 0.
    else:
        itgrnd = lambda xp: y1f(xp) * y2f(x-xp)
        y = quad(itgrnd, xmin, Imax)[0]

    return y
    
# parallelized brute force convolution
def conv_parallel(x1,y1,x2,y2,x):
    
    # Remove units from quantities
    yunit = y1.unit**2.*x.unit # unit for final convolution
    
    x1 = x1.to(x.unit).value
    x2 = x2.to(x.unit).value
    y2 = y2.to(y1.unit).value
    
    y1 = y1.value
    x = x.value
    
    # Setup parallel pool
    pool = ThreadPool(4)
    
    # Compute convolution in parallel
    fun = partial(conv_bruteforce,x1,y1,x2,y2)
    y = np.array(pool.map(fun,x))
    
    # Close parallel pool
    pool.close()
    pool.join()
    
    # Add units back on 
    return y*yunit


def create_dir(dir_name):
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name,exist_ok=True)

    return 