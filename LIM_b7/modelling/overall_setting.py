# SL: last update 01/17/2023

# IMPORT MAIN DIRECTORIES

import astropy.units as u
import astropy.constants as cu

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib as pl
import matplotlib.gridspec as gridspec
from matplotlib import rc
from getdist import plots as gp
import seaborn as sns

from copy import deepcopy
import inspect
import numpy as np
import os
import sys
import types
import importlib

from scipy.interpolate import interp1d, interp2d,RegularGridInterpolator
from scipy.special import sici,erf,legendre,j1
from scipy.stats import poisson
from scipy.integrate import quad
from getdist.gaussian_mixtures import GaussianND

from multiprocessing import Pool as ThreadPool
from functools import partial

import warnings
warnings.filterwarnings("ignore")

import camb 
import LIM_b7.modelling.run_acamb as axions

#PLOT SETTING

kitsune   = '#D9972F'
seiheki   = '#478384'
aonibi      = '#324356'
shion       = '#968ABD'
suoko     = '#B23E52'
chojizome = '#DDB87E'

colors = [suoko,seiheki,kitsune,'#79a43a',shion,'#828282']

colors_low = [suoko, chojizome]
colors_up = [shion, seiheki]
colors_lcdm = [aonibi, aonibi]

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['legend.fontsize'] = 17 
plt.rcParams['figure.figsize'] = (14,7)
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 20 
plt.rcParams['xtick.labelsize'] = 17 
plt.rcParams['ytick.labelsize'] = 17 
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Helvetica'
