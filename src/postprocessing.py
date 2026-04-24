'''========================================================================================
This Python code load data TOMAS

Ali Akherati - CSU - July 2018
========================================================================================'''


# CLEARING iPython
# =========================================================================================
from IPython import get_ipython
get_ipython().magic('reset -sf') 

# IMPORTING LIBRARIES
# =========================================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.patches as mpatch
from matplotlib.font_manager import FontProperties

# Changing font
# =========================================================================================
mpl.rcParams['font.family'] = 'Helvetica'
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# READING
# =========================================================================================
nbins = 36
directory = '/Users/Ali/Desktop/project/SOM_TOMAS/cosov1.2'
filename = 'burn_test_2'
ext = ['gc', 'aemass', 'noconc']

df = pd.read_csv('%s/%s_%s.dat'%(directory, filename, ext[2]),
                 delim_whitespace=True, header=None)
Nk = np.array(df)[:,1:] # data from TOMAS after resolving T0M

df = pd.read_csv('%s/%s_%s.dat'%(directory, filename, ext[1]),
                 delim_whitespace=True, header=None)
mk = np.array(df)[:,1:] # data from TOMAS after resolving T0M
Mk = np.zeros((Nk.shape[0], mk.shape[1], Nk.shape[1]))
for i in range(Nk.shape[0]):
    Mk[i,]

df = pd.read_csv('%s/%s_%s.dat'%(directory, filename, ext[0]),
                 delim_whitespace=True, header=None)
Gc = np.array(df)[:,1:] # data from TOMAS after resolving T0M
'''
directory = '/Users/Ali/Desktop/project/SOM_TOMAS/TOMAS_bian'
filename = 'burn37_160429_chemon2_0.943_2.493_1000000._5843.0_0.157_0.'

df = pd.read_csv('%s/%s.dat'%(directory, filename),
                 delim_whitespace=True, header=None)
Nk_T0M = np.array(df)[:,22:22+nbins] # data from old TOMAS w/ T0M
Gc_T0M = np.array(df)[:,5:22] # data from old TOMAS w/ T0M
Mk_T0M = np.array(df)[:,23:] # data from old TOMAS w/ T0M
'''

