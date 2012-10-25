'''
script to construct and save simulated photometric observations
for a set of stellar models

MUST BE RUN WITH STSDAS + PYTHON ENVIRONMENT
'''


# IMPORTS
import pysynphot as ps
import numpy as np
import matplotlib.pyplot as plt

# pull in all of the recommended models, copied from 
# http://www.stsci.edu/hst/HST_overview/documents/synphot/AppA_Catalogs4.html#48115
mods = np.loadtxt('recommended_models.txt', skiprows=2, dtype='S')
mod_names = mods[:,0]
mod_temps = np.array([row[-1].split('_')[1].split('[')[0] for row in mods]).astype('float')
mod_gravs = np.array([row[-1].split('_')[1].split('g')[1].strip(']') for row in mods]).astype('float') * .1

# get zeropoint corrections from http://www.stsci.edu/hst/HST_overview/documents/synphot/c034.html#335184
# for [y, B, V, R]
vegamag_corrections = [.038, .036, .026, .038]

# define the filters desired
filt_names = ['sdss,u', 'sdss,g', 'sdss,r', 'sdss,i', 'sdss,z', 'stromgren,y', 'B', 'V', 'R', 'J', 'H', 'K']
filts = []
for fn in filt_names:
    filts.append( ps.ObsBandpass(fn) )

# 2MASS has it's own internal mag system, so we need those zeropoints.
# from http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
#  all in Jy
tm_zps = [1594., 1024., 666.7]

# go through models and build all of the observations
big_list = []
fails = []
for i,model in enumerate(mod_names):
    spec = ps.Icat('ck04models', mod_temps[i], 0., mod_gravs[i]) #using Z=0. for all models
    little_list = []
    try:
        for j,filt in enumerate(filts):
            ob = ps.Observation(spec, filt)
            # get abmags from sloan filters, and vega mags from y and BVR (incorporating correction), and 2MASS mags from JHK
            if j < 5:
                mag = ob.effstim('abmag')
            elif 5 <= j < 9:
                correction = vegamag_corrections[ j-5 ]
                mag = ob.effstim('vegamag') + correction
            elif 9 <= j:
                zeropoint = tm_zps[ j-9 ]
                flux = ob.effstim('jy')
                mag = -2.5 * np.log10(flux/zeropoint)
            little_list.append(mag)
    except:
        fails.append(i)
    if little_list:
        big_list.append(little_list)

# save to a numpy savefile, for use later
np.save( open('all_models.npy', 'w'), big_list)
            
print 'failed on:'
for val in fails:
    print mod_names[val]
    
    