'''
A synthetic photometry constructor.
This script produces a file called all_models.npy,
 which is used to produce catalogs of SEDs by the 
 functions in get_SEDs.py

NOTE:
Requires the STSCI pysynphot environment, with the CBDB
 database present, as well as a subfolder named 
 transmission_functions/ containing the relative spectral
 response functions for 2Mass and the PS1 y band - see below.
'''


################################################
import numpy as np
import matplotlib.pyplot as plt
import pysynphot as ps
import os, re


################################################
# PASSBANDS
pb_dir = 'transmission_functions/'

# SDSS passbands
u = ps.ObsBandpass('sdss,u')
g = ps.ObsBandpass('sdss,g')
r = ps.ObsBandpass('sdss,r')
i = ps.ObsBandpass('sdss,i')
z = ps.ObsBandpass('sdss,z')

# DECam passbands
#u_d = ps.FileBandpass(pb_dir+'u_DECAM_RSR.dat')
#g_d = ps.FileBandpass(pb_dir+'g_DECAM_RSR.dat')
#r_d = ps.FileBandpass(pb_dir+'r_DECAM_RSR.dat')
#i_d = ps.FileBandpass(pb_dir+'i_DECAM_RSR.dat')
#z_d = ps.FileBandpass(pb_dir+'z_DECAM_RSR.dat')
y_d = ps.FileBandpass(pb_dir+'y_DECAM_RSR.dat')

# USNOB1 Johnson passbands
B = ps.ObsBandpass('B')
R = ps.ObsBandpass('R')
# Additional Johnson passbands
V = ps.ObsBandpass('V')
I = ps.ObsBandpass('I')

# 2Mass passbands
J = ps.FileBandpass(pb_dir+'J_2Mass_RSR.dat')
H = ps.FileBandpass(pb_dir+'H_2Mass_RSR.dat')
K = ps.FileBandpass(pb_dir+'K_2Mass_RSR.dat')



passbands = [u,g,r,i,z, y_d, B,V,R,I, J,H,K]

################################################
# SPECTRA

# Using the Pickles stellar atlas, described in detail at:
#  http://www.stsci.edu/hst/observatory/cdbs/pickles_atlas.html
directory = '/usr/stsci/cdbs/grid/pickles/dat_uvk/'

fnames = os.listdir( directory )
spectra = [fn for fn in fnames if '.fits' in fn]
spectra.remove('pickles_uk.fits')

save_spectra = False

################################################
# PERFORM SYNTHETIC PHOTOMETRY
print 'here we go...'
count = 0
all_colors = []

# first get the effective wavelengths for all passbands
single = [0.]
for bp in passbands:
    single.append( bp.avgwave() )
all_colors.append(single)   

#now go through all spectra and perform photometry
# Note: fills in 0K for temperature, only matters in plots, but if
#  I keep Pickles should make this into spectral type (have types for all models,
#  but don't have temperatures for all models)
for fn in spectra:
    p_index = float(re.findall( '\d+', fn )[0])
    
    #impose limits on what models to use
    if (45 < p_index < 60) or (105 < p_index): continue
    
    print fn
    spec = ps.FileSpectrum( directory+fn )
    if save_spectra:
        np.save( 'model_spectra/pickles/{}.npy'.format(fn.split('.')[0]),
                  np.vstack( (spec.wave, spec.flux)) )
    single_phot = []
    for bp in passbands:
        ob = ps.Observation(spec, bp)
        
        # Use ABmags for ugrizy, but Vegamags for BRJHK
        if bp in [B,R, J,H,K, V,I]:
            mag = ob.effstim('vegamag')
        else:
            mag = ob.effstim('abmag')
        single_phot.append( mag )

    # Save the model with the 0th index as index
    single_colors = [p_index] + single_phot
    all_colors.append( single_colors )
    
    # ((obsolete)) Save the model relative to K=0, with the 0th index as the temperature (C&K) or number (Pickles)
    #single_colors = [p_index]+[ val-single_phot[-1] for val in single_phot ]
    

np.save( open('all_models_trim.npy', 'w'), np.array(all_colors) )
            