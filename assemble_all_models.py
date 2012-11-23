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
import os


################################################
# PASSBANDS
pb_dir = 'transmission_functions/'

# SDSS passbands
u = ps.ObsBandpass('sdss,u')
g = ps.ObsBandpass('sdss,g')
r = ps.ObsBandpass('sdss,r')
i = ps.ObsBandpass('sdss,i')
z = ps.ObsBandpass('sdss,z')

# y-band from Pan-STARRS1
y = ps.FileBandpass(pb_dir+'y_PS1_RSR.dat')

# USNOB1 passbands
B = ps.ObsBandpass('B')
R = ps.ObsBandpass('R')

# 2Mass passbands
J = ps.FileBandpass(pb_dir+'J_2Mass_RSR.dat')
H = ps.FileBandpass(pb_dir+'H_2Mass_RSR.dat')
K = ps.FileBandpass(pb_dir+'K_2Mass_RSR.dat')

passbands = [u,g,r,i,z,y,B,R,J,H,K]

################################################
# SPECTRA

# Using the Castelli & Kurucz models for stars with temperatures
#  ranging from ~3500K to ~50,000K, following suggestions from
#  http://www.stsci.edu/hst/HST_overview/documents/synphot/AppA_Catalogs4.html#48115
logg = 4.5
metal = 0.0
Temps = range(3500, 10000, 100) + range(10000, 20000, 500) + range(20000, 50000, 1000)

save_spectra = True

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

#now go through all temperatures and perform photometry
for T in Temps:
    print T
    spec = ps.Icat( 'ck04models', T, metal, logg )
    if save_spectra:
        np.save( 'model_spectra/{}_{}Z_{}g.npy'.format(T, str(metal).replace('.',''), str(logg).replace('.','')),
                  np.vstack( (spec.wave, spec.flux)) )
    single_phot = []
    for bp in passbands:
        ob = ps.Observation(spec, bp)
        # Use ABmags for ugrizy, but Vegamags for BRJHK
        if bp in [u,g,r,i,z,y]:
            mag = ob.effstim('abmag')
        else:
            mag = ob.effstim('vegamag')
        single_phot.append(mag)
    # Save the model relative to K=0, with the 0th index as the temperature
    single_colors = [T]+[ val-single_phot[-1] for val in single_phot ]
    all_colors.append(single_colors)

np.save( open('all_models.npy', 'w'), np.array(all_colors) )
            