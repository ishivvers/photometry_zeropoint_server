from get_SEDs import _find_field, _split_field  #need to import helpers explicitly
from get_SEDs import *
import numpy as np
import matplotlib.pyplot as plt



def compare_catalogs( source_cat, synth_cat, bands=['g','r','i','z','y'], err_cut=.5 ):
    '''
    Plot histograms of the errors when estimating source magnitudes
    in grizy bands using USNOB+2MASS.
    '''
    data = np.loadtxt( source_cat )
    # trim out sources with z>20., since there's no hope
    #  that they'll be in USNOB or 2MASS
    source_coords, source_obs = [],[]
    for i,row in enumerate(data[:,3:8]):
        if row[3]<20.:
            source_coords.append( data[i,:2] )
            source_obs.append( row )
    source_coords = np.array(source_coords)
    source_obs = np.array(source_obs)
    # obs in format [g,r,i,z,y] with unobserved=99.
    
    data = np.loadtxt( synth_cat )
    synth_coords = data[:,:2]
    synth_obs = data[:,3:8]
    # obs in format [g,r,i,z,y]
    synth_errs = data[:,-1] #average error of model fit, should correlate with quality
    
    # split fields to make cross-matching faster
    field_center, field_width = _find_field( synth_coords )
    centers, tile_width = _split_field( field_center, field_width, 900., object_coords=source_coords )
    tile_width = 2.778e-4*tile_width #convert to degrees
    
    matched_mags, matched_errs = [],[]
    for i,center in enumerate(centers[:2]):
        print 'field',i,'of',len(centers)

        # find all synths within the field
        rmin, rmax = center[0]-tile_width, center[0]+tile_width
        dmin, dmax = center[1]-tile_width, center[1]+tile_width
        in_field_mask = [1 if ((rmin<coord[0]<rmax) and (dmin<coord[1]<dmax)) else 0 for coord in synth_coords]
        in_field_mask = np.array(in_field_mask).astype(bool)
        test_coords = synth_coords[ in_field_mask ]
        
        matches = identify_matches( test_coords, source_coords )
        
        test_obs = synth_obs[ in_field_mask ]
        test_errs = synth_errs[ in_field_mask ]
        for i, synth in enumerate(test_obs):
            if matches[i] != None:
                imatch = matches[i]
                matched_mags.append( (synth, source_obs[imatch]) )
                matched_errs.append( test_errs[i] )
    
    # now go through and build lists of all of the errors
    errors = [ [],[],[],[],[] ]
    model_fit_errs = []
    band_map = {0:'g', 1:'r', 2:'i', 3:'z', 4:'y'}
    
    for i,match in enumerate(matched_mags):
        #if matched_errs[i] > err_cut:
        #    # fails error cut --- ignore it.
        #    continue
        for j,band in enumerate(errors):
            if match[1][j] == 99.:
                # wasn't actually observed
                band.append( 99. )
            else:
                band.append( match[0][j] - match[1][j] )
            # synthetic magnitude - observed magnitude
        model_fit_errs.append( matched_errs[i] )

            
    # now, we plot them
    prefix = 'catalog_proc/CK_models/cut_hists/'
    mybins = np.linspace(-2, 2, 50)
    for i,band in enumerate(errors):
        if not band: continue  # no observations in that band at all
        plt.figure()
        plt.hist( band, bins=mybins, alpha=.7 )
        plt.xlabel( band_map[i]+': synthetic - observed')
        plt.ylabel('count')
        plt.title( band_map[i] + '--- Mean error: ' + str(np.mean(band)) )
        #plt.savefig( prefix + band_map[i] + source_cat.split('/')[-1].split('.')[0] + '.png')
    
        #now plot scatterplot of model errors vs. true errors
        plt.figure()
        band = np.array(band)
        model_fit_errs = np.array(model_fit_errs)
        plt.scatter( band[band!=99.], model_fit_errs[band!=99.], alpha=.5 )
        plt.xlabel( 'true errors' )
        plt.ylabel( 'model fit errors' )
        plt.title( '{}'.format(band_map[i]) )
        plt.show()
    
    
    
    
"""
if __name__ == '__main__':
    '''
    Command-line tool to run this on fed-in files, saving the output.
    
    right now runs everthing in the prefix folder that matches the input
    string.  (for running pseudoparallel)
    '''
    import os, sys
    prefix = 'catalogs/'
    for fname in os.listdir(prefix):
        if sys.argv[1] not in fname or 'sources' not in fname and 'synth' not in fname:
            continue
        print '*'*20
        print 'working on', fname
        print '*'*20
        data = np.loadtxt(prefix+fname)
        coords = data[:,:2]
        dec_obs = np.zeros_like( data[:,:12] ) #prepare the array
        dec_obs[:,::2] = data[:,2:8]
        dec_obs[:,1::2] = data[:,8:14]
        
        good_coords = []
        for i,row in enumerate(dec_obs):
            if row[8]<20.:
                good_coords.append( coords[i] )
        good_coords = np.array(good_coords)
        
        field_center, field_width = _find_field( good_coords )
        
        # now make and save a catalog for that location
        out_fname = prefix+fname.split('.')[0]+'_synth.catalog'
        catalog( field_center, field_width, object_coords=good_coords, savefile=out_fname, max_size=1800.)
"""
if __name__ == '__main__':
    '''
    run compare_catalogs on files in the catalogs/ directory,
    saving histograms of the determined errors.
    '''
    import os
    directory = 'catalog_proc/CK_models/'
    files = os.listdir(directory)
    synths = [f for f in files if 'synth' in f]
    sources = [f for f in files if 'sources_new.catalog' in f]
    
    for i in range(len(synths)):
        print '*'*20
        print 'working on',sources[i], synths[i]
        print '*'*20
        compare_catalogs( directory+sources[i], directory+synths[i] )
