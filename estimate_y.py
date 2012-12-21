'''
Script to use the get_SEDs library to estimate the y-band
magnitude for objects in catalogs.

# reading in DEC_obs and formatting them:
dat = np.loadtxt('catalogs/C1_stars_new.catalog')
dec_obs = np.zeros_like( dat[:,:8] )
dec_obs[:,::2] = dat[:,3:7]
dec_obs[:,1::2] = dat[:,9:13]

# reading in from sources, not just stars
dat = np.loadtxt('catalogs/C1_sources_new.catalog')
dec_obs = np.zeros_like( dat[:,:8] )
dec_obs[:,::2] = dat[:,3:7]
dec_obs[:,1::2] = dat[:,9:13]

good_coords = dat[:,:2][ dec_obs[:,-2] < 20. ]
good_obs = dec_obs[ dec_obs[:,-2] < 20. ]
'''

# first, produce photometry for only those that cross-match
#  to 2MASS

from get_SEDs import _find_field, _split_field, _get_reddening  #need to import helpers explicitly
from get_SEDs import *
import numpy as np
import matplotlib.pyplot as plt


def get_Y_wth_match( object_coords, DEC_obs, err_cut=.5, redden=False, max_size=900., plot=False):
    '''
    Estimate y-band photometry for sources that cross-match to 2MASS
    
    DEC_obs: array of [[g, g_err, r, r_err, i, i_err, z, e_err], ...]
    '''
    field_center, field_width = _find_field( object_coords )
    # split field and query/match for each tile
    #  works even if we don't need to split
    centers, tile_width = _split_field( field_center, field_width, max_size, object_coords=object_coords )
    
    print 'Matching to 2MASS...'
    # go through each tile and accumulate the results:
    obj_coords, obj_mags, = [],[]
    for center in centers:
        mass = query_2mass( center[0], center[1], boxsize=tile_width )
        if mass != None:
            # produce list of matched objects
            matches = identify_matches( mass[:,:2], object_coords )
            if matches != None:
                for i,obj in enumerate(mass):
                    if matches[i] != None:
                        i_match = matches[i]
                        obj_mags.append( np.hstack( (DEC_obs[i_match], obj[2:]) ) )
                        obj_coords.append( object_coords[i_match] )
    
    # now go through the matched sources and model each, predicting the Ymag
    out_coords, out_mags, Ys, errors = [],[],[],[]
    print 'Modeling matches...'
    for i, obs in enumerate(obj_mags):
        if i%100==0: print i
        
        mask = np.array( [0,1,1,1,1, 0, 0,0, 1,1,1] ).astype(bool) #expected observations
        # if any of these are not observed (value == 99.), mask them as well
        ##  DOES NOT WORK DOES NOT WORK
        mask[mask][ obs[::2] == 99. ] = False
        
        if redden:
            reddening = _get_reddening( obj_coords[i][0],obj_coords[i][1], ALL_FILTERS )
            # de-redden the observations before comparing
            #  to the models
            obs[::2] -= reddening[mask]
        
        model, T, err = choose_model( obs, mask )

        if err > err_cut: #apply quality-of-fit cut
            continue
        else:
            if redden:
                # re-redden the model and observations
                obs[::2] += reddening[mask]
                model += reddening
            
            # return Ymag, and keep track of errors
            Ys.append(model[5])
            errors.append( err )
            out_coords.append( obj_coords[i] )
            out_mags.append( obs )

            if i<10 and plot:
                plt.scatter( MODELS[0][1:][mask], obs[::2], c='g' )
                plt.scatter( MODELS[0][1:], model, c='k' )
                plt.scatter( MODELS[0][1:][5], model[5], c='r' )
                plt.xlabel('wavelength')
                plt.ylabel('mag')
                plt.title('Temp: {} --- Error: {}'.format(T,err))
                plt.gca().invert_yaxis()
                plt.show()
                
    return np.array(out_coords), np.array(out_mags), np.array(Ys), np.array(errors)


def get_Y_all( coords, obs, err_cut=.5, plot=False):
    '''
    Estimate all missing photometry for all sources, whether they cross-match
    or no.  (Do not bother with crossmatch!)
    '''
    out_mags = np.zeros( (len(obs), 6) )
    out_errs = np.zeros( (len(obs), 6) )
    for i in range(len(obs)):
        if not i%100: print i
            
        good_obs = []
        mask = [1,1,1,1,1, 0, 0,0, 0,0,0] #expected
        # if any of these are not observed (value == 99.), mask them as well
        #  this also builds up an array of good observations and associated errors
        for j,val in enumerate(obs[i][::2]):
            if val == 99.:
                mask[j] = 0
            else:
                good_obs += [val, obs[i][1::2][j]]
        good_obs = np.array(good_obs)
        mask = np.array(mask).astype(bool)
        
        model, T, err = choose_model( good_obs, mask )
        
        out_mags[i][:-1] = obs[i][::2]
        for j,val in enumerate(mask[:6]):
            if not val:
                # insert the modeled value if we didn't observe it
                out_mags[i][j] = model[j]
        out_errs[i][:-1] = obs[i][1::2]
        for j,val in enumerate(mask[:6]):
            if not val:
                # insert the modeled error if this wasn't observed
                out_errs[i][j] = err
    
        # now plot up a few, to see how they look:
        if plot:
            if i < plot:
                plt.scatter( MODELS[0][1:][mask], out_mags[i][mask[:6]], c='g', alpha=.7 )
                plt.scatter( MODELS[0][1:], model, c='k', alpha=.7 )
                plt.xlabel('wavelength')
                plt.ylabel('mag')
                plt.title('Error: {}'.format(out_errs[i][0]))
                plt.gca().invert_yaxis()
                plt.show()
    
    return coords, out_mags, out_errs
    
    
if __name__ == '__main__':
    '''
    Command-line tool to run this on fed-in files, saving the output.
    
    right now runs everthing in the prefix folder that matches the input
    string.  (for running pseudoparallel)
    '''
    import os, sys
    prefix = 'catalogs/'
    for fname in os.listdir(prefix):
        if sys.argv[1] not in fname:
            continue
        print '*'*20
        print 'working on', fname
        print '*'*20
        data = np.loadtxt(prefix+fname)
        data[ data==0. ] = .00001  # hack to fix reported zero errors
        # dec_obs = [u, u_err, g, g_err, r, r_err, i, i_err, z, z_err]
        dec_obs = np.zeros_like( data[:,:10] ) #prepare the array
        dec_obs[:,::2] = data[:,2:7]
        dec_obs[:,1::2] = data[:,8:13]
        
        '''
        maskA = data[:,6] < 20.  #only bother to match sources with z<20
        good_coords = data[:,:2][ maskA ]
        good_obs = dec_obs[ maskA ]
        '''
        
        out_coords, obj_mags, errors = get_Y_all( data[:,:2], dec_obs )
        
        #save to file
        fff = open(prefix+fname.split('.')[0]+'_all.catalog', 'w')
        hhh = '# RA   DEC  mag(u g r i z y)   mag_error(u g r i z y)'
        fff.write(hhh+'\n')
        for iii,coord in enumerate(out_coords):
            ra,dec = coord
            fff.write('%.6f '*2 %(ra, dec) )
            u,g,r,i,z,y = obj_mags[iii]
            fff.write('%.5f '*6 %( u,g,r,i,z,y ) )
            ue,ge,re,ie,ze,ye = errors[iii]
            fff.write('%.5f '*6 %( ue,ge,re,ie,ze,ye ) + '\n')
        fff.close()
        
        