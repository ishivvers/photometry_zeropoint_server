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

from get_SEDs import _find_field, _split_field  #need to import helpers explicitly
from get_SEDs import *
import numpy as np
import matplotlib.pyplot as plt


def get_Y_match( object_coords, DEC_obs, max_size=900., plot=False):
    '''
    Estimate y-band photometry for sources that cross-match to 2MASS
    
    DEC_obs: array of [[u, u_err, g, g_err, r, r_err, i, i_err, z, z_err, y, y_err], ...]
    '''
    
    # only bother with sources z<20., since any dimmer than that will not be in
    #  the 2MASS catalog
    good_obs, good_coords = [],[]
    for i,row in enumerate(DEC_obs):
        if row[8]<20.:
            good_obs.append( (row,i) ) #keep track of modeled sources original indices
            good_coords.append( object_coords[i] )
    good_coords = np.array(good_coords)  #keep good_obs as list, because of the indexer
    
    field_center, field_width = _find_field( good_coords )
    # split field and query/match for each tile
    #  works even if we don't need to split
    centers, tile_width = _split_field( field_center, field_width, max_size, object_coords=good_coords )
    
    print 'Matching to 2MASS...'
    # go through each tile and accumulate the results:
    matched_coords, matched_obs, matched_indexer = [],[],[]
    for center in centers:
        mass = query_2mass( center[0], center[1], boxsize=tile_width )
        if mass != None:
            # produce list of matched objects
            #  If can figure a better way to track indices, could make this faster by
            #  only attempting to match to sources in this queried field.
            matches = identify_matches( mass[:,:2], good_coords )
            if matches != None:
                for i,obj in enumerate(mass):
                    if matches[i] != None:
                        i_match = matches[i]
                        i_orig = good_obs[i_match][1]
                        if i_orig in matched_indexer:
                            # we already matched this source. continue for now,
                            #  perhaps there's a more graceful way to handle this?
                            continue
                        matched_indexer.append( i_orig )
                        mobs = np.hstack( (good_obs[i_match][0], obj[2:]) )
                        matched_obs.append( mobs )
                        matched_coords.append( good_coords[i_match] )
    
    # now go through the matched sources and model each, predicting the missing mags
    modeled_mags, modeled_errs = [],[]
    print 'Modeling matches...'
    for i, obs in enumerate(matched_obs):
        if i%10==0: print i
        
        good_obs = []
        mask = [1,1,1,1,1, 0, 0,0, 1,1,1] #expected
        # if any of the orignal obs are not actually observed (value == 99.), mask them as well
        #  this builds up an array of good observations and associated errors
        for j,val in enumerate(obs[::2]):
            if val == 99.:
                mask[j] = 0 # this is irrelevant for the 2MASS obs, i.e. j>6
            else:
                good_obs += [val, obs[1::2][j]]
        good_obs = np.array(good_obs)
        mask = np.array(mask).astype(bool)
        
        model, T, err = choose_model( good_obs, mask )

        final_mags, final_errs = [],[]
        for j,val in enumerate(obs[::2][:6]): #trim 2MASS mags here, with the :6 cut
            if val == 99.:
                #replace this value with modeled
                final_mags.append( model[j] )
                final_errs.append( err )
            else:
                #keep the observed value
                final_mags.append( val )
                final_errs.append( obs[1::2][j] )
        
        final_mags = np.array(final_mags)
        final_errs = np.array(final_errs)
        modeled_mags.append( final_mags )
        modeled_errs.append( final_errs )
        
        # now plot up a few, to see how they look:
        if plot:
            if i<plot:
                plt.scatter( MODELS[0][1:][mask], good_obs[::2], c='r' )
                plt.scatter( MODELS[0][1:], model, c='k' )
                plt.xlabel('wavelength')
                plt.ylabel('mag')
                plt.title('Temp: {} --- Error: {}'.format(T,err))
                plt.gca().invert_yaxis()
                plt.show()
    
    #piece together a complete list of all sources again,
    #  inserting the modeled ones where relevant
    out_obs, out_errs = [],[]
    for i,row in enumerate( DEC_obs ):
        if i in matched_indexer:
            # replace this row with modeled row
            i_modeled = matched_indexer.index(i)
            out_obs.append( modeled_mags[i_modeled] )
            out_errs.append( modeled_errs[i_modeled] )
        else:
            out_obs.append( row[::2] )
            out_errs.append( row[1::2] )
                
    return object_coords, out_obs, out_errs


def get_Y_all( coords, obs, plot=False):
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
                plt.title('Temp: {} --- Error: {}'.format(T, err))
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
        if sys.argv[1] not in fname or 'sources' not in fname:
            continue
        print '*'*20
        print 'working on', fname
        print '*'*20
        data = np.loadtxt(prefix+fname)
        data[ data==0. ] = .00001  # hack to fix reported zero errors
        '''
        #for get_Y_all:
        # dec_obs = [u, u_err, g, g_err, r, r_err, i, i_err, z, z_err]
        dec_obs = np.zeros_like( data[:,:10] ) #prepare the array
        dec_obs[:,::2] = data[:,2:7]
        dec_obs[:,1::2] = data[:,8:13]
        '''
        #for get_Y_match:
        # dec_obs = [u, u_err, g, g_err, r, r_err, i, i_err, z, z_err, y, y_err]
        dec_obs = np.zeros_like( data[:,:12] ) #prepare the array
        dec_obs[:,::2] = data[:,2:8]
        dec_obs[:,1::2] = data[:,8:14]
        
        #out_coords, obj_mags, errors = get_Y_all( data[:,:2], dec_obs )
        out_coords, obj_mags, errors = get_Y_match( data[:,:2][:1000], dec_obs[:1000], plot=10 )
        
        #save to file
        fff = open(prefix+fname.split('.')[0]+'_match.catalog', 'w')
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
        
        