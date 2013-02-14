'''
Library to produce a catalog of fully-populated SEDs and calculating
zeropoints for any arbitrary location on the sky, using a combination 
of online catalogs (USNOB1, 2MASS, SDSS) and synthetic photometry.

Requires:
- all_models_P.npy, a file produced by assemble_models.py

To do:
- can improve speed by having each fitting task do a tree search amongst
  models, instead of testing every single model.e
'''


############################################
# IMPORTS
############################################
import numpy as np
from subprocess import Popen, PIPE
from scipy.optimize import fmin_bfgs
from os.path import isfile
from threading import Thread
from time import time, strftime

import multiprocessing as mp
N_CORES = mp.cpu_count()  # use all the cpus you have

try:
    MODELS = np.load( open('all_models_yVega_P.npy','r') )
    # rezero so that K=0 for all models (makes fitting faster)
    for row in MODELS[1:]:
        row[1:] = row[1:] - row[-1]
except:
    raise IOError('cannot find models file')


ALL_FILTERS = ['u','g','r','i','z','y','B','R','J','H','K']
# band: (central wavelength (AA), zeropoint (erg/s/cm^2/AA), catalog index)
FILTER_PARAMS =  {'u': (3551., 8.5864e-9, 0), 'g': (4686., 4.8918e-9, 1),
                  'r': (6165., 2.8473e-9, 2), 'i': (7481., 1.9367e-9, 3),
                  'z': (8931., 1.3564e-9, 4), 'y': (10091., 1.0696e-9, 5),
                  'B': (4400., 6.6000e-9, 6),  'R':(6500., 2.1900e-9, 7),
                  'J':(12350., 3.1353e-10, 8), 'H':(16620., 1.1121e-10, 9),
                  'K':(21590., 4.2909e-11, 10)}

############################################
# ONLINE CATALOG STUFF
############################################

class online_catalog_query():
    '''
    A class to handle all queries of remote catalogs.
    The main function, query_all(), queries all catalogs in a
     multithreaded manner, to minimize lag from I/O communications.
     
    Standard usage:
     q = online_catalog_query( ra, dec, field_size ) #ra,dec in decimal degrees, field_size in arcsec
     Mass, SDSS, USNOB1 = q.query_all()
     # OR #
     Mass = q.query_2mass() #same for all catalogs
    '''
    def __init__(self, ra, dec, boxsize=10. ):
        self.coords = (ra, dec)
        self.boxsize = boxsize
    
    
    def _parse_sdss( self, s ):
        '''
        Parse findsdss8 output.
        
        s: a string as returned by findsdss8 using the flag "-e0,"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, u, u_sigma, g, g_sigma, r, r_sigma, i, i_sigma, z, z_sigma]
        ''' 
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = [ val for val in line[4:].split(' ') if val ] #drop first part to avoid inconsistency
                # RA and DEC
                if '+' in line[0]:
                    char = '+'
                else: 
                    char = '-'
                ra  = float(line[0].split(char)[0])
                dec = float(char + line[0].split(char)[1])
                # magnitudes and errors
                #  in order: u, u_sig, g, g_sig, r, r_sig, i, i_sig, z, z_sig 
                tmp = [ra, dec]
                for band in line[3:8]:
                    tmp += map( float, band.split(':') )
                out.append( tmp )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        return np.array(out)
    
    
    def query_sdss( self, container=None, cont_index=1, trim_mag=21. ):
        '''
        Query sdss8 server for sources found in a box of width boxsize (arcsecs)
        around ra,dec.
        Returns an array (if objects found) or None (if not)
        
        ra,dec: coordinates in decimal degrees
        boxsize: width of box in which to query
        trim_mag: do not return sources with r > trim_mag
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for SDSS objects around coordinates with
        #  defined box size, return only basic parameters, and
        #  sort by distance from coordinates, and return only those
        #  sources brighter than trim_mag
        request = 'findsdss8 -c "{} {}" -bs {} -lmr 0,{} -e0 -sr -m 50000'.format( ra, dec, boxsize, trim_mag )
        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()
        # parse the response
        sdss_objects = self._parse_sdss(o)
        if len(sdss_objects) == 0:
            # no matches
            output = None
        else:
            sdss = np.array( [obj for obj in sdss_objects if obj[6] < trim_mag] )
            output = sdss
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def _parse_2mass( self, s ):
        '''
        parse find2mass output.
        
        s: a string as returned by find2mass using the flag "-eb"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, J, J_sigma, H, H_sigma, K, K_sigma]
        '''
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = line.split('|')
                # RA and DEC
                ra, dec = map( float, line[0].split(' ') )
                # magnitudes and errors
                #  in order: J, J_sig, H, H_sig, K, K_sig
                tmp = [ra, dec]
                for band in line[2:5]:
                    mag, err = [val for val in band.split(' ') if val]
                    mag = float(mag)
                    if '-' not in err:
                        # some objects do not have reported errors - for these,
                        #  assume a conservative error of 0.25mag
                        err = float(err)
                    else:
                        err = .25
                    tmp += [mag, err]
                out.append( tmp )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        return np.array(out)
    
    
    def query_2mass( self, container=None, cont_index=0 ):
        '''
        Query 2mass server for sources found in a box of width boxsize (arcsecs)
        around ra,dec.
        Returns an array (if objects found) or None (if not)
        
        ra,dec: coordinates in decimal degrees
        boxsize: width of box in which to query
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for 2Mass point sources around coordinates with
        #  defined box size, return only basic parameters, and 
        #  sort by distance from coordinates, and return a 
        #  maximum of 10000 sources (i.e. return everything)
        request = 'find2mass -c {} {} -bs {} -eb -sr -m 50000'.format( ra, dec, boxsize )
        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()
        # parse the response
        mass_objects = self._parse_2mass(o)
        if len(mass_objects) == 0:
            # no matches
            output = None
        else:
            output = mass_objects
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def _parse_usnob1( self, s ):
        '''
        Parse findusnob1 output.
        The photometric errors in USNOB1 are pretty terrible (~.3mag for each observation),
          but most sources have more than one observation, so I average all results
          and return the average magnitudes and the error of the mean.
        
        s: a string as returned by usnob1 using the flag "-eb"
        returns: a 2-d array, with each row containing the results for one object:
          [ra, dec, avg_B, B_sigma, avg_R, R_sigma]
        '''
        # parse the header to see how many of each magnitude are reported
        header = s.split('\n')[3]
        obs_count = header.count('Bmag')
        # just as a sanity check, make sure it reports the same number of B,R mags
        assert( obs_count == header.count('Rmag') )
        
        out = []
        lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
        for line in lines:
            try:
                line = line.split('|')
                # RA and DEC
                tmp = [ val for val in line[0].split(' ') if val ]
                if '+' in tmp[1]:
                    char = '+'
                else: 
                    char = '-'
                ra  = float(tmp[1].split(char)[0])
                dec = float(char + tmp[1].split(char)[1])
            
                # magnitudes and errors
                #  in order: B, B_sigma, R, R_sigma
                Bs, Rs = [], []
                for i in range(obs_count):
                    tmp = line[1 + 2*i]
                    if '-' not in tmp:
                        Bs.append( float(tmp) )
                    tmp = line[2 + 2*i]
                    if '-' not in tmp:
                        Rs.append( float(tmp) )
                # ignore sources that don't include at least one B and one R observation
                if not Bs or not Rs: continue
                B = np.mean(Bs)
                R = np.mean(Rs)
                # each measure has an error of about .3 mag
                B_err = 0.3/np.sqrt(len(Bs))
                R_err = 0.3/np.sqrt(len(Rs))
            
                out.append( [ra, dec, B, B_err, R, R_err] )
            except:
                # silently fail on sources that are not formatted properly,
                #  they are probably anomalous anyways.
                pass
        
        return np.array(out)
    
    
    def query_usnob1( self, container=None, cont_index=2 ):
        '''
        Query usnob1 server for sources found in a box of width boxsize (arcsecs)
        around ra,dec.
        Returns an array (if objects found) or None (if not)
        
        ra,dec: coordinates in decimal degrees
        boxsize: width of box in which to query
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # search for USNOB1 point sources around coordinates with
        #  defined box size, return only basic parameters, and 
        #  sort by distance from coordinates, and return a maximum
        #  of 10000 objects (i.e. return everything)
        request = 'findusnob1 -c {} {} -bs {} -eb -sr -m 50000'.format( ra, dec, boxsize )
        out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
        o,e = out.communicate()
        usnob1_objects = self._parse_usnob1(o)
        if len(usnob1_objects) == 0:
            # no matches
            output = None
        else:
            output = usnob1_objects
        if container == None:
            return output
        else:
            container[ cont_index ] = output
    
    
    def query_all( self ):
        '''
        Query all sources, with an independent thread for each
         so that the communications happen concurrently, to save
         time.
        
        returns: [2Mass, SDSS, USNOB1]
        '''
        ra,dec  = self.coords
        boxsize = self.boxsize
        # results is a container into which the threads will put their responses
        results = [None]*3
        t0 = Thread(target=self.query_2mass, args=(results,))
        t1 = Thread(target=self.query_sdss, args=(results,))
        t2 = Thread(target=self.query_usnob1, args=(results,))    
        threads = [t0, t1, t2]
        for t in threads:
            t.start()
        for t in threads:
            # join each thread and wait until each is completed
            t.join()
        return results
    


def identify_matches( queried_stars, found_stars, match_radius=1., chunk_size=300 ):
    '''
    Match queried stars to found stars.
    
    Returns list of indices in found_stars for corresponding
     star in queried_stars, or None if no match found.
     
    queried_stars: an array of coordinates of stars to match
    found_stars: an array of coordinates of located stars
    match_radius: the maximum offset between queried and found star to 
      call a match, in arcseconds
    chunk_size: size (in arcseconds) of each chunk to process at a time.
      (since matching time goes as N^2)
    '''
    mr = 2.778e-4*match_radius # convert arcseconds into degrees
    match_radius_squared = (mr)**2
    
    # find the field that encloses all of the queried sources
    field_center, field_width = find_field( queried_stars )
    # split into seperate chunks
    chunks, size = split_field( field_center, field_width, chunk_size, object_coords=queried_stars )
    # convert size to degrees, for use below
    size = size*2.778e-4
    
    # go through each chunk, filling in matches as they're found
    matches = [None]*len(queried_stars)
    for chunk in chunks:
        ra_range = (chunk[0]-size-mr, chunk[0]+size+mr)
        dec_range = (chunk[1]-size-mr, chunk[1]+size+mr)
        
        # find all the stars in this chunk, keeping track of their original indices
        qs = [(i,star) for i,star in enumerate(queried_stars) if (ra_range[0] < star[0] < ra_range[1]) and (dec_range[0] < star[1] < dec_range[1])]
        fs = [(i,star) for i,star in enumerate(found_stars) if (ra_range[0] < star[0] < ra_range[1]) and (dec_range[0] < star[1] < dec_range[1])]
        
        # now go through and actually find the matches
        for q_index,star in qs:
            diffs_squared = [ (star[0]-other[0])**2 + (star[1]-other[1])**2 for tmp,other in fs ]
            if min(diffs_squared) < match_radius_squared:
                i_best = np.argmin(diffs_squared)
                matches[q_index] = fs[i_best][0] # record, in the proper spot, the original index of the matched star
            else:
                pass
    return matches


def find_field( star_coords, extend=.0015 ):
    '''
    determine the best field for a list of star coordinates,
    so we can perform only a single query.
    
    star_coords: a list of star coordinates, in decimal degrees
    extend: the buffer beyond the requested coordinates to add to the field
      (also in decimal degrees)
    
    returns: (coordinates of center in decimal degrees), (RA_width_of_box, DEC_width_of_box) (in arcseconds)
    '''
    ras = star_coords[:,0]
    decs = star_coords[:,1]
    width_ra = (max(ras)-min(ras) + 2*extend)
    center_ra = np.mean(ras)
    width_dec = (max(decs)-min(decs) + 2*extend)
    center_dec = np.mean(decs)
    return (center_ra, center_dec), (width_ra*3600, width_dec*3600.)


def split_field( field_center, field_width, max_size, object_coords=None ):
    '''
    split a single field into multiple smaller chunks.
    
    field_center: (ra, dec) in decimal degrees
    field_width: each side of box, in arcseconds (RA_width, DEC_width)
    max_size: maximum size of tile, in arcseconds
    object_coords: array of coordinates - if given, will only return chunks
        that contain at least one of these sources
        
    Returns: list of new field centers, width of new fields (in arcseconds)
    '''
    a2d = 2.778e-4 #conversion between arcseconds & degrees
    
    R_n_tile = np.ceil(field_width[0]/max_size).astype(int) # number of tiles needed in RA
    D_n_tile = np.ceil(field_width[1]/max_size).astype(int) # number of tiles needed in DEC
    # use the width implied by larger of the two dimensions
    w_tile = max( (field_width[0]/R_n_tile,field_width[1]/D_n_tile) )*a2d # width of each tile in degrees
    
    ra,dec = field_center
    centers = []
    rr = ra - a2d*field_width[0]/2. + w_tile/2.  # the beginning positions of the tiling
    for i in range(R_n_tile):
        dd = dec - a2d*field_width[1]/2. + w_tile/2.
        for j in range(D_n_tile):
            centers.append( (rr, dd) )
            dd += w_tile
        rr += w_tile
        
    if object_coords != None:
        # keep only those that contain sources
        RAs = object_coords[:,0]
        DECs = object_coords[:,1]
        good_centers = []
        for cent in centers:
            tf_array = (np.abs(RAs - cent[0]) < w_tile/2.) & (np.abs(DECs - cent[1]) < w_tile/2.)
            if any(tf_array):
                good_centers.append(cent)
        centers = good_centers
    
    return centers, w_tile/a2d


############################################
# MODEL-FITTING FUNCTIONS
############################################

def _error_C(C, model, obs, weights):
    '''
    C: number, a constant akin to distance modulus
    model: array-like, model mags
    real: array-like, observed mags
    
    returns: sum squared errors
    '''
    nm = model+C
    return np.sum( weights*(nm-obs)**2 )


def choose_model( obs, mask, models=MODELS ):
    '''
    Find and return the best model for obs.
    Do this by fitting to all magnitudes weighted by error.
    
    Returns: model, temperature, quality_parameter
    
    obs: an array of the observed magnitudes and errors, in
           order as defined by the mode key (see below)
    mask: defines what colors to use in the fit, i.e. what observations exist
           in order: [u,g,r,i,z,y,B,R,J,H,K]
    models: an array of modeled SEDs, where 0th entry is temperature
             of the model, and the rest are magnitudes
    '''
    # mask is an array used to choose which modeled magnitudes
    #  correspond to the included observations.
    #  Note: I append a leading zero to ignore the temperature of the model
    #   (saved in position 0 for all models)
    mask = np.hstack( (np.array([False]), np.array(mask).astype(bool)) )
    
    mags = obs[::2]
    Jmag = mags[-1]
    zerod_mags = mags - Jmag # recenter to compare to models
    weights = 1./obs[1::2]
    
    # Go through all models and choose the one with the most similar SED
    #  Keep track of the sum_squared error, as returned by _error_C()
    sum_sqrs, Cs = [], []
    for model in models[1:]:
        res = fmin_bfgs( _error_C, 0., args=(model[mask], zerod_mags, weights), full_output=True, disp=False )
        Cs.append(res[0][0])
        sum_sqrs.append(res[1])
    
    i_best = np.argmin(sum_sqrs)
    best_model = models[1:][ i_best ]
    
    # if there's one point more than <max_diff> away from best model, try again without that point
    #  (set that weight to zero, and return the index of the cut value)
    i_cut = None
    if max( np.abs(best_model[mask] - zerod_mags) ) > 1.:
        i_cut = np.argmax(np.abs(best_model[mask] - zerod_mags))
        weights[i_cut] = 0.
        # Go through all models again, with new weights
        sum_sqrs, Cs = [], []
        for model in models[1:]:
            res = fmin_bfgs( _error_C, 0., args=(model[mask], zerod_mags, weights), full_output=True, disp=False )
            Cs.append(res[0][0])
            sum_sqrs.append(res[1])
            
        i_best = np.argmin(sum_sqrs)
        best_model = models[1:][ i_best ]
        
    # now add back in the Jmag value to get a model for the non-zeroed observations
    C = Cs[i_best] + Jmag
    
    # return all magnitudes for best model, the offset C, the index, and a quality metric for the best fit
    #  The quality metric is the average error between the best model and the observations
    weights = np.ones_like(mags)
    if i_cut != None:
        # handle the errors properly if a value was dropped
        weights[i_cut] = 0.
        sum_sqr_err = _error_C( C, best_model[mask], mags, weights )
        metric = (sum_sqr_err/(len(mags)-1))**.5
    else:
        sum_sqr_err = _error_C( C, best_model[mask], mags, weights )
        metric = (sum_sqr_err/len(mags))**.5
    return (best_model[1:] + C, C, best_model[0], metric, i_cut )


def fit_sources( inn ):
    '''
    Wrapper function for choose_model to facilitate multiprocessor
     use ( in catalog.produce_catalog() )
    '''
    # the masks show, in order, which bands are included
    #  order: u,g,r,i,z, y, B,R, J,H,K
    mode, obs = inn
    if mode == 0: # sdss+2mass
        mask = [1,1,1,1,1, 0, 0,0, 1,1,1]
    elif mode == 1: # usnob+2mass
        mask = [0,0,0,0,0, 0, 1,1, 1,1,1]
    mask = np.array(mask).astype(bool)
    model, C, index, err, i_cut = choose_model( obs, mask )
    
    sed = np.empty(11)
    full_errs = np.empty(11)
    # keep the real observations
    sed[mask] = obs[::2]
    full_errs[mask] = obs[1::2]
    # fill in rest with modeled magnitudes
    sed[~mask] = model[~mask]
    full_errs[~mask] = err
    
    # if a value was cut while fitting, return the modeled magnitude instead of the observed
    if i_cut != None:
        # do some gymnastics to get the cut passband since it's behind a mask
        cut_band = ALL_FILTERS.index( np.array(ALL_FILTERS)[mask][i_cut] )
        sed[cut_band] = model[mask][i_cut]
        full_errs[cut_band] = err
        
    # put in the hack to recenter g,r (if in mode 1)
    if mode == 1:
        sed[1] += .227 #g
        sed[2] += .114 #r
        
    return ( sed, full_errs, err, index )


############################################
# SYNTHETIC CATALOG STUFF
############################################

class catalog():
    '''
    A class to handle all catalog creation.
    The main function produces, for any requested field,
     a catalog of cross-matched sources from
     2-MASS and (SDSS or USNOB1), with any missing
     bands filled in through modeling.
    
    Standard usage:
     c = catalog( (ra, dec), field_size )  # with ra, dec in degrees and field_size in arcseconds
     catalog_coords = c.coords
     catalog_SEDs = c.SEDs
    
    Optional keywords:
     ignore_sdss: if True, will not use any SDSS magnitudes in modeling. Used for testing.
    '''
    MAX_SIZE = 7200 # max size of largest single query
    ERR_CUT  = .5   # maximum fitting error allowed
    
    def __init__( self, field_center, field_width, ignore_sdss=False ):
        self.field_center = field_center
        self.field_width = field_width
        self.coords = []
        self.SEDs = []
        self.full_errors = []
        self.model_errors = []
        self.models = []
        self.modes = []
        
        if field_width > self.MAX_SIZE:
            # simply don't allow queries that are too large
            raise ValueError( 'Field is too large. Max allowed: {}"'.format(self.MAX_SIZE) )
        else:
            self.produce_catalog( ignore_sdss=ignore_sdss )
    
    
    def produce_catalog( self, ignore_sdss=False ):
        '''
        Create a catalog of all objects found in field.
        Requires records in 2MASS + (SDSS and/or USNOB1).
        '''
        ra,dec = self.field_center
        q = online_catalog_query( ra, dec, self.field_width )
        mass, sdss, usnob = q.query_all()
        
        object_mags = []
        modes = []
        object_coords = []
        if mass != None:
            # match sdss, usnob objects to 2mass objects
            if sdss != None and not ignore_sdss:
                sdss_matches = identify_matches( mass[:,:2], sdss[:,:2] )
            else:
                sdss_matches = [None]*len(mass)
            if usnob != None:
                usnob_matches = identify_matches( mass[:,:2], usnob[:,:2] )
            else:
                usnob_matches = [None]*len(mass)
            
            # Go through 2mass objects and assemble a catalog
            #  of all objects present in 2mass and (sdss or usnob)
            #  Use 2mass+sdss OR 2mass+usnob (ignore usnob if sdss present)
            for i,obj in enumerate(mass):
                if sdss_matches[i] != None:
                    i_sdss = sdss_matches[i]
                    obs = np.hstack( (sdss[i_sdss][2:], obj[2:]) )
                    mode = 0
                elif sdss_matches[i] == None and usnob_matches[i] != None:
                    i_usnob = usnob_matches[i]
                    obs = np.hstack( (usnob[i_usnob][2:], obj[2:]) )
                    mode = 1
                elif sdss_matches[i] == None and usnob_matches[i] == None:
                    continue
                object_mags.append( obs )
                modes.append( mode )
                object_coords.append( obj[:2] )
        else:
            raise ValueError( "No 2-MASS sources in this field!" )
        
        # send all of these matches to the CPU pool to get modeled
        objects = zip( modes, object_mags )
        pool = mp.Pool( processes=N_CORES )
        results = pool.map( fit_sources, objects )
        pool.close()
        
        # now go through results and construct the final values
        for i,row in enumerate(results):
            # each row is (sed, full_errs, model_err, index)
            if row[2] > self.ERR_CUT: #apply quality-of-fit cut
                pass
            else:
                self.coords.append( object_coords[i] )
                self.SEDs.append( row[0] )
                self.full_errors.append( row[1] )
                self.model_errors.append( row[2] )
                self.models.append( row[3] )
                self.modes.append( modes[i] )
        # keep the multi-dimensional data in numpy arrays
        self.coords = np.array(self.coords)
        self.SEDs = np.array(self.SEDs)
        self.full_errors = np.array(self.full_errors)
    
    
    def save_catalog( self, file_name ):
        save_catalog( self.coords, self.SEDs, self.full_errors, self.modes, file_name )
    
    


def save_catalog( coordinates, seds, errors, modes, file_name ):
    '''
    Save an output ASCII file of the catalog.
    file_name: output file to create
    '''
    fff = open(file_name,'w')
    fff.write('# Observed/modeled SEDs produced by get_SEDs.py \n' +
              '# Generated: {}\n'.format(strftime("%H:%M %B %d, %Y")) +
              '# modes: 0 -> SDSS+2MASS; 1 -> USNOB1+2MASS\n' +
              "# " + "{}      {}       ".format("RA","DEC") + (' '*6).join(ALL_FILTERS) + \
              " "*6 + ' '.join([val+"_err" for val in ALL_FILTERS]) + "  Mode\n")
    for i,row in enumerate(seds):
        row_txt = " ".join(map(lambda x: "%.6f"%x, coordinates[i]))+" "+ \
                  " ".join(map(lambda x: "%.3f"%x, row))+" "+ \
                  " ".join(map(lambda x: "%.3f"%x, errors[i]))+" "+ \
                  str(modes[i])+"\n"
        fff.write( row_txt )
    fff.close()


############################################
# ZEROPOINT STUFF
############################################

def calc_zeropoint( input_coords, catalog_coords, input_mags, catalog_mags, clip=True, sig_clip=3., max_iter=5, convergence=.02, return_zps=False ):
    '''
    Calculate the zeropoint for a set of input stars and set of catalog stars.
    
    input_coords: a 2D array of [ [RA, DEC], [..., ...] ... ]
    catalog_coords: similar array for objects created with catalog()
    input_mags: a 1D array of instrumental magnitudes
    catalog_mags: a similar array of true magnitudes as created with catalog()
    sig_clip: iteratively trim values more than sig_clip*std away from median
    max_iter: stop the above sigma clipping after max_iter iterations
    convergence: the fractional convergence required
    return_zps: return all zeropoint estimates, not just the median
    
    Returns: zeropoint (mags), the median average deviation, and a list of matched indices for input and catalog sources.
    '''
    matches = identify_matches( input_coords, catalog_coords )
    matched_inputs = [input_mags[i] for i in range(len(matches)) if matches[i] != None ]
    matched_catalogs = [catalog_mags[ matches[i] ] for i in range(len(matches)) if matches[i] != None ]
    
    zp_estimates = []
    for i,inst_mag in enumerate(matched_inputs):
        zp_estimates.append( matched_catalogs[i] - inst_mag )
    zp = np.array(zp_estimates)
    if clip:
        # perform iterative sigma clipping
        for count in range(max_iter):
            in_len = len(zp)
            zp = zp[ np.abs(zp-np.median(zp)) < sig_clip*np.std(zp) ]
            if float(in_len-len(zp))/in_len < convergence:
                break
    mad = np.median( np.abs( zp-np.median(zp) ) )
    if return_zps:
        return np.median(zp), mad, matches, zp
    else:
        return np.median(zp), mad, matches
        


def zeropoint( input_file, band, output_file=None ):
    '''
    Calculate <band> zeropoint for stars in <input_file>.
    
    Expects a space-or-tab-delimited ascii input file with the
     first column RA, the second DEC, and the third instrumental magnitude.
     Header/comments should be #-demarcated, and all non-commented rows in the
     file should be numbers only.
    If an output_file name is given, it saves the entire catalog to that file.
    '''
    # load the data and produce a catalog
    in_data = np.loadtxt( input_file )
    input_coords = in_data[:, :2]
    input_mags = in_data[:, 2]
    field_center, field_width = find_field( input_coords )
    c = catalog( field_center, max(field_width) )
    
    # identify which band this is for and calculate the zeropoint
    cat_index = FILTER_PARAMS[band][2]
    cat_mags = c.SEDs[:, cat_index]
    zp, mad, matches = calc_zeropoint( input_coords, c.coords, input_mags, cat_mags )
    
    if output_file:
        # save matched catalog to file
        oc, os, oe, om = [],[],[],[]
        for i,match in enumerate(matches):
            oc.append( input_coords[i] )
            if match:
                os.append( c.SEDs[match] )
                oe.append( c.full_errors[match] )
                om.append( c.modes[match] )
            else:
                os.append( [99]*len(ALL_FILTERS) )
                oe.append( [9]*len(ALL_FILTERS) )
                om.append( -1 )
        save_catalog( oc, os, oe, om, output_file )
    
    return zp, mad

