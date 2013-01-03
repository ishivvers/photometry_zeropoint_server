'''
Library to produce a catalog of fully-populated SEDs for
any arbitrary location on the sky, using a combination 
of online catalogs (USNOB1, 2MASS, SDSS) and synthetic photometry.

Requires:
- all_models.npy, a file produced by assemble_models.py
- Schlegel et al. dust extiction maps in a folder named dust_maps

TO DO:
- split field in calculate_zeropoint to decrease running time
- try running everything in Jy instead of flam, to see if that helps with fitting/precision errors
'''


############################################
# IMPORTS
############################################
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from scipy.optimize import fmin_bfgs, fmin_l_bfgs_b
import pyfits
from os.path import isfile
from threading import Thread

try:
    MODELS = np.load( open('all_models_P.npy','r') )
except:
    raise IOError('cannot find models file')

try:
    from _libastro import eq_gal
except:
    raise Exception('cannot find _libastro.so')

# open the dust extinction maps
try:
    hdu = pyfits.open('dust_maps/SFD_dust_1024_ngp.fits')[0]
    nsgp_n = hdu.header['LAM_NSGP']
    scale_n = hdu.header['LAM_SCAL']
    map_n = hdu.data
    hdu = pyfits.open('dust_maps/SFD_dust_1024_sgp.fits')[0]
    nsgp_s = hdu.header['LAM_NSGP']
    scale_s = hdu.header['LAM_SCAL']
    map_s = hdu.data
    MAP_DICT = {}
    MAP_DICT['N'] = [map_n, nsgp_n, scale_n]
    MAP_DICT['S'] = [map_s, nsgp_s, scale_s]
except:
    raise IOError('cannot find/open dust maps')

# Include the E(B-V)-to-mag mapping for all filters,
#  from Schlafly & Finkbeiner 2011
#  Using UKIRT JHK for 2MASS
MAP_DICT['Filters'] = {'u':4.239, 'g':3.303, 'r':2.285, 'i':1.698,
                       'z':1.263, 'y':1.058, 'B':4.315, 'R':2.673,
                       'J':0.709, 'H':0.449, 'K':0.302}

ALL_FILTERS = np.array(['u','g','r','i','z','y','B','R','J','H','K'])
# central wavelength (AA) and flux zeropoint (erg/s/cm^2/A, i.e. flam)
FILTER_PARAMS =  {'u': (3551., 8.5864e-9), 'g': (4686., 4.8918e-9),
                  'r': (6165., 2.8473e-9), 'i': (7481., 1.9367e-9),
                  'z': (8931., 1.3564e-9), 'y': (10091., 1.0696e-9),
                  'B': (4400., 6.6000e-9),  'R':(6500., 2.1900e-9),
                  'J':(12350., 3.1353e-10), 'H':(16620., 1.1121e-10),
                  'K':(21590., 4.2909e-11)}

############################################
# CATALOG INTERFACE FUNCTIONS
############################################

def _parse_sdss( s ):
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


def query_sdss( ra, dec, boxsize=10., container=None, cont_index=1, trim_mag=21. ):
    '''
    Query sdss8 server for sources found in a box of width boxsize (arcsecs)
    around ra,dec.
    Returns an array (if objects found) or None (if not)
    
    ra,dec: coordinates in decimal degrees
    boxsize: width of box in which to query
    trim_mag: do not return sources with r > trim_mag
    '''
    # search for SDSS objects around coordinates with
    #  defined box size, return only basic parameters, and
    #  sort by distance from coordinates, and return a maximum
    #  of 10000 sources (i.e. return everything)
    request = 'findsdss8 -c "{} {}" -bs {} -e0 -sr -m 10000'.format( ra, dec, boxsize )
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    sdss_objects = _parse_sdss(o)
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


def _parse_2mass( s ):
    '''
    parse find2mass output.
    
    s: a string as returned by find2mass using the flag "-eb"
    returns: a 2-d array, with each row containing the results for one object:
      [ra, dec, J, J_sigma, H, H_sigma, K, K_sigma]
    '''
    out = []
    lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
    for line in lines:
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
    '''
        except:
            # silently fail on sources that are not formatted properly,
            #  they are probably anomalous anyways.
            pass
    '''
    
    return np.array(out)


def query_2mass( ra, dec, boxsize=10., container=None, cont_index=0 ):
    '''
    Query 2mass server for sources found in a box of width boxsize (arcsecs)
    around ra,dec.
    Returns an array (if objects found) or None (if not)
    
    ra,dec: coordinates in decimal degrees
    boxsize: width of box in which to query
    '''
    # search for 2Mass point sources around coordinates with
    #  defined box size, return only basic parameters, and 
    #  sort by distance from coordinates, and return a 
    #  maximum of 10000 sources (i.e. return everything)
    request = 'find2mass -c {} {} -bs {} -eb -sr -m 10000'.format( ra, dec, boxsize )
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    mass_objects = _parse_2mass(o)
    if len(mass_objects) == 0:
        # no matches
        output = None
    else:
        output = mass_objects
    if container == None:
        return output
    else:
        container[ cont_index ] = output


def _parse_usnob1( s ):
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


def query_usnob1( ra, dec, boxsize=10., container=None, cont_index=2 ):
    '''
    Query usnob1 server for sources found in a box of width boxsize (arcsecs)
    around ra,dec.
    Returns an array (if objects found) or None (if not)
    
    ra,dec: coordinates in decimal degrees
    boxsize: width of box in which to query
    '''
    # search for USNOB1 point sources around coordinates with
    #  defined box size, return only basic parameters, and 
    #  sort by distance from coordinates, and return a maximum
    #  of 10000 objects (i.e. return everything)
    request = 'findusnob1 -c {} {} -bs {} -eb -sr -m 10000'.format( ra, dec, boxsize )
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    usnob1_objects = _parse_usnob1(o)
    if len(usnob1_objects) == 0:
        # no matches
        output = None
    else:
        output = usnob1_objects
    if container == None:
        return output
    else:
        container[ cont_index ] = output


def query_all( ra, dec, boxsize=10. ):
    '''
    Query all sources, with an independent thread for each
     so that the communications happen concurrently, to save
     time.
    
    returns: [2Mass, SDSS, USNOB1]
    '''
    # results is a container into which the threads will put their responses
    results = [None]*3
    t0 = Thread(target=query_2mass, args=(ra, dec, boxsize, results))
    t1 = Thread(target=query_sdss, args=(ra, dec, boxsize, results))
    t2 = Thread(target=query_usnob1, args=(ra, dec, boxsize, results))    
    threads = [t0, t1, t2]
    for t in threads:
        t.start()
    for t in threads:
        # join each thread and wait until each is completed
        t.join()
    return results


def identify_matches( queried_stars, found_stars, match_radius=1. ):
    '''
    Match queried stars to found stars.
    
    Returns list of indices in found_stars for corresponding
     star in queried_stars, or None if no match found.
     
    queried_stars: an array of coordinates of stars to match
    found_stars: an array of coordinates of located stars
    match_radius: the maximum offset between queried and found star to 
      call a match, in arcseconds
    '''
    match_radius = 2.778e-4*match_radius  # convert arcseconds into degrees
    # find the box that encloses all queried stars, and only match to
    #  stars within that box
    ra_range = (np.min(queried_stars[:,0]-match_radius), np.max(queried_stars[:,0]+match_radius))
    dec_range = (np.min(queried_stars[:,1]-match_radius), np.max(queried_stars[:,1])+match_radius)
    
    match_radius_squared = match_radius**2 
    matches = []
    for star in queried_stars:
        # calculate the distance to each star, but only if it's within the box:
        diffs_squared = []
        for other in found_stars:
            # don't bother with stars outside the box
            if not (ra_range[0] < other[0] < ra_range[1]) \
                    or not (dec_range[0] < other[1] < dec_range[1]):
                diffs_squared.append(999.)
                continue
            ds = (star[0]-other[0])**2 + (star[1]-other[1])**2
            diffs_squared.append(ds)
        
        if min(diffs_squared) < match_radius_squared:
            i_best = np.argmin(diffs_squared)
            matches.append(i_best)
        else:
            matches.append(None) 
    return matches


def produce_catalog( field_center, field_width, err_cut=.5, redden=True, return_models=False ):
    '''
    Create a catalog of all objects found in field.
    Requires records in 2MASS + (SDSS and/or USNOB1).
    
    field_center: (ra, dec) in decimal degrees
    field_width: full width of field box, in arcseconds
    err_cut: require models to have fitting parameters lower than this value,
              return None instead of SED in final_seds if the model is a poor fit.
    redden: T/F, include galactic reddening in model fit?
    return_model: If True, returns all modeled magnitudes.
                  If False, returns modeled magnitudes supplemented
                  by observed.
    '''
    ra, dec = field_center #in decimal degrees
    mass, sdss, usnob = query_all(ra, dec, boxsize=field_width)
    
    object_mags = []
    modes = []
    object_coords = []
    if mass != None:
        # match sdss, usnob objects to 2mass objects
        if sdss != None:
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
    
    # now fit a model to each object, and construct the final SED,
    #  filling in missing observations with synthetic photometry.
    final_seds, out_coords, out_modes, out_errs, models_used = [], [], [], [], []
    for i, obs in enumerate(object_mags):
        mode = modes[i]
        # the masks show, in order, which bands are included
        #  order: u,g,r,i,z, y, B,R, J,H,K
        if mode == 0: # sdss+2mass
            mask = [1,1,1,1,1, 0, 0,0, 1,1,1]
        elif mode == 1: # usnob+2mass
            mask = [0,0,0,0,0, 0, 1,1, 1,1,1]
        mask = np.array(mask).astype(bool)
        
        if redden:
            reddening = get_reddening( ra,dec, ALL_FILTERS )
            # de-redden the observations before comparing
            #  to the models
            obs[::2] -= reddening[mask]
        
        model, C, T, err = choose_model( obs, mask )
        
        if err > err_cut: #apply quality-of-fit cut
            continue
        else:
            if redden:
                # re-redden the model and observations
                obs[::2] += reddening[mask]
                model += reddening
                
            sed = np.empty(11)
            # keep the real observations
            sed[mask] = obs[::2]
            # fill in rest with modeled magnitudes
            sed[~mask] = model[~mask]
            final_seds.append( sed )
            out_coords.append( object_coords[i] )
            out_modes.append( mode )
            out_errs.append( err )
            models_used.append( (T,C) ) #( index or temp, offset )
    
    if return_models:
        return out_coords, final_seds, out_modes, out_errs, models_used
    else:
        return out_coords, final_seds, out_modes, out_errs


def _split_field( field_center, field_width, max_size, object_coords=None ):
    '''
    split a single field into many smaller chunks, to save time matching sources.
    
    field_center: (ra, dec) in decimal degrees
    field_width: each side of box, in arcseconds (RA_width, DEC_width)
    max_size: maximum size of tile, in arcseconds
    
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


def _find_field( star_coords, extend=.0015 ):
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


############################################
# MODEL FITTING FUNCTIONS
############################################

def get_reddening( ra, dec, filters, dust_map=MAP_DICT ):
    '''
    returns reddening for filters at ra,dec as 
    defined by the dust_map_dict
    
    example:
     de_reddened_mags = mags - reddening( ra,dec,filters )
    '''
    # coordinate-to-pixel mapping from the dust map fits header
    X_pix = lambda l,b,pole: np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.cos(l)*dust_map[pole][2]
    Y_pix = lambda l,b,pole: -dust_map[pole][1]*np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.sin(l)*dust_map[pole][2]
    
    # get galactic coordinates with eq_gal, which does everything in radians
    ra_rad = ra*(np.pi/180.)
    dec_rad = dec*(np.pi/180.)
    l,b = eq_gal( 2000., ra_rad, dec_rad )
    if b>0:
        pole = 'N'
    else:
        pole = 'S'
    
    # get E(B-V) for these coordinates
    X = int(round( X_pix(l,b,pole) ))
    Y = int(round( Y_pix(l,b,pole) ))
    EBV = dust_map[pole][0][X,Y]
    
    # get reddening value for each filter in filters
    reddening = []
    for filt in filters:
        reddening.append( dust_map['Filters'][filt]*EBV )
    
    return np.array( reddening )


def _error_C(C, model, obs, weights ):
    '''
    C: number, a constant
    model, obs, weights: numpy arrays
    returns: sum squared errors
    '''
    nm = model*C
    return np.sum( weights*(nm-obs)**2 )

def mag2flam( magnitudes, bands ):
    '''
    Converts magnitudes to flam (erg/s/cm^2/A).
    '''
    flam = np.empty_like(magnitudes)
    for i,b in enumerate(bands):
        f0 = FILTER_PARAMS[b][1]
        flam[i] = f0*10**(-.4*magnitudes[i])
    return flam

def flam2mag( flams, bands ):
    '''
    Converts flam back to magnitudes.
    '''
    mags = np.empty_like(flams)
    for i,b in enumerate(bands):
        f0 = FILTER_PARAMS[b][1]
        mags[i] = -2.5*np.log10( flams[i]/f0 )
    return mags

def choose_model( obs, mask, models=MODELS ):
    '''
    Find and return the best model for obs.
    Do this by fitting to all magnitudes weighted by error.
    
    Returns: model, C, temperature, quality_parameter
    
    obs: an array of the observations and errors (in mags), in
           order as defined by the mode key (see below)
    mask: defines what colors to use in the fit, i.e. what observations exist
           in order: [u,g,r,i,z,y,B,R,J,H,K]
    models: an array of modeled SEDs, where 0th entry is temperature or index
             of the model, and the rest are brightnesses (in flam)
    '''
    # mask is an array used to choose which modeled magnitudes
    #  correspond to the included observations.
    #  Note: I prepend a zero to ignore the temperature/index of the model
    #   (saved in position 0 for all models)
    model_mask = np.hstack( (np.array([False]), np.array(mask).astype(bool)) )
    
    flams = mag2flam( obs[::2], ALL_FILTERS[mask] )
    weights = 1./obs[1::2]
    
    # Go through all models and choose the one with the most similar SED
    #  Keep track of the sum_squared error, as returned by _error_C()
    sum_sqrs, Cs = [], []
    # inside the loop, we'll use the below to estimate a value for C
    f0 = FILTER_PARAMS[ ALL_FILTERS[mask][0] ][1]
    for model in models[1:]:
        # estimate C
        f_mod = model[1:][mask][-2]
        C_est = (f0/f_mod)*10**(-.4*obs[::2][-2]) # using Mag = -2.5Log[ (C*f_model)/f0 ] for Hband
        res = fmin_bfgs( _error_C, C_est, args=(model[1:][mask], flams, weights), full_output=True, disp=False )
        Cs.append(res[0][0])
        sum_sqrs.append(res[1])
        import pdb; pdb.set_trace()
        
    i_best = np.argmin(sum_sqrs)
    C = Cs[i_best]
    best_model = models[1:][i_best]
    best_model_mags = flam2mag( best_model[1:]*C, ALL_FILTERS )
    # return all magnitudes for best model, the offset C, the temperature,
    #  and a quality metric for the best fit.
    #  The quality metric is the average error (in mags) between the best model and the observations
    sum_sqr_mag_err = _error_C( C, best_model_mags[mask], obs[::2], np.ones_like(obs[::2]) )
    return (best_model_mags, C, best_model[0], (sum_sqr_mag_err/len(obs[::2]))**.5 )


############################################
# MAIN FUNCTIONS
############################################

def catalog( field_center, field_width, object_coords=None, 
               redden=False, savefile=None, max_size=1800., return_models=False):
    '''
    Main cataloging function, this produces a catalog of all objects found in field.
    Requires records in 2MASS + (SDSS and/or USNOB1).
    
    field_center: (ra, dec) in decimal degrees
    field_width: full width of field box, in arcseconds
    object_coords: found objects we will compare to; if not none,
               the script will not catalog fields in which there are no objects
    redden: boolean; account for galactic reddening
    savefile: optional; saves to specified file if present, otherwise returns answer
    return_models: if True, returns the index (or temperature) of model used for each source
    
    NOTE: if requested field_width is greater than max_size (in arcsec),
          this splits up the request into 900-arcsec chunks, to save time.
    '''
    if max(field_width) > max_size:
        # split field up into smaller chunks, to run more quickly
        centers, tile_width = _split_field( field_center, field_width, max_size, object_coords=object_coords )
        
        # go through each tile and accumulate the results:
        object_coords, final_seds, modes, errors, models = [],[],[],[],[]
        for i,center in enumerate(centers):
            print i,'of',len(centers)
            oc, fs, ms, ers, mods = produce_catalog( center, tile_width, redden=redden, return_models=True )
            object_coords += oc
            final_seds += fs
            modes += ms
            errors += ers
            models += mods
    else:
        object_coords, final_seds, modes, errors, models = \
                         produce_catalog( field_center, max(field_width), redden=redden, return_models=True )
    
    # Done! Save to file, or return SEDs and coordinates
    if savefile:
        f_format = lambda x: str(round(x, 5)) # a quick function to format the output
        fff = open(savefile,'w')
        fff.write('# Produced by get_SEDs.py \n# Catalog of objects in field of ' +
                  'size (%.f, %.f) (RA,DEC in arcsec) centered at (%.4f, %.4f).\n' %( field_width[0], field_width[1], field_center[0], field_center[1]) +
                  '# modes: 0 -> SDSS+2MASS; 1 -> USNOB1+2MASS\n'
                  '# RA\tDEC\t' + '\t'.join(ALL_FILTERS) + '\tmode\terror\n')
        for i,row in enumerate(final_seds):
            fff.write( '\t'.join( map(str, object_coords[i]) ) + '\t' + '\t'.join( map(f_format, row) ) + '\t{}\t{}\n'.format(modes[i], f_format(errors[i])) )
        fff.close()
    else:
        if return_models:
            return np.array(object_coords), np.array(final_seds), models
        else:
            return np.array(object_coords), np.array(final_seds)
    


def calc_zeropoint( input_coords, catalog_coords, input_mags, catalog_mags, sigma_cut=2. ):
    '''
    Calculate the zeropoint for a set of input stars and set of catalog stars.
    
    input_coords: a 2D array of [ [RA, DEC], [..., ...] ... ]
    catalog_coords: similar array for objects created with catalog()
    input_mags: a 1D array of instrumental magnitudes
    catalog_mags: a similar array of true magnitudes as created with catalog()
    sigma_cut: trim values beyond SC*sigma from mean
    
    Returns: zeropoint (mags) and an estimate of quality
    
    TO DO: split these fields before matching.  See how it is done in produce_catalog.
    '''
    matches = identify_matches( input_coords, catalog_coords )
    matched_inputs = [input_mags[i] for i in range(len(matches)) if matches[i] != None ]
    matched_catalogs = [catalog_mags[ matches[i] ] for i in range(len(matches)) if matches[i] != None ]
    
    zp_estimates = []
    for i,inst_mag in enumerate(matched_inputs):
        zp_estimates.append( matched_catalogs[i] - inst_mag )
    zp = np.array(zp_estimates)
    zp_cut = zp[ np.abs(zp-np.mean(zp)) < 2*np.std(zp) ]
    return np.mean(zp_cut)


def zeropoint( input_file ):
    '''
    Calculate zeropoint for stars in an input text file.
    
    Maybe run a crontab to delete these files in datadir daily?
    
    Expects input file of format:
     IDSTRING_band_XXX.txt  (ascii, np.loadtxt-readable)
    '''
    # internal definitions
    data_dir = 'data/'
    sdss_map = { 'u':2, 'g':4, 'r':6, 'i':8, 'z':10 }  # index of band in query_sdss() response
    catalog_map = { 'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5 }  # index of band in catalog() response
    
    # load the data and pick a field
    in_data = np.loadtxt( input_file )
    input_coords = in_data[:, :2]
    input_mags = in_data[:, 2]
    field_center, field_width = _find_field( input_coords )
    idstring, band = input_file.split('_')[:2]
    
    # quick test whether field is in SDSS:
    ra,dec = field_center
    sdss = query_sdss( ra,dec, boxsize=50., trim_mag=30. )
    if sdss != None:
        in_sdss = True
    else:
        in_sdss = False
    
    # if field is in SDSS, see whether we've already queried and saved the objects
    fmt = lambda x: str( round(x,2) )
    if in_sdss:
        if band.lower() == 'y':
            if isfile( data_dir+idstring+'_catalog_coords.npy' ):
                catalog_coords = np.load( data_dir+idstring+'_catalog_coords.npy' )
                catalog_seds = np.load( data_dir+idstring+'_catalog_seds.npy')
            else:
                catalog_coords, catalog_seds = catalog( field_center, field_width, object_coords=input_coords )
                np.save( data_dir+idstring+'_catalog_coords.npy', catalog_coords )
                np.save( data_dir+idstring+'_catalog_seds.npy', catalog_seds )
            catalog_mags = catalog_seds[:, catalog_map['y'] ]
            
        else:
            if isfile( data_dir+idstring+'_sdss_query.npy' ):
                sdss = np.load( data_dir+idstring+'_sdss_query.npy' )
            else:
                sdss = query_sdss( ra, dec, boxsize=field_width )
                np.save( data_dir+idstring+'_sdss_query.npy', sdss )
            catalog_coords = sdss[:,:2]
            catalog_mags = sdss[:, sdss_map[band.lower()] ]
    
    else: # hit here if field is not in SDSS, but check to see whether we've already built the catalog
        if isfile( data_dir+idstring+'_catalog_coords.npy' ):
            catalog_coords = np.load( data_dir+idstring+'_catalog_coords.npy' )
            catalog_seds = np.load( data_dir+idstring+'_catalog_seds.npy' )
        else:
            catalog_coords, catalog_seds = catalog( field_center, field_width, object_coords=input_coords )
            np.save( data_dir+idstring+'_catalog_coords.npy', catalog_coords )
            np.save( data_dir+idstring+'_catalog_seds.npy', catalog_seds )
        catalog_mags = catalog_seds[:, catalog_map[band.lower()] ]
    
    zp = calc_zeropoint( input_coords, catalog_coords, input_mags, catalog_mags )
    return zp

