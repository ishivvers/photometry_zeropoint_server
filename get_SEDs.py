'''
Library to produce a catalog of fully-populated SEDs for
any arbitrary location on the sky, using a combination 
of online catalogs (USNOB1, 2MASS, SDSS) and synthetic photometry.

Requires:
- all_models.npy, a file produced by assemble_models.py
- Schlegel et al. dust extiction maps in a folder named dust_maps

'''


############################################
# IMPORTS
############################################
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from scipy.optimize import minimize
import pyfits, ephem
from os.path import isfile
from threading import Thread

try:
    MODELS = np.load( open('all_models.npy','r') )
except:
    raise IOError('cannot find models file')

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

ALL_FILTERS = ['u','g','r','i','z','y','B','R','J','H','K']

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
    match_radius_squared = (2.778e-4*match_radius)**2 # convert arcseconds into degrees
    matches = []
    #matched = np.zeros(len(found_stars)).astype(bool) # mask to track which stars have already been matched
    for star in queried_stars:
        # calculate the distance to each located star
        diffs_squared = [((star[0]-other[0])**2 + (star[1]-other[1])**2) for other in found_stars] #[~matched]]
        if min(diffs_squared) < match_radius_squared:
            i_best = np.argmin(diffs_squared)
            matches.append(i_best)
            #matched[i_best] = True
        else:
            matches.append(None) 
    return matches


def produce_catalog( field_center, field_width, err_cut=2., redden=True, return_model=False ):
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
    final_seds, out_coords, out_modes, out_errs = [], [], [], []
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
            reddening = _get_reddening( ra,dec, ALL_FILTERS )
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
                
            if return_model:
                sed = model
            else:
                sed = np.empty(11)
                # keep the real observations
                sed[mask] = obs[::2]
                # fill in rest with modeled magnitudes
                sed[~mask] = model[~mask]
            final_seds.append( sed )
            out_coords.append( object_coords[i] )
            out_modes.append( mode )
            out_errs.append( err )
    
    return out_coords, final_seds, out_modes, out_errs


def _split_field( field_center, field_width, max_size ):
    '''
    split a single field into many smaller chunks, to save time matching sources.
    
    field_center: (ra, dec) in decimal degrees
    field_width: one side of box, in arcseconds
    max_size: maximum size of tile, in arcseconds
    
    Returns: list of new field centers, width of new fields (in arcseconds)
    '''
    a2d = 2.778e-4 #conversion between arcseconds & degrees
    
    n_tile = np.ceil(field_width/max_size).astype(int) # number of tiles needed in each dimension
    w_tile = (field_width/n_tile)*a2d # width of each tile in degrees
    
    ra,dec = field_center
    centers = []
    rr = ra - a2d*field_width/2. + w_tile/2.  # the beginning positions of the tiling
    for i in range(n_tile):
        dd = dec - a2d*field_width/2. + w_tile/2.
        for j in range(n_tile):
            centers.append( (rr, dd) )
            dd += w_tile
        rr += w_tile
    
    return centers, w_tile/a2d


def _find_field( star_coords, extend=.0015 ):
    '''
    determine the best field for a list of star coordinates,
    so we can perform only a single query.
    
    star_coords: a list of star coordinates, in decimal degrees
    extend: the buffer beyond the requested coordinates to add to the field
      (also in decimal degrees)
    
    returns: (coordinates of center in decimal degrees), width_of_box (in arcseconds)
    '''
    ras = star_coords[:,0]
    decs = star_coords[:,1]
    width_ra = (max(ras)-min(ras) + 2*extend)
    center_ra = np.mean(ras)
    width_dec = (max(decs)-min(decs) + 2*extend)
    center_dec = np.mean(decs)
    return (center_ra, center_dec), max(width_ra, width_dec)*3600.


############################################
# MODEL FITTING FUNCTIONS
############################################

def _get_reddening( ra, dec, filters, dust_map=MAP_DICT ):
    '''
    returns reddening for filters at ra,dec as 
    defined by the dust_map_dict
    
    example:
     de_reddened_mags = mags - reddening( ra,dec,filters )
    '''
    # coordinate-to-pixel mapping from the dust map fits header
    X_pix = lambda l,b,pole: np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.cos(l)*dust_map[pole][2]
    Y_pix = lambda l,b,pole: -dust_map[pole][1]*np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.sin(l)*dust_map[pole][2]
    
    # get galactic coordinates with pyephem, which does everything in radians
    ra_rad = ra*(np.pi/180.)
    dec_rad = dec*(np.pi/180.)
    coords_EQ = ephem.Equatorial(ra_rad,dec_rad)
    coords_GA = ephem.Galactic(coords_EQ)
    l = coords_GA.lon
    b = coords_GA.lat
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


def _error_C(C, model, obs):
    '''
    C: number, a constant akin to distance modulus
    model: array-like, model mags
    real: array-like, observed mags
    
    returns: sum squared errors
    '''
    nm = model+C
    return np.sum( (nm-obs)**2 )


def _error_C_reddening(pars, model, reddening, obs):
    '''
    C: number, a constant akin to distance modulus
    R: scalar by which to add reddening
    model: array-like, model mags
    reddening: the reddening array for model (mags)
    real: array-like, observed mags
    
    returns: sum squared errors
    '''
    C,R = pars
    nm = model + C + R*reddening
    return np.sum( (nm-obs)**2 )


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
    zerod_mags = mags - min(mags) # recenter to compare to models
    weights = 1./obs[1::2]
    #weights = np.ones(len(obs[1::2])) #try flat weights
    
    # Go through all models and choose the one with the most similar SED
    #  Keep track of the sum_squared error
    sum_sqrs, Cs = [], []
    for model in models[1:]:
        res = minimize( _error_C, 0., args=(model[mask], zerod_mags) )
        Cs.append(res.x[0])
        sum_sqrs.append(res.values()[4])
    i_best = np.argmin(sum_sqrs)
    best_model = models[1:][ i_best ] 
    # now add back in the zeropoint to get a model for the non-zerod observations
    C = Cs[i_best] + min(mags)
    # return all magnitudes for best model, the temperature, and a quality metric for the best fit
    return (best_model[1:] + C, best_model[0], sum_sqrs[i_best])


def choose_model_reddening( obs, mask, reddening, models=MODELS ):
    '''
    Find and return the best model for obs.
    Do this by fitting to all magnitudes weighted by error.
    
    Returns: model, temperature, quality_parameter
    
    obs: an array of the observed magnitudes and errors, in
           order as defined by the mode key (see below)
    mask: defines what colors to use in the fit, i.e. what observations exist
           in order: [u,g,r,i,z,y,B,R,J,H,K]
    reddening: an array of reddening corrections for these coordinates (mags)
    models: an array of modeled SEDs, where 0th entry is temperature
             of the model, and the rest are magnitudes
    '''
    # mask is an array used to choose which modeled magnitudes
    #  correspond to the included observations.
    #  Note: I append a leading zero to ignore the temperature of the model
    #   (saved in position 0 for all models)
    mask = np.hstack( (np.array([False]), np.array(mask).astype(bool)) )
    
    mags = obs[::2]
    zerod_mags = mags - min(mags) # recenter to compare to models
    weights = 1./obs[1::2]
    #weights = np.ones(len(obs[1::2])) #try flat weights
    
    # Go through all models and choose the one with the most similar SED
    #  Keep track of the sum_squared error
    sum_sqrs, Cs, Rs = [], [], []
    for model in models[1:]:
        res = minimize( _error_C_reddening, (0.,1.), args=(model[mask], reddening[mask[1:]], zerod_mags),
                        method='TNC', bounds=( (None,None), (0., 1.2)) )
        Cs.append(res.x[0])
        Rs.append(res.x[1])
        sum_sqrs.append(res.values()[3])
    i_best = np.argmin(sum_sqrs)
    best_model = models[1:][ i_best ] 
    # now add back in the zeropoint to get a model for the non-zerod observations
    C = Cs[i_best] + min(mags)
    R = Rs[i_best]
    # return all magnitudes for best model, the temperature, and a quality metric for the best fit
    return (best_model[1:] + C + R*reddening, best_model[0], sum_sqrs[i_best])


############################################
# MAIN FUNCTIONS
############################################

def catalog( field_center, field_width, redden=False, savefile=None, max_size=1800.):
    '''
    Main cataloging function, this produces a catalog of all objects found in field.
    Requires records in 2MASS + (SDSS and/or USNOB1).
    
    field_center: (ra, dec) in decimal degrees
    field_width: full width of field box, in arcseconds
    redden: boolean; account for galactic reddening
    savefile: optional; saves to specified file if present, otherwise returns answer
    
    NOTE: if requested field_width is greater than max_size (in arcsec),
          this splits up the request into 900-arcsec chunks, to save time.
    '''
    if field_width > max_size:
        # split field up into smaller chunks, to run more quickly
        centers, tile_width = _split_field( field_center, field_width, max_size )
        
        # go through each tile and accumulate the results:
        object_coords, final_seds, modes, errors = [],[],[],[]
        for i,center in enumerate(centers):
            oc, fs, ms, ers = produce_catalog( center, tile_width, redden=redden )
            object_coords += oc
            final_seds += fs
            modes += ms
            errors += ers
    else:
        object_coords, final_seds, modes, errors = produce_catalog( field_center, field_width, redden=redden )
    
    # Done! Save to file, or return SEDs and coordinates
    if savefile:
        format = lambda x: str(round(x, 3)) # a quick function to format the output
        fff = open(savefile,'w')
        fff.write('# Produced by get_SEDs.py \n# Catalog of objects in field of ' +
                  'size {} (arcsec) centered at {}.\n'.format( field_width, field_center) +
                  '# modes: 0 -> SDSS+2MASS; 1 -> USNOB1+2MASS\n'
                  '# RA\tDEC\t' + '\t'.join(ALL_FILTERS) + '\tmode\tdimensionless_error\n')
        for i,row in enumerate(final_seds):
            fff.write( '\t'.join( map(str, object_coords[i]) ) + '\t' + '\t'.join( map(format, row) ) + '\t{}\t{}\n'.format(modes[i], format(errors[i])) )
        fff.close()
    else:
        return np.array(object_coords), np.array(final_seds)
    


def _zeropoint( input_coords, catalog_coords, input_mags, catalog_mags, sigma_cut=2. ):
    '''
    Calculate the zeropoint for a set of input stars and set of catalog stars.
    
    input_coords: a 2D array of [ [RA, DEC], [..., ...] ... ]
    catalog_coords: similar array for objects created with catalog()
    input_mags: a 1D array of instrumental magnitudes
    catalog_mags: a similar array of true magnitudes as created with catalog()
    sigma_cut: trim values beyond SC*sigma from mean
    
    Returns: zeropoint (mags) and an estimate of quality
    '''
    matches = identify_matches( input_coords, catalog_coords )
    matched_inputs = [input_mags[i] for i in range(len(matches)) if matches[i] != None ]
    matched_catalogs = [catalog_mags[ matches[i] ] for i in range(len(matches)) if matches[i] != None ]
    
    zp_estimates = []
    for i,inst_mag in enumerate(matched_inputs):
        zp_estimates.append( matched_catalogs[i] - inst_mag )
    zp = np.array(zp_estimates)
    zp_cut = e[ np.abs(zp-np.mean(zp)) < 2*np.std(zp) ]
    
    return np.mean(zp_cut)


def zeropoint( input_file, band ):
    '''
    Calculate zeropoint for stars in an input text file.
    '''
    # internal definitions
    data_dir = 'data'
    sdss_map = { 'u':2, 'g':4, 'r':6, 'i':8, 'z':10 }  # index of band in query_sdss() response
    catalog_map = { 'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5 }  # index of band in catalog() response
    
    # load the data and pick a field
    in_data = np.loadtxt( input_file )
    input_coords = in_data[:, :2]
    input_mags = in_data[:, 2]
    field_center, field_width = _find_field( input_coords )
    
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
            if isfile( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ) ):
                catalog_coords = np.load( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ) )
                catalog_seds = np.load( data_dir+'/{}_{}_catalog_seds.np'.format( fmt(ra), fmt(dec) ) )
            else:
                catalog_coords, catalog_seds = catalog( field_center, field_width )
                np.save( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ), catalog_coords
                np.save( data_dir+'/{}_{}_catalog_seds.np'.format( fmt(ra), fmt(dec) ), catalog_seds )
            catalog_mags = catalog_seds[:, catalog_map['y'] ]
            
        else:
            if isfile( data_dir+'/{}_{}_sdss_query.np'.format( fmt(ra), fmt(dec) ) ):
                sdss = np.load( data_dir+'/{}_{}_sdss_query.np'.format( fmt(ra), fmt(dec) ) )
            else:
                sdss = query_sdss( ra, dec, boxsize=field_width )
                np.save( data_dir+'/{}_{}_sdss_query.np'.format( fmt(ra), fmt(dec) ), sdss )
            catalog_coords = sdss[:,:2]
            catalog_mags = sdss[:, sdss_map[band.lower()] ]
    
    else: # hit here if field is not in SDSS, but check to see whether we've already built the catalog
        if isfile( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ) ):
            catalog_coords = np.load( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ) )
            catalog_seds = np.load( data_dir+'/{}_{}_catalog_seds.np'.format( fmt(ra), fmt(dec) ) )
        else:
            catalog_coords, catalog_seds = catalog( field_center, field_width )
            np.save( data_dir+'/{}_{}_catalog_coords.np'.format( fmt(ra), fmt(dec) ), catalog_coords
            np.save( data_dir+'/{}_{}_catalog_seds.np'.format( fmt(ra), fmt(dec) ), catalog_seds )
        catalog_mags = catalog_seds[:, catalog_map[band.lower()] ]
    
    zp = _zeropoint( input_coords, catalog_coords, input_mags, catalog_mags )
    return zp

