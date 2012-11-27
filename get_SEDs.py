'''
Library to produce a catalog of fully-populated SEDs for
any arbitrary location on the sky, using a combination 
of online catalogs (USNOB1, 2MASS, SDSS) and synthetic photometry.

Requires:
- all_models.npy, a file produced by assemble_models.py
- Schlegel et al. dust extiction maps in a folder named dust_maps

'''


############################################
# imports
############################################
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from scipy.optimize import minimize
from sys import exit
import pyfits, ephem
from threading import Thread

try:
    MODELS = np.load( open('all_models.npy','r') )
except:
    print 'cannot find models file'
    exit()

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
    print 'cannot find/open dust maps'
    exit()
# Include the E(B-V)-to-mag mapping for all filters,
#  from http://www.astro.princeton.edu/~schlegel/dust/data/filter.txt
#  except for y-band, which is from Tonry et al. 2012
#  NOTE: using values for UKIRT JHK filters instead of 2MASS - close enough!
MAP_DICT['Filters'] = {'u':5.155, 'g':3.793, 'r':2.751, 'i':2.086,
                       'z':1.479, 'y':1.251, 'B':4.315, 'R':2.673,
                       'J':0.902, 'H':0.576, 'K':0.367}


ALL_FILTERS = ['u','g','r','i','z','y','B','R','J','H','K']

############################################
# CATALOG INTERFACE FUNCTIONS
############################################

def parse_sdss( s ):
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


def query_sdss( ra, dec, boxsize=10., container=None, cont_index=1 ):
    '''
    Query sdss8 server for sources found in a box of width boxsize (arcsecs)
    around ra,dec.
    Returns an array (if objects found) or None (if not)

    ra,dec: coordinates in decimal degrees
    boxsize: width of box in which to query
    '''
    # search for SDSS objects around coordinates with
    #  defined box size, return only basic parameters, and
    #  sort by distance from coordinates, and return a maximum
    #  of 10000 sources (i.e. return everything)
    request = 'findsdss8 -c "{} {}" -bs {} -e0 -sr -m 10000'.format( ra, dec, boxsize )
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    sdss_objects = parse_sdss(o)
    if len(sdss_objects) == 0:
        # no matches
        output = None
    else:
        output = sdss_objects
    if container == None:
        return output
    else:
        container[ cont_index ] = output


def parse_2mass( s ):
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
    mass_objects = parse_2mass(o)
    if len(mass_objects) == 0:
        # no matches
        output = None
    else:
        output = mass_objects
    if container == None:
        return output
    else:
        container[ cont_index ] = output


def parse_usnob1( s ):
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
    usnob1_objects = parse_usnob1(o)
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


def identify_matches( queried_stars, found_stars, match_radius=.5 ):
    '''
    Match queried stars to found stars.
    
    Returns list of indices in found_stars for corresponding
     star in queried_stars, or None if no match found.

    queried_stars: an array of coordinates of stars to match
    found_stars: an array of coordinates of located stars
    match_radius: the maximum offset between queried and found star to 
      call a match, in arcseconds
    '''
    match_radius = 2.778e-4*match_radius # convert arcseconds into degrees
    matches = []
    for star in queried_stars:
        # calculate the distance to each located star
        diffs = [((star[0]-other[0])**2 + (star[1]-other[1])**2)**.5 for other in found_stars]
        if min(diffs) < match_radius:
            matches.append(np.argmin(diffs))
        else:
            matches.append(None) 
    return matches


def produce_catalog( field_center, field_width, redden=True, return_model=False ):
    '''
    Create a catalog of all objects found in field.
    Requires records in 2MASS + (SDSS and/or USNOB1).

    field_center: (ra, dec) in decimal degrees
    field_width: full width of field box, in arcseconds
    redden: T/F, include galactic reddening in model fit?
    return_model: If True, returns all modeled magnitudes.
                  If False, returns modeled magnitudes supplemented
                  by observed.
    '''
    ra, dec = field_center #in decimal degrees
    mass, sdss, usnob = query_all(ra, dec, boxsize=field_width)

    # before matching, prune the SDSS catalog down to only
    #  the bright sources (saves time)
    if sdss != None:
        mask = []
        for obj in sdss:
            if np.mean( obj[2::2] ) > 20:
                # if the average SDSS magnitude is greater than 20, ignore this object
                mask.append(0)
            else:
                mask.append(1)
        mask = np.array(mask).astype(bool)
        sdss = sdss[mask]

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
    final_seds = []
    for i, obs in enumerate(object_mags):
        mode = modes[i]
        # the masks show, in order, which bands are included
        #  order: u,g,r,i,z,y,B,R,J,H,K
        if mode == 0: # sdss+2mass
            mask = [1,1,1,1,1,0,0,0,1,1,1]
        elif mode == 1: # usnob+2mass
            mask = [0,0,0,0,0,0,1,1,1,1,1]
        mask = np.array(mask).astype(bool)

        if redden:
            reddening = get_reddening( ra,dec, ALL_FILTERS )
            # de-redden the observations before comparing
            #  to the models
            obs[::2] -= reddening[mask]

        model, T = choose_model( obs, mask )

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
        final_seds.append(sed)

    return object_coords, final_seds, modes


def split_field( field_center, field_width, max_size ):
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


def error_C(C, model, obs):
    '''
    C: number, a constant akin to distance modulus
    model: array-like, model mags
    real: array-like, observed mags

    returns: sum squared errors
    '''
    nm = model+C
    return np.sum( (nm-obs)**2 )


def choose_model( obs, mask, models=MODELS ):
    '''
    Find and return the best model for obs.
    Do this by fitting to all magnitudes weighted by error.
    
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
        res = minimize( error_C, 0., args=(model[mask], zerod_mags) )
        C = res.x[0]
        sum_sqrs.append(np.sum( weights*(zerod_mags - model[mask]+C)**2 ))
    best_model = models[1:][ np.argmin(sum_sqrs) ]
    # now determine the C that best fits the non-zerod observations
    guess = obs[-2] # the J-Mag, since all models are centered at J=0
    res = minimize( error_C, guess, args=(best_model[mask], mags))
    C = res.x[0]
    # return all magnitudes for best model and the temperature
    return (best_model[1:] + C, best_model[0])


############################################
# MAIN FUNCTION
############################################

def catalog( field_center, field_width, redden=True, savefile=None, max_size=1800.):
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
        centers, tile_width = split_field( field_center, field_width, max_size )

        # go through each tile and accumulate the results:
        object_coords, final_seds, modes = [],[],[]
        for i,center in enumerate(centers):
            oc, fs, ms = produce_catalog( center, tile_width, redden=redden )
            object_coords += oc
            final_seds += fs
            modes += ms
    else:
        object_coords, final_seds, modes = produce_catalog( field_center, field_width, redden=redden )

    # Done! Save to file, or return SEDs and coordinates
    if savefile:
        format = lambda x: str(round(x, 3)) # a quick function to format the output
        fff = open(savefile,'w')
        fff.write('# Produced by get_SEDs.py \n# Catalog of objects in field of ' +
                  'size {} (arcsec) centered at {}.\n'.format( field_width, field_center) +
                  '# modes: 0 -> SDSS+2MASS; 1 -> USNOB1+2MASS\n'
                  '# RA\tDEC\t' + '\t'.join(ALL_FILTERS) + '\tmode\n')
        for i,row in enumerate(final_seds):
            fff.write( '\t'.join( map(str, object_coords[i]) ) + '\t' + '\t'.join( map(format, row) ) + '\t{}\n'.format(modes[i]) )
        fff.close()
    else:
        return np.array(object_coords), np.array(final_seds)


############################################
# TEST FUNCTIONS
############################################


def test_SDSS_errors( ra, dec, band_name='z', redden=True, size=900. ):
    '''
    Test the error accrued for all sources in a field when estimating 
     SDSS z-band photometry from all modes.  y-band photometry is 
     expected to be similar or better than this.
    
    Run this on any field within the SDSS footprint.

    ra,dec: coordinates in decimal degrees
    band: SDSS passband to derive errors for
    redden: if True, account for galactic reddening when modeling photometry
    size: size of field to query around ra, dec (arseconds)
    '''
    if band_name == 'u':
        sdss_mask = [0,0,1,1,1,1,1,1,1,1]
        i_band = 0
    elif band_name == 'g':
        sdss_mask = [1,1,0,0,1,1,1,1,1,1]
        i_band = 1
    elif band_name == 'r':
        sdss_mask = [1,1,1,1,0,0,1,1,1,1]
        i_band = 2
    elif band_name == 'i':
        sdss_mask = [1,1,1,1,1,1,0,0,1,1]
        i_band = 3
    elif band_name == 'z':
        sdss_mask = [1,1,1,1,1,1,1,1,0,0]
        i_band = 4
    else:
        print 'Incorrect band keyword!'
        exit()
    sdss_mask = np.array(sdss_mask).astype(bool)

    mass, sdss, usnob = query_all(ra, dec, boxsize=size)

    ##### the code between 5 hashes is modified from produce_catalog, above #####
    object_mags = []
    band_mags = []
    modes = []
    object_coords = []

    # match sdss, usnob objects to 2mass objects
    if sdss != None:
        # before matching, prune the SDSS catalog down to only
        #  the bright sources (saves time)
        mask = []
        for obj in sdss:
            if np.mean( obj[2::2] ) > 20:
                # if the average SDSS magnitude is greater than 20, ignore this object
                mask.append(0)
            else:
                mask.append(1)
        mask = np.array(mask).astype(bool)
        sdss = sdss[mask]
        
        mass_matches = identify_matches( sdss[:,:2], mass[:,:2] )
        usnob_matches = identify_matches( sdss[:,:2], usnob[:,:2] )
    else:
        print 'Must run this for coordinates in SDSS footprint!'
        exit()

    # Go through SDSS objects and assemble a catalog
    #  of all objects present in multiple catalogs
    for i,obj in enumerate(sdss):
        if mass_matches[i] != None and usnob_matches[i] !=None:
            i_mass = mass_matches[i]
            i_usnob = usnob_matches[i]
            # fit to SDSS+2MASS
            obs = np.hstack( (obj[2:][sdss_mask], mass[i_mass][2:]) )
            band = obj[2:][~sdss_mask][0] # keep track of the true band mag
            object_mags.append( obs )
            modes.append( 0 )
            band_mags.append( band )
            # also fit to USNOB+2MASS
            obs = np.hstack( (usnob[i_usnob][2:], mass[i_mass][2:]) )
            object_mags.append( obs )
            modes.append( 1 )
            band_mags.append( band )
        elif mass_matches[i] != None and usnob_matches[i] == None:
            i_mass = mass_matches[i]
            # fit to SDSS+2MASS
            obs = np.hstack( (obj[2:][sdss_mask], mass[i_mass][2:]) )
            band = obj[2:][~sdss_mask][0]
            object_mags.append( obs )
            modes.append( 0 )
            band_mags.append( band )

    # now fit a model to each object, and construct the final SED,
    #  without including the z-band.  Determine errors between
    #  predicted z-band and actual. Build an array of sample SEDs of
    #  the first 9 sources for each type of fit.
    # Modes are defined as:
    #  0 -> SDSS+2MASS; 1 -> USNOB+2MASS
    pltsize = 3 # adjust this parameter to show more/fewer SEDs in figures 1,2
    f_0, axs_0 = plt.subplots( pltsize, pltsize, sharex=True, figsize=(15,10))
    f_1, axs_1 = plt.subplots( pltsize, pltsize, sharex=True, figsize=(15,10))
    axs_0 = axs_0.flatten()
    axs_1 = axs_1.flatten()
    i_ax0, i_ax1 = 0,0
    
    errors_0, errors_1 = [],[]
    for i, obs in enumerate(object_mags):
        mode = modes[i]
        if mode == 0:
            mask = list(sdss_mask[::2])+[0,0,0,1,1,1]
        elif mode == 1:
            mask = [0,0,0,0,0,0,1,1,1,1,1]
        mask = np.array(mask).astype(bool)

        if redden:
            reddening = get_reddening( ra,dec, ALL_FILTERS )
            # de-redden the observations before comparing
            #  to the models
            obs[::2] -= reddening[mask]

        model, T = choose_model( obs, mask )
        if redden:
            # re-redden the model and observations
            obs[::2] += reddening[mask]
            model += reddening

        # compare calculated z-mag to observed
        true_band = band_mags[i]
        guess_band = model[i_band]
        error = true_band - guess_band
        if mode == 0:
            errors_0.append(error)
        elif mode == 1:
            errors_1.append(error)
            
        # plot up the SEDs themselves
        ax = None
        if mode == 0 and (i_ax0 < len(axs_0)):
            i_ax = i_ax0
            ax = axs_0[i_ax]
            if i_ax == 1: ax.set_title('SEDs as fit by SDSS+2MASS (excluding {})'.format(band_name))
            i_ax0 +=1
        elif mode == 1 and (i_ax1 < len(axs_1)):
            i_ax = i_ax1
            ax = axs_1[i_ax]
            if i_ax == 1: ax.set_title('SEDs as fit by USNOB1+2MASS')
            i_ax1 +=1
        if ax != None:
            ax.scatter( MODELS[0][1:][mask], obs[::2], c='k', marker='D', s=20, label='observations' )
            ax.scatter( MODELS[0][1:], model, c='b', marker='o', s=50, alpha=.5, label='model' )
            ax.scatter( MODELS[0][1:][i_band], true_band, c='r', marker='D', s=20, label='SDSS-{}'.format(band_name) )
            ax.invert_yaxis()
            if i_ax%pltsize == 0:
                ax.set_ylabel('Mag')
            if len(axs_0)-i_ax <= pltsize:
                ax.set_xlabel('Wavelength (A)')
            if i_ax == pltsize-1:
                ax.legend(loc=4)
    
    # now plot a histogram for each type
    plt.figure(3)
    alph = .5
    bns = map( lambda x: round(x,2), np.linspace(-2, 2, 50) )
    plt.hist( errors_0, bins=bns, alpha=alph, normed=True, color='g', label='SDSS+2MASS' )
    plt.hist( errors_1, bins=bns, alpha=alph, normed=True, color='b', label='USNOB1+2MASS' )
    plt.legend(loc='best')
    plt.ylabel('Normalized count')
    plt.xlabel('Error in {}-band (mag)'.format(band_name))
    plt.title('SDSS+2MASS: {} --- USNOB1+2MASS: {}'.format(len(errors_0), len(errors_1)) )
    plt.show()
    
    return errors_0, errors_1
        

# example: ra, dec = (314.136483, -6.081352)
def construct_SED( ra, dec, redden=True ):
    '''
    Construct the SED for a single object using SDSS
    magnitudes as well as synthetic photometry from a 
    modeled spectrum.
    
    REQUIRES OBJECT BE IN SDSS+2MASS.
    '''
    # simply assume the first object returned by each catalog
    #  is the relevant object
    mass, sdss, usnob = query_all(ra, dec, boxsize=1.)
    obs = np.hstack( (sdss[0][2:], mass[0][2:]) )
    
    if redden:
        # de-redden the observations
        filts = ['u','g','r','i','z','J','H','K']
        obs[::2] -= get_reddening( ra, dec, filts )
    
    # fit a model to the observations
    model, T = choose_model( obs, [1,1,1,1,1,0,0,0,1,1,1] )

    if redden:
        # redden the model and re-redden the obs
        model += get_reddening( ra, dec, ALL_FILTERS )
        obs[::2] += get_reddening( ra, dec, filts )
        
    mask = np.array([0,1,1,1,1,1,0,1,1,1,1,1]).astype(bool)
    full_obs = np.hstack( (sdss[0][2:], usnob[0][2:], mass[0][2:]) )
    plt.scatter( MODELS[0][mask], full_obs[::2], c='g', label='observations' )
    plt.scatter( MODELS[0][1:], model, c='b', marker='D', label='model' )
    plt.legend(loc='best')
    plt.xlabel('Wavelength (A)')
    plt.ylabel('Mag')
    plt.title( 'Complete SED, using model: {}K'.format(round(T)) )
    plt.show()


# example: ra, dec = (314.136483, -6.081352)
def test_sdss_interp( ra, dec, redden=True ):
    '''
    Construct the SED of a single object (present in all 3 catalogs)
    without using the SDSS mags, and see how well we estimate
    SDSS mags from USNOB+2MASS
    
    REQUIRES OBJECT TO BE IN ALL THREE CATALOGS.
    '''
    # query the catalogs
    mass, sdss, usnob = query_all(ra, dec, boxsize=1.)
    obs = np.hstack( (usnob[0][2:], mass[0][2:]) )
    
    if redden:
        # de-redden the observations
        filts = ['B','R','J','H','K']
        obs[::2] -= get_reddening( ra, dec, filts )
    
    # fit a model to the observations
    model, T = choose_model( obs, [0,0,0,0,0,0,1,1,1,1,1] )
    
    if redden:
        # redden the model and re-redden the observations
        model += get_reddening( ra, dec, ALL_FILTERS )
        filts = ['B','R','J','H','K']
        obs[::2] += get_reddening( ra, dec, filts )
        
    # plot it up
    mask = np.array([0,1,1,1,1,1,0,1,1,1,1,1]).astype(bool)
    full_obs = np.hstack( (sdss[0][2:], obs) )
    plt.scatter( MODELS[0][mask], full_obs[::2], c='g', label='observations' )
    plt.scatter( MODELS[0][1:], model, c='b', marker='D', label='model' )
    plt.legend(loc='best')
    plt.xlabel('Wavelength (A)')
    plt.ylabel('Mag')
    plt.title('SED as determined by fit to USNOB+2MASS, using model: {}K'.format(round(T)) )
    plt.show()


def show_SED( ra, dec, redden=True):
    '''
    show best-fit SED and observations for any object.
    '''
    # first, get all observations locally
    mass, sdss, usnob = query_all(ra, dec, boxsize=1.)
    obs = np.array([])
    mask = [0]
    if sdss != None:
        obs = np.hstack( (obs, sdss[0][2::2]) )
        mask += [1,1,1,1,1,0]
    else:
        mask += [0,0,0,0,0,0]
    if usnob != None:
        obs = np.hstack( (obs, usnob[0][2::2]) )
        mask += [1,1]
    else:
        mask += [0,0]
    if mass != None:
        obs = np.hstack( (obs, mass[0][2::2]) )
        mask += [1,1,1]
    else:
        mask += [0,0,0]
    mask = np.array(mask).astype(bool)
    
    # now, use the program above to get the modeled SEDs
    oc, fs, modes = produce_catalog( (ra, dec), 1., redden=redden, return_model=True )
    
    # plot up the object
    plt.scatter( MODELS[0][1:], fs[0], c='b', marker='o', s=50, alpha=.5, label='model' )
    plt.scatter( MODELS[0][mask], obs, c='r', marker='D', s=20, label='observations')
    plt.gca().invert_yaxis()
    plt.legend(loc=4)
    plt.xlabel('Wavelength (A)')
    plt.ylabel('Mag')
    if modes[0] == 0:
        mmm = 'SDSS+2MASS'
    else:
        mmm = 'USNOB1+2MASS'
    plt.title('SED determined with mode {}'.format(mmm) )
    plt.show()
    

