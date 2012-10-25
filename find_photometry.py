'''
OUTLINE:
a code that accepts:
  region, list of coordinates of stars
and returns:
  magnitude estimates for each item in the above list
 
 
TODO:
 - account for differences in BVR passband types?
 - account for difference between actual response of full system
   (SDSS throughputs) and simple response of filters (BVRYJHK)?
    - I think this is automatically accounted for in the photometry system
      definitions.
 - confirm y-band info; Y (stromgren vega mag), or Y abmag?

'''

import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from scipy.optimize import minimize
from sys import exit
import pyfits, ephem

try:
    models = np.load( open('all_models.npy','r') )
except:
    print 'cannot find models file'
    exit()


# open the dust extinction maps
try:
    hdu = pyfits.open('dust_maps/SFD_dust_1024_ngp.fits')[0]
    NSGP_n = hdu.header['LAM_NSGP']
    SCALE_n = hdu.header['LAM_SCAL']
    map_n = hdu.data
    hdu = pyfits.open('dust_maps/SFD_dust_1024_sgp.fits')[0]
    NSGP_s = hdu.header['LAM_NSGP']
    SCALE_s = hdu.header['LAM_SCAL']
    map_s = hdu.data
    map_dict = {}
    map_dict['N'] = [map_n, NSGP_n, SCALE_n]
    map_dict['S'] = [map_s, NSGP_s, SCALE_s]
except:
    print 'cannot find/open dust maps'
    exit()

############################################
# supporting functions
############################################

def parse_sdss8( s ):
    '''
    parse findsdss8 output
    
    s: a string as returned by findsdss8 using the flag "-e0,"
    returns: a list of objects, each line containing:
      [ra, dec, u, u_sigma, g, g_sigma, r, r_sigma, i, i_sigma, z, z_sigma]
    ''' 
    out = []
    lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
    for line in lines:
        lout = line.split(',')
        # RA and DEC
        if '+' in lout[1]:
            char = '+'
        else: 
            char = '-'
        ra  = float(lout[1].split(char)[0])
        dec = float(char + lout[1].split(char)[1])
        # magnitudes and errors
        #  in order: u, u_sig, g, g_sig, r, r_sig, i, i_sig, z, z_sig 
        mags = []
        for val in lout[4:9]:
            mags.append( float(val.split(':')[0]) )
            mags.append( float(val.split(':')[1]) )
        out.append( [ra, dec] + mags )
    return np.array(out)
    
    
def identify_matches( queried_stars, found_stars, match_radius=.5 ):
    '''
    match queried stars to found stars
    
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


def query_sdss8( center, width, star_coords, flags='' ):
    '''
    Query sdss8 server for info on stars in a field.
    Returns a list of magnitude arrays (when object is found) or None (when not found)
    
    center: decimal coordinates of center-of-field for query
    width: field size to request (in decimal degrees)
    star_coords: list of decimal coordinates of stars to match
    flags: additional flags to pass to findsdss8
    '''
    print 'asking sdss8...'
    # query about the field
    width = width/2.778e-4  #convert decimal degrees to arcseconds
    request = 'findsdss8 -c "{} {}" -bs {} -e0,'.format(center[0], center[1], width)
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    sdss_objects = parse_sdss8(o)
    if len(sdss_objects) == 0:
        # no matches at all in this field
        return [None]*len(star_coords)
    else:
        # identify matches
        matches = identify_matches(star_coords, sdss_objects[:,:2])
        out_list = []
        for i, match in enumerate(matches):
            if match != None:
                # return only the magnitudes
                out_list.append( sdss_objects[matches[i]][2::2] )
            else:
                # keep a None as a placeholder
                out_list.append(None)
        return out_list


def parse_nomad1( s ):
    '''
    parse nomad1 output
    
    s: a string as returned by nomad1
    returns: a numpy array of objects, each line containing:
      [ra, dec, B, V, R, J, H, K] (using Nones as placeholders)
    '''
    out = []
    lines = [lll for lll in s.split('\n') if lll and lll[0]!='#']
    for line in lines:
        lout = line.split('|')
        # RA and DEC
        if '+' in lout[2]:
            char = '+'
        else: 
            char = '-'
        ra  = float(lout[2].split(' ')[0].split(char)[0])
        dec = float(char + lout[2].split(' ')[0].split(char)[1])
        # magnitudes
        #  in order: [B, V, R, J, H, K] 
        mags = []
        # BVR + JHK
        for val in [v for v in lout[5].split(' ') if v] + [v for v in lout[6].split(' ') if v]:
            try:
                # strip off letter at end of magnitude
                #  FIXME: this ignores the differences between types of V bands, etcetera
                mags.append(float(val[:-1]))
            except:
                # land here if not detected or observed in this band
                mags.append(None)
        out.append( [ra, dec] + mags )
    return np.array(out)
    
    
def query_nomad1( center, width, star_coords, flags='' ):
    '''
    Query nomad1 server for info on stars in a field.
    Returns a list of magnitude arrays (when object is found) or None (when not found)
    width: field size to request (in decimal degrees)    
    center: decimal coordinates of center-of-field for query
    star_coords: list of decimal coordinates of stars to match
    flags: additional flags to pass to findnomad1
    '''
    print 'asking nomad...'
    width = width/2.778e-4  #convert decimal degrees to arcseconds
    request = 'findnomad1 -c "{} {}" -bs {}'.format(center[0], center[1], width)
    out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    nomad_objects = parse_nomad1(o)
    if len(nomad_objects) == 0:
        # no objects at all in this field
        return [None]*len(star_coords)
    else:
        # identify matches
        matches = identify_matches(star_coords, nomad_objects[:,:2])
        out_list = []
        for i, match in enumerate(matches):
            if match != None:
                # return the magnitudes
                out_list.append( nomad_objects[matches[i]][2:] )
            else:
                # keep a None as a placeholder
                out_list.append(None)
        return out_list

    
def find_field( star_coords, extend=.0015 ):
    '''
    determine the best field for a list of star coordinates,
    so we can perform only a single query.
    
    star_coords: a list of star coordinates, in decimal degrees
    extend: the buffer beyond the requested coordinates to add to the field
      also in decimal degrees
    
    returns: (coordinates of center), width_of_box
    '''
    ras, decs = zip(*star_coords)
    width_ra = (max(ras)-min(ras) + 2*extend)
    center_ra = np.mean(ras)
    width_dec = (max(decs)-min(decs) + 2*extend)
    center_dec = np.mean(decs)
    return (center_ra, center_dec), max(width_ra, width_dec)


def error_C(C, model, obs):
    '''
    C: number, a constant akin to distance modulus
    model: array-like, model absolute mags
    real: array-like, true observed mags

    returns: sum squared errors
    '''
    nm = model+C
    return np.sum( (nm-obs)**2 )


def choose_model( obs, mask, reddening=None, models=models ):
    '''
    find and return the best model for obs, where mask
    describes which of the modeled passbands are in obs.

    '''
    if reddening == None:
        reddening = np.zeros(len(models[0]))
    # rezero both the model and the obs, for comparision of SED shapes
    norm_obs = np.array(obs) - np.min( np.array(obs) )
    # go through all models and choose the one with the most similar shape
    sum_sqrs = []
    for i,model in enumerate(models):
        mod = model + reddening
        norm_mod = mod[mask] - np.min( mod[mask] )
        sum_sqrs.append(np.sum( (norm_mod-norm_obs)**2 ))
    best_model = models[ np.argmin(sum_sqrs) ] + reddening
    # fit for the additive constant (akin to distance modulus)
    res = minimize( error_C, 50., args=(best_model[mask], obs) )
    C = res.values()[5][0]
    # return all magnitudes for best model
    return best_model + C
    

def get_reddening( ra, dec, dust_map=map_dict ):
    '''
    return reddening values as set of corrections
    for modeled filters.
    
    assumes ra,dec are in decimal degrees.
    
    returns reddening for filters =
     [u, g, r, i, z, Y, B, V, R, J, H, K ]
    '''
    # coordinate-to-pixel mapping from the dust map fits header
    X_pix = lambda l,b,pole: np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.cos(l)*dust_map[pole][2]
    Y_pix = lambda l,b,pole: -dust_map[pole][1]*np.sqrt(1.-dust_map[pole][1]*np.sin(b))*np.sin(l)*dust_map[pole][2]
    # E(B-V)-to-mag mapping from http://www.astro.princeton.edu/~schlegel/dust/data/filter.txt
    filt_factor = np.array([5.155, 3.793, 2.751, 2.086, 1.479, 3.277, 4.315, 3.315, 2.673, 0.902, 0.576, 0.367])
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
    reddening = filt_factor*EBV
    return reddening
    
    
############################################
# main function
############################################

def full_query( star_coords ):
    '''
    Return the magnitude in all DECam bands for a list of objects.

    star_coords: coordinates of the stars for which you want magnitudes (must be a list)
    
    returns:
    an array with one line for each object in star_coords
     zeros if object not found in SDSS or NOMAD
     ([g,r,i,z,Y])  if object found in SDSS or NOMAD
    '''
    center, width = find_field( star_coords )
    # get estimate of reddening for this field
    reddening = get_reddening( *center )
    s_result = query_sdss8( center, width, star_coords )
    # return SDSS mags & Y for those that have SDSS
    final_result = np.zeros( (len(star_coords), 5) )
    nomad_requests = []
    for i,line in enumerate(s_result):
        if line == None:
             #object not in SDSS
            nomad_requests.append(i)
            continue
        mask = np.array([True]*5+[False]*7)
        model = choose_model( np.array(line), mask, reddening=reddening )
        # keep SDSS mags and stromgren Y
        final_row = np.append( line[1:], model[5] )
        final_result[i] = final_row
        
    # go though and fill in NOMAD results for objects not in SDSS
    if len(nomad_requests) > 0:
        # query NOMAD if we need to
        trim_star_coords = [star_coords[j] for j in nomad_requests]
        n_result = query_nomad1( center, width, trim_star_coords )
        # return estimated SDSS mags & Y
        for i,line in enumerate(n_result):
            i_final = nomad_requests[i] # relates final output array index to working index
            if line == None:
                # not in NOMAD either, so report nothing
                continue
            elif not any(line[:3]) or not any(line[3:]):
                # demand at least one optical and one IR observation to fit
                continue
            # otherwise, fit for a model &etc
            obs = np.array([val for val in line if val])
            mask = np.array( [False]*6 + [True if val else False for val in line] )
            model = choose_model( obs, mask, reddening=reddening )
            # return SDSS mags and stromgren Y
            final_row = model[1:6]
            final_result[i_final] = final_row
    
    return final_result



################################################
# tests
################################################


def test_me():
    '''
    should return a list with two succesful SDSS gets,
    two complete fails.
    '''
    query_locations = [ (314.111725, -6.051116), (314.113992, -6.049693),\
    (314.108888, -6.15000), (314.1504417, -06.0502722) ]
    result = full_query( query_locations )
    print '\nresult of test:\n\n'
    print result
    
def test_single_model_fit():
    ''' queries a single star, and shows how the interpolation goes '''
    query_locations = [(314.128304, -06.070742)]
    center, width = find_field( query_locations )
    s_result = query_sdss8( center, width, query_locations )[0]
    n_result = query_nomad1( center, width, query_locations )[0]
    reddening = get_reddening( *center )
    # find the best-fit model, fitting to NOMAD 
    model_mask = np.array( [False]*6 + [True if val else False for val in n_result])
    nomad_mask = np.array([True if val else False for val in n_result])
    best_model = choose_model( n_result[nomad_mask], model_mask, reddening=reddening )
    nored_model = choose_model( n_result[nomad_mask], model_mask)
    nomad_wl = [445., 551., 658., 1250., 1650., 2150.] #nm
    sdss_wl  = [354., 475., 622., 763., 905.] #nm
    y_wl = [1020.] #nm
    
    plt.scatter( sdss_wl+y_wl+nomad_wl, best_model, c='k', marker='x', label='model')
    plt.scatter( sdss_wl+y_wl+nomad_wl, nored_model, c='k', marker='+', label='model (no reddening)')
    plt.scatter( np.array(nomad_wl)[nomad_mask], n_result[nomad_mask] , c='g', label='nomad')
    plt.scatter( sdss_wl, s_result, c='r', label='sdss')
    
    plt.xlabel('wavelength (nm)')
    plt.ylabel('magnitude')
    plt.legend(loc=4)
    plt.gca().invert_yaxis()
    plt.show()


def test_model_fits( n=10 ):
    '''
    queries n stars, and shows how the interpolation goes
    NOTE: run get_standards.py first!
    '''
    nomad_wl = [445., 551., 658., 1250., 1650., 2150.] #nm
    sdss_wl  = [354., 475., 622., 763., 905.] #nm
    y_wl = [1020.] #nm
    
    query_locations = np.load(open('matches.npy','r'))
    center, width = find_field( query_locations )
    reddening = get_reddening( *center )
    
    s_results = query_sdss8( center, width, query_locations )
    n_results = query_nomad1( center, width, query_locations )
    
    for i,s_result in enumerate(s_results):
        if i>n: break
        n_result = n_results[i]
    
        # find the best-fit model, fitting to NOMAD 
        model_mask = np.array( [False]*6 + [True if val else False for val in n_result])
        nomad_mask = np.array([True if val else False for val in n_result])
        best_model = choose_model( n_result[nomad_mask], model_mask, reddening=reddening )
        nored_model = choose_model( n_result[nomad_mask], model_mask, reddening=None )

        plt.scatter( sdss_wl+y_wl+nomad_wl, best_model, c='k', marker='x', label='model')
        plt.scatter( sdss_wl+y_wl+nomad_wl, nored_model, c='k', marker='+', label='model (no reddening)')
        plt.scatter( np.array(nomad_wl)[nomad_mask], n_result[nomad_mask] , c='g', label='nomad')
        plt.scatter( sdss_wl, s_result, c='r', label='sdss')

        plt.title('object {}'.format(i))
        plt.xlabel('wavelength (nm)')
        plt.ylabel('magnitude')
        plt.legend(loc=4)
        plt.gca().invert_yaxis()
        plt.show()

