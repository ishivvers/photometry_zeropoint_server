'''
OUTLINE:
a code that accepts:
  region, list of coordinates of stars
and returns:
  sdss magnitude estimates for each item in the above list
 
 
TODO:
 - account for differences in BVR passband types in parse_nomad1
 - assume (and fit to) real SEDs instead of interpolating when
    determining SDSS magnitudes
'''

import numpy as np
from subprocess import Popen, PIPE
from scipy.interpolate import interp1d



############################################
# supporting functions
############################################

def parse_sdss8(s):
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
        diffs = [(abs(star[0]-other[0])**2 + abs(star[1]-other[1])**2)**.5 for other in found_stars]
        if min(diffs) < match_radius:
            matches.append(np.argmin(diffs))
        else:
            matches.append(None) 
    return matches


def query_sdss8( center, star_coords, fieldsize=275, flags='' ):
    '''
    Query sdss8 server for info on stars in a field.
    Returns a list of magnitude arrays (when object is found) or None (when not found)
    
    center: decimal coordinates of center-of-field for query
    star_coords: list of decimal coordinates of stars to match
    fieldsize: size of field to request in arcseconds (default for DECam: 275")
    flags: additional flags to pass to findsdss8
    '''
    print 'asking sdss8...'
    # query about the field
    request = 'findsdss8 -c "{} {}" -bs {} -e0,'.format(center[0], center[1], fieldsize)
    out = Popen(request.format(center[0], center[1], fieldsize), shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    sdss_objects = parse_sdss8(o)
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
    
    
def query_nomad1( center, star_coords, fieldsize=275, flags='' ):
    '''
    Query nomad1 server for info on stars in a field.
    Returns a list of magnitude arrays (when object is found) or None (when not found)
    
    center: decimal coordinates of center-of-field for query
    star_coords: list of decimal coordinates of stars to match
    fieldsize: size of field to request in arcseconds (default for DECam: 275")
    flags: additional flags to pass to findnomad1
    '''
    print 'asking nomad...'
    request = 'findnomad1 -c "{} {}" -bs {}'.format(center[0], center[1], fieldsize)
    out = Popen(request.format(center[0], center[1], fieldsize), shell=True, stdout=PIPE, stderr=PIPE)
    o,e = out.communicate()
    # parse the response
    nomad_objects = parse_nomad1(o)
    # identify matches
    matches = identify_matches(star_coords, nomad_objects[:,:2])
    out_list = []
    for i, match in enumerate(matches):
        if match != None:
            # return the magnitues
            out_list.append( nomad_objects[matches[i]][2:] )
        else:
            # keep a None as a placeholder
            out_list.append(None)
    return out_list


def Nomad2Jansky( mag, band ):
    '''
    convert NOMAD magnitude to Janskies
    conversions from
     http://astro.wku.edu/strolger/UNITS.txt
     http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
    '''
    if band == 'B':
        Fnu0 = 4260.  #Jy
    elif band == 'V':
        Fnu0 = 3540.
    elif band == 'R':
        Fnu0 = 3080.
    elif band == 'J':
        Fnu0 = 1594.
    elif band == 'H':
        Fnu0 = 1024.
    elif band == 'K':
        Fnu0 = 666.7
    return Fnu0 * 10**(-.4*mag)


def Jansky2mag( f ):
    '''
    convert flux density (in Janskies, at central wavelength of SDSS band) to magnitude.
    I assume SDSS filters are true AB mags, and get conversions from
     http://www.sdss.org/dr5/algorithms/fluxcal.html
    '''
    f0 = 3631.
    return -2.5 * np.log10(f/f0)


def Nomad2SDSS( nomad_mags ):
    ''' interpolate what SDSS magnitudes you can from a set of NOMAD magnitudes '''
    # [B, V, R, J, H, K], ignoring differences between types of BVR
    bands = ['B', 'V', 'R', 'J', 'H', 'K']
    nomad_wl = [445., 551., 658., 1250., 1650., 2150.] #nm
    # [u, g, r, i, z]
    sdss_wl  = [354., 475., 622., 763., 905.] #nm
    # convert to janskies
    x,y = [],[]
    for i,mag in enumerate(nomad_mags):
        if mag == None:
            continue
        else:
            x.append( nomad_wl[i] )
            y.append( Nomad2Jansky(mag, bands[i]) )
    # interpolate in Janskies
    if len(x) < 3:
        kind = 'linear'
    else:
        kind = 'cubic'
    f = interp1d(x, y, kind=kind)
    # turn Janskies into SDSS mags
    sdss_mags = []
    for i, wl in enumerate(sdss_wl):
        if not x[0] < wl < x[-1]:
            # asking for a wl outside our interpolated region
            sdss_mags.append(None)
        else:
            flux = f(wl)
            sdss_mags.append( round(Jansky2mag(flux), 2) )
    return np.array(sdss_mags)
    

############################################
# main function
############################################

def full_query( center, star_coords ):
    '''
    return the magnitude in sdss bands for objects in field.

    center: center of the chip (in decimal coordinates)
    star_coords: coordinates of the stars for which you want magnitudes (must be a list)
    
    returns:
    a list containing one item for each object in star_coords:
     None - if object not found in either SDSS or NOMAD
     [u,g,r,i,z] - if object found in SDSS
     interpolated [u,g,r,i,z] - if object found in NOMAD but not SDSS
    '''
    result = query_sdss8( center, star_coords )
    # go though and fill in NOMAD results
    nomad_requests = []
    for i,line in enumerate(result):
        if line == None:
            nomad_requests.append( star_coords[i] )
    if len(nomad_requests) > 0:
        # query NOMAD if we need to
        n_result = query_nomad1( center, nomad_requests )
        # replace unknown SDSS values with values interpolated from NOMAD
        for i,line in enumerate(result):
            if line == None:
                nomad = n_result.pop(0)
                if nomad == None:
                    # found nothing, report nothing
                    result[i] = None
                    continue
                # otherwise, convert nomad-reported magnitudes to SDSS
                result[i] = Nomad2SDSS( nomad )
    return result





################################################
# tests
################################################


def test_me():
    ''' should return a list with two succesful SDSS gets, one fail, and one interpolated from NOMAD '''
    field_center    = (314.118, -6.0526)
    query_locations = [(314.111725, -6.051116), (314.113992, -6.049693),\
    (314.108888, -6.15000), (314.1504417, -06.0502722) ]
    result = full_query( field_center, query_locations)
    print '\nresult of test:\n\n'
    print result
    
def test_interpolation():
    field_center    = (314.118, -6.0526)
    query_locations = [(314.1252969, -6.0715939)]
    s_result = query_sdss8( field_center, query_locations )[0]
    n_result = query_nomad1( field_center, query_locations )[0]
    s_interp = Nomad2SDSS( n_result )
    print '\ntest result:\n\n'
    print 'SDSS:        ', s_result
    print 'interpolated:', s_interp
    print 'difference:  ', [round(s_result[i]-s_interp[i],2) if s_interp[i]!=None else None for i in range(len(s_result))]
    
