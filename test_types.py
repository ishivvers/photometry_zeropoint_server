'''
A script that produces a confusion matrix between true and modeled
 classifications for stars in the Skiff catalog
'''
import matplotlib.pyplot as plt
import numpy as np
import get_SEDs as gs
import re
ALL_FILTERS = ['u','g','r','i','z','y','B','V','R','I','J','H','K']
TYPES = ['O','B','A','F','G','K','M']
# make a dictionary of the basic types for all of the models we use
d = np.load('all_models.npy')
t = np.loadtxt('pickles_uk.ascii',dtype=str)
types_dict = {}
for mod in map(int, d[1:,0]):
    print mod, t[mod]
    types_dict[mod] = re.findall('[A-Z]', t[mod][1])[0]


def parse_ra( inn ):
    '''
    Parse input RA string, either decimal degrees or sexagesimal HH:MM:SS.SS (or similar variants).
    Returns decimal degrees.
    '''
    # if simple float, assume decimal degrees
    try:
        ra = float(inn)
        return ra
    except:
        # try to parse with phmsdms:
        res = phmsdms(inn)
        ra = 15.*( res['vals'][0] + res['vals'][1]/60. + res['vals'][2]/3600. )
        return ra


def parse_dec( inn ):
    '''
    Parse input Dec string, either decimal degrees or sexagesimal DD:MM:SS.SS (or similar variants).
    Returns decimal degrees.
    '''
    # if simple float, assume decimal degrees
    try:
        dec = float(inn)
        return dec
    except:
        # try to parse with phmsdms:
        res = phmsdms(inn)
        dec = res['sign']*( res['vals'][0] + res['vals'][1]/60. + res['vals'][2]/3600. )
        return dec


def phmsdms(hmsdms):
    """
    +++ Pulled from python package 'angles' +++
    Parse a string containing a sexageismal number.
    
    This can handle several types of delimiters and will process
    reasonably valid strings. See examples.
    
    Parameters
    ----------
    hmsdms : str
        String containing a sexagesimal number.
    
    Returns
    -------
    d : dict
    
        parts : a 3 element list of floats
            The three parts of the sexagesimal number that were
            identified.
        vals : 3 element list of floats
            The numerical values of the three parts of the sexagesimal
            number.
        sign : int
            Sign of the sexagesimal number; 1 for positive and -1 for
            negative.
        units : {"degrees", "hours"}
            The units of the sexagesimal number. This is infered from
            the characters present in the string. If it a pure number
            then units is "degrees".
    """
    units = None
    sign = None
    # Floating point regex:
    # http://www.regular-expressions.info/floatingpoint.html
    #
    # pattern1: find a decimal number (int or float) and any
    # characters following it upto the next decimal number.  [^0-9\-+]*
    # => keep gathering elements until we get to a digit, a - or a
    # +. These three indicates the possible start of the next number.
    pattern1 = re.compile(r"([-+]?[0-9]*\.?[0-9]+[^0-9\-+]*)")
    # pattern2: find decimal number (int or float) in string.
    pattern2 = re.compile(r"([-+]?[0-9]*\.?[0-9]+)")
    hmsdms = hmsdms.lower()
    hdlist = pattern1.findall(hmsdms)
    parts = [None, None, None]
    
    def _fill_right_not_none():
        # Find the pos. where parts is not None. Next value must
        # be inserted to the right of this. If this is 2 then we have
        # already filled seconds part, raise exception. If this is 1
        # then fill 2. If this is 0 fill 1. If none of these then fill
        # 0.
        rp = reversed(parts)
        for i, j in enumerate(rp):
            if j is not None:
                break
        if  i == 0:
            # Seconds part already filled.
            raise ValueError("Invalid string.")
        elif i == 1:
            parts[2] = v
        elif i == 2:
            # Either parts[0] is None so fill it, or it is filled
            # and hence fill parts[1].
            if parts[0] is None:
                parts[0] = v
            else:
                parts[1] = v
                
    for valun in hdlist:
        try:
            # See if this is pure number.
            v = float(valun)
            # Sexagesimal part cannot be determined. So guess it by
            # seeing which all parts have already been identified.
            _fill_right_not_none()
        except ValueError:
            # Not a pure number. Infer sexagesimal part from the
            # suffix.
            if "hh" in valun or "h" in valun:
                m = pattern2.search(valun)
                parts[0] = float(valun[m.start():m.end()])
                units = "hours"
            if "dd" in valun or "d" in valun:
                m = pattern2.search(valun)
                parts[0] = float(valun[m.start():m.end()])
                units = "degrees"
            if "mm" in valun or "m" in valun:
                m = pattern2.search(valun)
                parts[1] = float(valun[m.start():m.end()])
            if "ss" in valun or "s" in valun:
                m = pattern2.search(valun)
                parts[2] = float(valun[m.start():m.end()])
            if "'" in valun:
                m = pattern2.search(valun)
                parts[1] = float(valun[m.start():m.end()])
            if '"' in valun:
                m = pattern2.search(valun)
                parts[2] = float(valun[m.start():m.end()])
            if ":" in valun:
                # Sexagesimal part cannot be determined. So guess it by
                # seeing which all parts have already been identified.
                v = valun.replace(":", "")
                v = float(v)
                _fill_right_not_none()
        if not units:
            units = "degrees"
            
    # Find sign. Only the first identified part can have a -ve sign.
    for i in parts:
        if i and i < 0.0:
            if sign is None:
                sign = -1
            else:
                raise ValueError("Only one number can be negative.")
                
    if sign is None:  # None of these are negative.
        sign = 1
        
    vals = [abs(i) if i is not None else 0.0 for i in parts]
    return dict(sign=sign, units=units, vals=vals, parts=parts)


# first calculate the type for all that we can
fskiff = '/o/ishivvers/zeropoint_code/data/mktypes.dat'
lines = open(fskiff,'r').readlines()
objs = []
for l in lines:
    try:
        ra = parse_ra( l[49:60] )
        dec = parse_dec( l[61:72] )
        true_type = l[83]
        if true_type not in TYPES:
            continue
        band = l[80]
        true_mag = float(l[75:80])
    except:
        print 'passing on line:\n',l
        continue
    try:
        c = gs.catalog( ra, dec, size=10 )
    except:
        print 'no match found at',ra,dec
        continue
    dists = np.sqrt( (ra-c.coords[:,0])**2 + (dec-c.coords[:,1])**2 )
    match = np.argmin(dists)
    mod_type = types_dict[c.models[match]]
    mod_mag = c.SEDs[match][ ALL_FILTERS.index(band) ]
    if np.abs( mod_mag-true_mag ) < 1.0:
        objs.append( [ra,dec,true_type,mod_type] )
    else:
        print 'discarding mismatch:',mod_mag,true_mag

# now make the plot
confusion = np.zeros( (len(TYPES), len(TYPES)) )
for obj in objs:
    # true type
    i = TYPES.index( obj[2] )
    # model type
    j = TYPES.index( obj[3] )
    confusion[i,j] += 1
np.save('confusion.np',confusion)

plt.imshow( confusion, interpolation='nearest', cmap='Greys' )
plt.colorbar()
plt.xticks( range(len(TYPES)), TYPES )
plt.yticks( range(len(TYPES)), TYPES )
plt.xlabel('True Type')
plt.ylabel('Model Type')
plt.savefig('confusion_matrix.png')
plt.show()
