'''
A series of test functions writting to go along with 
 get_SEDs.py
 
USAGE:
$ python
> from get_SEDs import *
> [test function of choice]

'''
from get_SEDs import *

############################################
# TEST FUNCTIONS
############################################


def test_SDSS_errors( ra, dec, band_name='z', redden=False, size=900., plot=True ):
    '''
    Test the error accrued for all sources in a field when estimating 
     SDSS _-band photometry from all modes.  Errors in y-band photometry
     are expected to be similar to z-band errors.
    
    Run this on any field within the SDSS footprint.

    ra,dec: coordinates in decimal degrees
    band: SDSS passband to derive errors for
    redden: if True, account for galactic reddening when modeling photometry
    size: size of field to query around ra, dec (arseconds)
    plot: show/build plots?
    '''
    if band_name == 'u':
        sdss_mask = [0,0,1,1,1,1,1,1,1,1]  #these mask both mags and errors
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
        raise Exception('Incorrect band keyword!')
    sdss_mask = np.array(sdss_mask).astype(bool)

    mass, sdss, usnob = query_all(ra, dec, boxsize=size)

    object_mags = []
    band_mags = []
    sdss_mags = []
    #jmk0,jmk1,gmr0,gmr1 = [],[],[],[]  # used to plot model predictions as a function of colors
    errs0, errs1 = [],[]  # used to plot model predictions as a function of model fit
    modes = []
    object_coords = []

    # match sdss, usnob objects to 2mass objects
    if sdss != None:
        sdss_matches = identify_matches( mass[:,:2], sdss[:,:2] )
        usnob_matches = identify_matches( mass[:,:2], usnob[:,:2] )
    else:
        raise Exception('Must run this for coordinates in SDSS footprint!')

    # Go through 2MASS objects and assemble a catalog
    #  of all objects present in multiple catalogs
    for i,obj in enumerate(mass):
        if sdss_matches[i] != None and usnob_matches[i] !=None:
            i_sdss = sdss_matches[i]
            i_usnob = usnob_matches[i]
            
            # fit to SDSS+2MASS
            obs = np.hstack( (sdss[i_sdss][2:][sdss_mask], obj[2:]) )
            band = sdss[i_sdss][2:][~sdss_mask][0] # keep track of the true band mag
            object_mags.append( obs )
            sdss_mags.append( sdss[i_sdss][2:][sdss_mask])
            object_coords.append( obj[:2] )
            modes.append( 0 )
            band_mags.append( band )
            
            #gmr0.append( obj[2] - obj[4] )
            #jmk0.append( mass[i_mass][2] - mass[i_mass][-2] )

            # also fit to USNOB+2MASS
            obs = np.hstack( (usnob[i_usnob][2:], obj[2:]) )
            object_mags.append( obs )
            sdss_mags.append( sdss[i_sdss][2:][sdss_mask])
            object_coords.append( obj[:2] )
            modes.append( 1 )
            band = sdss[i_sdss][2:][~sdss_mask][0] # keep track of the true band mag
            band_mags.append( band )
            
            #gmr1.append( obj[2] - obj[4] )
            #jmk1.append( mass[i_mass][2] - mass[i_mass][-2] )
            
        elif sdss_matches[i] != None and usnob_matches[i] == None:

            i_sdss = sdss_matches[i]
            # fit to SDSS+2MASS
            obs = np.hstack( (sdss[i_sdss][2:][sdss_mask], obj[2:]) )
            band = sdss[i_sdss][2:][~sdss_mask][0]
            object_mags.append( obs )
            sdss_mags.append( sdss[i_sdss][2:][sdss_mask])
            object_coords.append( obj[:2] )
            modes.append( 0 )
            band_mags.append( band )
            
            #gmr0.append( obj[4] - obj[6] )
            #jmk0.append( mass[i_mass][2] - mass[i_mass][-2] )


    # now fit a model to each object, and construct the final SED,
    #  without including the relevant band.  Determine errors between
    #  predicted band and actual. Build an array of sample SEDs of
    #  the first 9 sources for each type of fit.
    # Modes are defined as:
    #  0 -> SDSS+2MASS; 1 -> USNOB+2MASS
    if plot:
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
            reddening = _get_reddening( ra,dec, ALL_FILTERS )
            model, T, err = choose_model_reddening( obs, mask, reddening )
        else:
            model, T, err = choose_model( obs, mask )
        if err > 2.: continue # impose a quality-of-fit cut
        
        # compare calculated z-mag to observed
        true_band = band_mags[i]
        guess_band = model[i_band]
        error = true_band - guess_band
        if mode == 0:
            errors_0.append(error)
            errs0.append(err)
        elif mode == 1:
            errors_1.append(error)
            errs1.append(err)
            
        if plot:
            # if this is a particularly egregrious one, plot it up
            if abs(error) > 1.:
                plt.figure()
                ax = plt.subplot(111)
                ax.scatter( MODELS[0][1:][mask], obs[::2], c='k', marker='D', s=20, label='observations' )
                ax.scatter( MODELS[0][1:], model, c='b', marker='o', s=50, alpha=.5, label='model' )
                ax.scatter( MODELS[0][1:][i_band], true_band, c='r', marker='D', s=20, label='SDSS-{}'.format(band_name) )
                ax.scatter( MODELS[0][1:6][sdss_mask[::2]], sdss_mags[i][::2], c='r', marker='x', s=20 )
                ax.invert_yaxis()
                ax.set_title( 'Terrible fit at {}'.format(object_coords[i]) )
            
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
    
    if plot:
        # now plot a histogram for each type
        plt.figure()
        alph = .5
        bns = map( lambda x: round(x,2), np.linspace(-2, 2, 50) )
        plt.hist( errors_0, bins=bns, alpha=alph, normed=True, color='g', label='SDSS+2MASS' )
        plt.hist( errors_1, bins=bns, alpha=alph, normed=True, color='b', label='USNOB1+2MASS' )
        plt.legend(loc='best')
        plt.ylabel('Normalized count')
        plt.xlabel('Error in {}-band (mag)'.format(band_name))
        plt.title('SDSS+2MASS: {} --- USNOB1+2MASS: {}'.format(len(errors_0), len(errors_1)) )

        '''  # these were used to determine whether I should use a color cut - i.e. do color=0 stars fit better?  answer: No
        plt.figure(4)
        plt.scatter( errors_0, jmk0, color='r', alpha=.8, label='J-K' )
        plt.scatter( errors_0, gmr0, color='b', alpha=.8, label='g-r' )
        plt.legend(loc='best')
        plt.ylabel('Color')
        plt.xlabel('Error in {}-band (mag)'.format(band_name))
        plt.title('SDSS+2MASS: {}'.format(len(errors_0)) )

        plt.figure(5)
        plt.scatter( errors_1, jmk1, color='r', alpha=.8, label='J-K' )
        plt.scatter( errors_1, gmr1, color='b', alpha=.8, label='g-r' )
        plt.legend(loc='best')
        plt.ylabel('Color')
        plt.xlabel('Error in {}-band (mag)'.format(band_name))
        plt.title('USNOB1+2MASS: {}'.format(len(errors_1)) )
        plt.show()
        '''
    
        # these are used to determine whether I should use a quality-of-fit cut. answer: yes
        plt.figure()
        plt.scatter( errors_0, errs0, color='r', alpha=.8, label='SDSS+2MASS' )
        plt.scatter( errors_1, errs1, color='b', alpha=.8, label='USNOB+2MASS' )
        plt.legend(loc='best')
        plt.ylabel('quality parameter')
        plt.xlabel('Error in {}-band (mag)'.format(band_name))
        plt.title('SDSS+2MASS: {} --- USNOB+2MASS: {}'.format(len(errors_0), len(errors_1)) )
        plt.show()

    return errors_0, errors_1
        

# example: ra, dec = (314.136483, -6.081352)
def construct_SED( ra, dec, redden=True ):
    '''
    Construct the SED for a single object using SDSS or USNOB
     and 2MASS magnitudes as well as synthetic photometry from a 
     modeled spectrum.
    
    REQUIRES OBJECT BE IN SDSS+2MASS or USNOB+2MASS.
    '''
    # simply assume the first object returned by each catalog
    #  is the relevant object
    mass, sdss, usnob = query_all(ra, dec, boxsize=1.)
    if sdss != None and mass != None:
        obs = np.hstack( (sdss[0][2:], mass[0][2:]) )
        mask = [1,1,1,1,1,0,0,0,1,1,1]
    elif usnob != None and mass != None:
        obs = np.hstack( (usnob[0][2:], mass[0][2:]) )
        mask = [0,0,0,0,0,0,1,1,1,1,1]
    else:
        raise Exception('Cannot find source!')
    
    if redden:
        reddening = _get_reddening( ra,dec, ALL_FILTERS )
        model, T, err = choose_model_reddening( obs, mask, reddening )
    else:
        model, T, err = choose_model( obs, mask )
        
    plt.scatter( MODELS[0][1:][mask], obs[::2], c='g', label='observations' )
    plt.scatter( MODELS[0][1:], model, c='b', marker='D', label='model' )
    plt.legend(loc='best')
    plt.xlabel('Wavelength (A)')
    plt.ylabel('Mag')
    plt.title( 'SED, using model: {}K --- error param: {}'.format(round(T), round(err,2)) )
    plt.show()


