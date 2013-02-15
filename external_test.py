"""
Suite of tests to compare my code's predictions to APASS, UKIDSS, SDSS, and USNOB.

Must run from my_code folder.
"""

import numpy as np
import matplotlib.pyplot as plt

import get_SEDs as gs


def compare_to_SDSS( field_center, field_width, clip=True ):
    '''
    Compare models built from USNOB to true SDSS mags.
    '''
    q = gs.online_catalog_query( field_center[0], field_center[1], field_width )
    sdss = q.query_sdss()
    c = gs.catalog( field_center, field_width, ignore_sdss=True )
    
    matches = gs.identify_matches( sdss[:,:2], c.coords )
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match != None:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( sdss[i,2:][::2] )
            model_err.append( c.full_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    model_err = np.array(model_err)
    
    colors = ['purple','b','g','orange','r']
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['u','g','r','i','z']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip: err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.hist( err, bins=bins, linewidth=3, histtype='step', normed=True, label=filt+": %.2f -- %.2f" %(med, mad))
    plt.legend( loc='best' )
    plt.xlabel('SDSS - model')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    
    plt.figure()
    for i,filt in enumerate(['u','g','r','i','z']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err[:,i_filt], np.abs(err), color=colors[i], label=filt )
    plt.legend(loc='best')
    plt.xlabel('Reported error (mags)')
    plt.ylabel('True error (mags)')
    
    plt.show()


def compare_to_USNOB( field_center, field_width, clip=True ):
    '''
    Compare models built from SDSS to true USNOB mags.
    '''
    q = gs.online_catalog_query( field_center[0], field_center[1], field_width )
    usnob = q.query_usnob1()
    c = gs.catalog( field_center, field_width )
    
    matches = gs.identify_matches( usnob[:,:2], c.coords )
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match != None:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( usnob[i,2:][::2] )
            model_err.append( c.full_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    model_err = np.array(model_err)
    
    colors = ['b','r']
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['B','R']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip: err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.hist( err, bins=bins, linewidth=3, histtype='step', normed=True, label=filt+": %.2f -- %.2f" %(med, mad))
    plt.legend( loc='best' )
    plt.xlabel('USNOB1 - model')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    
    plt.figure()
    for i,filt in enumerate(['B','R']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err[:,i_filt], np.abs(err), color=colors[i], label=filt )
    plt.legend(loc='best')
    plt.xlabel('Reported error (mags)')
    plt.ylabel('True error (mags)')
    
    plt.show()
    
    
    
def compare_to_UKIDSS( field_center, field_size, fn='data/ukidss_in_SDSS.csv', clip=True ):
    '''
    # need to use model file with yband in vegamags
    #  MUST MAKE EXCHANGE IN get_SEDs.py: use 'all_models_yVega_P.npy'
    #  Note: yAB - yVega = 0.58693544 across all modelS
    #   For now, just add that in to Y below.
    '''
    
    data = np.loadtxt( fn, delimiter=',')[:ncomp]
    # go through and keep only the objects within our field
    fs = 2.778e-4*field_size
    new_data = []
    for row in data:
        if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[1] < field_center[1]+fs):
            new_dat.append(row)
    data = np.array(new_data)
    # get rid of any that don't have a yband
    dat = data[ data[:,2]>0 ]
    input_coords = dat[:,:2] # ra,dec
    Y = dat[:,2:4] + 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
    J = dat[:,4:6]           # mag, err
    H = dat[:,6:8]
    K = dat[:,8:10]
    
    field_center, field_width = gs.find_field( input_coords )
    c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=False )
    matches = gs.identify_matches( input_coords, c.coords, 3. )
    
    i_y = 5  #index of yband in full_errors and in SEDs
    model_match = []
    obs_match = []
    model_err = []
    for i,match in enumerate(matches):
        if match != None:
            model_match.append( c.SEDs[match,i_y] )
            obs_match.append( Y[i][0] )
            model_err.append( c.full_errors[match,i_y] )
    model_match = np.array(model_match)
    obs_match = np.array(obs_match)
    
    bins = np.linspace(-1, 1, 20)
    error = obs_match - model_match
    if clip:
        err = gs.clip_me( error )
    else:
        err = error
    med = np.median( err )
    mad = np.median( np.abs( err-np.median(err) ) )
    
    plt.figure()    
    plt.hist( err, bins=bins, linewidth=3, histtype='step', normed=True, label='y'+": %.2f -- %.2f" %(med, mad))
    plt.legend( loc='best' )
    plt.xlabel('UKIDSS - model')
    plt.title('Errors over {} sources'.format(len(obs_match)))
    
    plt.figure()
    plt.scatter( model_err, np.abs(error), color='g' )
    plt.xlabel('Reported error (mags)')
    plt.ylabel('True error (mags)')
    
    plt.show()


def compare_to_APASS( field_center, field_size, fn='data/apass_in_SDSS.csv', clip=True ):
    dat = np.loadtxt( fn, delimiter=',', skiprows=1 )
    # go through and keep only the objects within our field
    fs = 2.778e-4*field_size
    new_dat = []
    for row in dat:
        if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[2] < field_center[1]+fs):
            new_dat.append(row)
    dat = np.array(new_dat)
    input_coords = dat[:,:4:2]
    V = dat[:,5:7]
    B = dat[:,7:9]
    g = dat[:,9:11]
    r = dat[:,11:13]
    iband = dat[:,13:15]
    obs = np.hstack( (g,r,iband,B,V) )

    field_center, field_width = gs.find_field( input_coords )

    c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=False )
    matches = gs.identify_matches( input_coords, c.coords )
    
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match != None:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( obs[i][::2] )
            model_err.append( c.model_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    model_err = np.array(model_err)
    
    bins = np.linspace(-1, 1, 20)
    for i,color in enumerate(['g','r','i','B','V']):
        if color=='g': c='g'
        elif color=='r': c='r'
        elif color=='i': c='orange'
        elif color=='B': c='b'
        elif color=='V': c='grey'
        i_color = gs.FILTER_PARAMS[color][-1]
        
        error = obs_matches[:, i] - model_matches[:,i_color]
        if clip:
            err = gs.clip_me(error)
        else:
            err = error
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.figure(1)
        plt.hist( err, label=color+': %.2f -- %.2f' %(med, mad), color=c, bins=bins, linewidth=3, histtype='step', normed=True)
        plt.title('Errors over {} sources'.format(len(err)))
        plt.xlabel('APASS - model'.format(color,color))
        
        plt.figure(2)
        plt.scatter( model_err, np.abs(error), color=c, label=color )
        plt.xlabel('Reported error (mags)')
        plt.ylabel('True error (mags)')
    
    plt.figure(1); plt.legend(loc='best')
    plt.figure(2); plt.legend(loc='best')
    plt.show()


def check_APASS_phot(ncomp=500, fn='data/apass_in_SDSS.csv', clip=True):
    '''
    To check how close APASS and SDSS photometry are.
    '''
    dat = np.loadtxt( fn, delimiter=',', skiprows=1 )[:ncomp]
    input_coords = dat[:,:4:2]
    V = dat[:,5:7]
    B = dat[:,7:9]
    g = dat[:,9:11]
    r = dat[:,11:13]
    iband = dat[:,13:15]

    field_center, field_width = gs.find_field( input_coords )
    q = gs.online_catalog_query( field_center[0], field_center[1], max(field_width) )
    sdss = q.query_sdss()

    matches = gs.identify_matches( input_coords, sdss[:,:2] )

    sdss_match = []
    g_match, r_match, i_match = [],[],[]
    for i,match in enumerate(matches):
        if match == None: continue
        if (g[i,0] < 0) or (r[i,0] < 0) or (iband[i,0] < 0): continue # NA values
        sdss_match.append( sdss[match] )
        g_match.append( g[i] )
        r_match.append( r[i] )
        i_match.append( iband[i] )
    sdss = np.array(sdss_match)
    g = np.array(g_match)
    r = np.array(r_match)
    iband = np.array(i_match)

    plt.figure()
    bins = np.linspace(-1, 1, 20)
    for i in [4,6,8]:
        if i==4:
            l = c = 'g'
            band = g
        elif i==6:
            l = c = 'r'
            band = r
        elif i==8:
            l = 'i'; c = 'orange'
            band = iband
        err = sdss[:,i] - band[:,0]
        if clip:
            err = gs.clip_me( err )
        else:
            err = error
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        plt.hist( err, label=l+': %.2f -- %.2f' %(med, mad), color=c, bins=bins, linewidth=3, histtype='step', normed=True)
    plt.legend( loc='best' )
    plt.xlabel( 'sdss - model' )
    plt.title('comparing APASS and SDSS')
    plt.show()

