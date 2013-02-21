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
    
    matches,tmp = gs.identify_matches( sdss[:,:2], c.coords )
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( sdss[i,2:][::2] )
            model_err.append( c.full_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    model_err = np.array(model_err)
    
    colors = ['b','g','r','orange','gray']
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['u','g','r','i','z']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip: err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, label=filt+": %.2f -- %.2f" %(med, mad))
        
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('SDSS - model')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    plt.savefig("sdss_errs_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    
    plt.figure()
    for i,filt in enumerate(['u','g','r','i','z']):
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err[:,i_filt], np.abs(err), color=colors[i], label=filt )
    plt.legend(loc='best')
    plt.xlabel('Reported error (mags)')
    plt.ylabel('True error (mags)')
    
    plt.show()


def web_plots_SDSS( size=1800. ):
    '''
    produce web plot and error report for SDSS predicted from USNOB
    '''
    '''
    coords= [(210., 50.), (210., 20.), (185., 30.),
              (150., 50.), (150., 20.)]
    '''
    
    coords= [(335., .5), (350., .5), (1., .5),
              (15., .5), (30., .5)]
              
    medians = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    mads = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    bins = np.linspace(-1, 1, 20)
    for field_center in coords:
        print 'starting {}\n\n\n'.format(field_center)
        q = gs.online_catalog_query( field_center[0], field_center[1], size )
        sdss = q.query_sdss()
        c = gs.catalog( field_center, size, ignore_sdss=True )

        matches,tmp = gs.identify_matches( sdss[:,:2], c.coords )
        model_matches = []
        obs_matches = []
        for i,match in enumerate(matches):
          if match >= 0:
              model_matches.append( c.SEDs[match] )
              obs_matches.append( sdss[i,2:][::2] )
        model_matches = np.array(model_matches)
        obs_matches = np.array(obs_matches)
        
        for i,band in enumerate(['u','g','r','i','z']):
            i_filt = gs.FILTER_PARAMS[band][-1]
            err = gs.clip_me(obs_matches[:,i] - model_matches[:,i_filt])
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
        
            medians[band].append(med); mads[band].append(mad)
            
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, label="(%.2f, %.2f)" %(field_center[0], field_center[1]))
    
    for i,band in enumerate(['u','g','r','i','z']):
        plt.figure(i)
        plt.xlabel('Error in {}-band (mag)'.format(band))
        plt.ylabel('Normalized Count')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig("sdss_errs_%c.png" %(band), transparent=True)
    return medians, mads
    

    

def compare_to_USNOB( field_center, field_width, clip=True ):
    '''
    Compare models built from SDSS to true USNOB mags.
    '''
    q = gs.online_catalog_query( field_center[0], field_center[1], field_width )
    usnob = q.query_usnob1()
    c = gs.catalog( field_center, field_width )
    
    matches,tmp = gs.identify_matches( usnob[:,:2], c.coords )
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
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
    
    data = np.loadtxt( fn, delimiter=',')
    # go through and keep only the objects within our field
    fs = 2.778e-4*field_size
    new_data = []
    for row in data:
        if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[1] < field_center[1]+fs):
            new_data.append(row)
    data = np.array(new_data)
    # get rid of any that don't have a yband
    dat = data[ data[:,2]>0 ]
    input_coords = dat[:,:2] # ra,dec
    Y = dat[:,2:4] + 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
    J = dat[:,4:6]           # mag, err
    H = dat[:,6:8]
    K = dat[:,8:10]
    
    field_center, field_width = gs.find_field( input_coords )
    c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=True )
    matches,tmp = gs.identify_matches( input_coords, c.coords, 3. )
    
    i_y = 5  #index of yband in full_errors and in SEDs
    model_match = []
    obs_match = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
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
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('UKIDSS - model')
    plt.title('Errors over {} sources'.format(len(obs_match)))
    plt.savefig("ukidss_errs_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure()
    plt.scatter( model_err, np.abs(error), color='g' )
    plt.xlabel('Reported error (mags)')
    plt.ylabel('True error (mags)')
    
    plt.show()


def web_plots_UKIDSS( size=1800., ignore_sdss=True ):
    '''
    Produce web plots for UKIDSS y compared to model y
    '''
    coords= [(335., .5), (350., .5), (1., .5),
              (15., .5), (30., .5)]
    all_data = np.loadtxt( 'data/ukidss_all.csv', delimiter=',')
    medians, mads = [],[]
    bins = np.linspace(-1, 1, 20)
    for i_field,field_center in enumerate(coords):
        # go through and keep only the objects within our field
        fs = 2.778e-4*size
        data = []
        for row in all_data:
            if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[1] < field_center[1]+fs):
                data.append(row)
        data = np.array(data)
        # get rid of any that don't have a yband
        dat = data[ data[:,2]>0 ]
        input_coords = dat[:,:2] # ra,dec
        Y = dat[:,2:4] + 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
        J = dat[:,4:6]           # mag, err
        H = dat[:,6:8]
        K = dat[:,8:10]
        
        field_center, field_width = gs.find_field( input_coords )
        c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=ignore_sdss )
        matches,tmp = gs.identify_matches( input_coords, c.coords )
        
        i_y = gs.FILTER_PARAMS['y'][-1]
        model_match = []
        obs_match = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_match.append( c.SEDs[match,i_y] )
                obs_match.append( Y[i][0] )
        model_match = np.array(model_match)
        obs_match = np.array(obs_match)
        
        err = gs.clip_me( obs_match - model_match )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        medians.append(med)
        mads.append(mad)
        
        plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]))
    plt.xlabel('Error in y-band (mag)')
    plt.ylabel('Normalized Count')
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig("ukidss_errs_y.png", transparent=True)
    return medians, mads


'''
Error results over 5 fields: coords= [(335., .5), (350., .5), (1., .5),
          (15., .5), (30., .5)]

-- y no sdss --
median: .023
mad: .110

-- y sdss --
median: .021
mad: .071
'''

def compare_to_APASS( field_center, field_size, fn='data/apass_all.csv', ignore_sdss=True, clip=True ):
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
    c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=ignore_sdss )
    matches,tmp = gs.identify_matches( input_coords, c.coords )
    
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
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
    
    plt.figure(1)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig("apass_errs_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure(2)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.show()


def web_plots_APASS( size=1800, ignore_sdss=True ):
    '''
    Produce plots and errors from coparisions with APASS data.
    '''
    coords= [(335., .5), (350., .5), (1., .5),
              (15., .5), (30., .5)]
    all_data = np.loadtxt( 'data/apass_all.csv', delimiter=',')
    medians = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    mads = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    bins = np.linspace(-1, 1, 20)
    fs = 2.778e-4*size
    
    for i_field,field_center in enumerate(coords):
        data = []
        for row in all_data:
            if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[2] < field_center[1]+fs):
                data.append(row)
        data = np.array(data)
        input_coords = data[:,:4:2]
        V = data[:,5:7]
        B = data[:,7:9]
        g = data[:,9:11]
        r = data[:,11:13]
        iband = data[:,13:15]
        obs = np.hstack( (g,r,iband,B,V) )
        
        field_center, field_width = gs.find_field( input_coords )
        c = gs.catalog( field_center, max(field_width), input_coords=input_coords, ignore_sdss=ignore_sdss )
        matches,tmp = gs.identify_matches( input_coords, c.coords )
        
        model_matches = []
        obs_matches = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_matches.append( c.SEDs[match] )
                obs_matches.append( obs[i][::2] )
        model_matches = np.array(model_matches)
        obs_matches = np.array(obs_matches)
        
        for i,color in enumerate(['g','r','i','B','V']):
            i_color = gs.FILTER_PARAMS[color][-1]

            error = obs_matches[:, i] - model_matches[:,i_color]
            err = gs.clip_me(error)
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[color].append(med); mads[color].append(mad)
            
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]))

    for i,color in enumerate(['g','r','i','B','V']):
        plt.figure(i)
        plt.xlabel('Error in {}-band (mag)'.format(color))
        plt.ylabel('Normalized Count')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig("apass_errs_{}.png".format(color), transparent=True)
    return medians, mads


    
def check_APASS_phot(size=1800 ):
    '''
    To check how close APASS and SDSS photometry are.
    '''
    coords= [(335., .5), (350., .5), (1., .5),
              (15., .5), (30., .5)]
    all_data = np.loadtxt( 'data/apass_all.csv', delimiter=',')
    medians = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    mads = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    bins = np.linspace(-1, 1, 20)

    fs = 2.778e-4*size
    for i_field,field_center in enumerate(coords):
        data = []
        for row in all_data:
            if (field_center[0]-fs < row[0] < field_center[0]+fs) and (field_center[1]-fs < row[2] < field_center[1]+fs):
                data.append(row)
        data = np.array(data)
        input_coords = data[:,:4:2]
        V = data[:,5:7]
        B = data[:,7:9]
        g = data[:,9:11]
        r = data[:,11:13]
        iband = data[:,13:15]
        obs = np.hstack( (g,r,iband,B,V) )
        
        field_center, field_width = gs.find_field( input_coords )
        q = gs.online_catalog_query( field_center[0], field_center[1], max(field_width) )
        sdss = q.query_sdss()
        usnob = q.query_usnob1()
        sdss_matches,tmp = gs.identify_matches( input_coords, sdss[:,:2] )
        usnob_matches,tmp = gs.identify_matches( input_coords, usnob[:,:2] )
        
        sdss_matched = []
        s_apass_matched = []
        for i,match in enumerate(sdss_matches):
            if match >= 0:
                sdss_matched.append( sdss[match][2::2] )
                s_apass_matched.append( obs[i][::2] )
        sdss_matched = np.array(sdss_matched)
        s_apass_matched = np.array(s_apass_matched)
        
        usnob_matched = []
        u_apass_matched = []
        for i,match in enumerate(usnob_matches):
            if match >= 0:
                usnob_matched.append( usnob[match][2::2] )
                u_apass_matched.append( obs[i][::2] )
        usnob_matched = np.array(usnob_matched)
        u_apass_matched = np.array(u_apass_matched)
        
        
        for i,color in enumerate(['g','r','i','B']):
            i_color = gs.FILTER_PARAMS[color][-1]
            if color=='B':
                i_color = 0
                error = u_apass_matched[:, i] - usnob_matched[:,i_color]
            else:
                error = s_apass_matched[:, i] - sdss_matched[:,i_color]
            err = gs.clip_me(error)
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[color].append(med)
            mads[color].append(mad)
            
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]))
            plt.xlabel('Error in {}-band (mag)'.format(color))
            plt.ylabel('Normalized Count')
            leg = plt.legend(loc='best', fancybox=True)
            leg.get_frame().set_alpha(0.0)
            plt.savefig("apass_errs_{}.png".format(color), transparent=True)
            
    return medians, mads

