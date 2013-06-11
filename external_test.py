"""
Suite of tests to compare my code's predictions to APASS, UKIDSS, SDSS, and USNOB.

To Do:
 - add error dict for all colors from APASS
 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile
import get_SEDs as gs

COORDS = [(335., .5), (350., .5), (1., .5),
          (15., .5), (30., .5)]
COLORS = ['b','g','r','orange','gray']
imdir = 'images/'

def compare_to_SDSS( field_center, field_width, clip=True, colors=COLORS ):
    '''
    Compare models built from USNOB to true SDSS mags.
    '''
    q = gs.online_catalog_query( field_center[0], field_center[1], field_width )
    sdss = q.query_sdss()
    if sdss == None:
        raise ValueError('No SDSS sources in this field!')
    c = gs.catalog( field_center, field_width, ignore_sdss=True )
    matches,tmp = gs.identify_matches( sdss[:,:2], c.coords )
    
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( sdss[i,2:][::2] )
            model_err.append( c.model_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    matched_modes.append( c.modes[match] )
    
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['u','g','r','i','z']):
        ccc = colors[i]
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip:
            err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.figure(1)
        plt.hist( err, bins=bins, alpha=.8, color=ccc, linewidth=3, histtype='step', normed=True, label=filt+": %.2f -- %.2f (%s)" %(med, mad, 'usnob'))
        
        plt.figure(2)
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err, np.abs(err), color=ccc, alpha=.5, marker='+' )
        
        dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
        delt_bin = dbins[1:]-dbins[:-1]
        inbin = np.digitize( model_err, dbins )
        onesigs = []
        for iii in range(1,len(dbins)):
            try:
                mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
            except:
                mmm = onesigs[-1]
            onesigs.append( mmm )
        x = list(dbins[1:]-delt_bin/2)
        # extend the onesigs out to the edges for interpolation
        x = [dbins[0]] + x + [dbins[-1]]
        onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
        plt.scatter( x, onesigs, color=ccc, label=filt, s=50, alpha=.9 )
        f_med = interp1d(x, onesigs)
        x2 = np.linspace(min(x), max(x), 1000)
        plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)    
    
    plt.figure(1)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('SDSS - model')
    plt.ylabel('Normalized count')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    plt.savefig(imdir+"sdss_errs_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure(2)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('Model Chi^2 (mags)')
    plt.ylabel('True error (mags)')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    plt.savefig("sdss_errs_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.show()


def web_plots_SDSS( size=1800., coords=COORDS, ignore='sdss' ):
    '''
    produce web plot and error report for SDSS predicted from USNOB
    '''
        
    medians = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    mads = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    onesig_errs = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    
    bins = np.linspace(-1, 1, 20)
    axs = []
    for ifield,field_center in enumerate(coords):
        print 'starting {}\n\n\n'.format(field_center)
        ccc = COLORS[ifield]
        
        q = gs.online_catalog_query( field_center[0], field_center[1], size )
        sdss = q.query_sdss()
        c = gs.catalog( field_center, size, ignore=ignore )
        matches,tmp = gs.identify_matches( sdss[:,:2], c.coords )
        
        model_matches = []
        obs_matches = []
        model_err = []
        models = []
        matched_modes = []
        matched_J = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_matches.append( c.SEDs[match] )
                obs_matches.append( sdss[i,2:][::2] )
                model_err.append( c.model_errors[match] )
                models.append( c.models[match] )
                matched_modes.append( c.modes[match] )
                matched_J.append( c.SEDs[match][gs.ALL_FILTERS.index('J')] )
        model_matches = np.array(model_matches)
        obs_matches = np.array(obs_matches)
        matched_modes = np.array(matched_modes)
        
        # temporary: a plot-up of some of the results
        '''
        from time import sleep
        numplot = 40
        plt.ion()
        plt.figure(20)
        plt.show()
        all_wl = [ gs.FILTER_PARAMS[key][0] for key in gs.ALL_FILTERS ]
        all_zp = np.array([ gs.FILTER_PARAMS[key][1] for key in gs.ALL_FILTERS ])
        sdss_wl = [3551., 4686., 6165., 7481., 8931.]
        sdss_zp = np.array( [8.6387e-9, 4.9607e-9, 2.8660e-9, 1.9464e-9, 1.3657e-9] )
        for i,model in enumerate(model_matches[:numplot]):
            print 'mode:',matched_modes[i]
            
            fluxmodel = np.log10( all_zp*10**(-.4*model) )
            fluxsdss  = np.log10( sdss_zp*10**(-.4*obs_matches[i]) )
            plt.clf()
            plt.scatter( all_wl, fluxmodel, marker='x', c='g', label='model' )
            plt.scatter( sdss_wl, fluxsdss, marker='o', c='r', label='sdss' )
            plt.title( '{} of {}'.format(i, numplot) )
            plt.legend(loc=4)
            if i==0: plt.gca().invert_yaxis()
            plt.draw()
            #raw_input('\nenter to continue')
            sleep(.25)
        plt.close()
        plt.ioff()
        '''
        for i,band in enumerate(['u','g','r','i','z']):
            # the histogram
            i_filt = gs.FILTER_PARAMS[band][-1]
            
            err = gs.clip_me(obs_matches[:,i] - model_matches[:,i_filt])
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[band].append(med); mads[band].append(mad)
            
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=2, histtype='step', normed=True, color=ccc,\
                           label="(%.2f, %.2f)" %(field_center[0], field_center[1]))
            
            # the scatterplot
            err = obs_matches[:,i] - model_matches[:,i_filt]

            #EDIT: a 3d scatterplot with x=chi2, y=Jmag, z=error
            
            plt.figure(i+5)
            if not ifield:
                ax = fig.add_subplot(111, projection='3d')
                axs.append(ax)
            else:
                ax = axs[ifield]
            ax.scatter3D( model_err, matched_J, err, alpha=.5, c=ccc )
            #plt.scatter( model_err, np.abs(err), alpha=.5, marker='+', c=ccc )
            dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
            delt_bin = dbins[1:]-dbins[:-1]
            inbin = np.digitize( model_err, dbins )
            onesigs = []
            for iii in range(1,len(dbins)):
                try:
                    mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
                except:
                    mmm = onesigs[-1]
                onesigs.append( mmm )
            x = list(dbins[1:]-delt_bin/2)
            # extend the medians out to the edges for interpolation
            x = [dbins[0]] + x + [dbins[-1]]
            onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
            #plt.scatter( x, onesigs, label="(%.2f, %.2f)" %(field_center[0], field_center[1]),\
            #            s=50, alpha=.9, c=ccc )
            f_med = interp1d(x, onesigs)
            x2 = np.linspace(min(x), max(x), 1000)
            #plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
            
            onesig_errs[band].append(onesigs)
            
            # a histogram of errors per model
            plt.figure(12)
            plt.scatter( models, err, alpha=.25, marker='.', c=ccc )

    
    
    for i,band in enumerate(['u','g','r','i','z']):
        plt.figure(i)
        plt.xlabel('Error in {}-band (mag)'.format(band))
        plt.ylabel('Normalized Count')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"sdss_errs_hist_%c.png" %(band), transparent=True)
        
        plt.figure(i+5)
        plt.xlabel('Model X^2, {}-band (mag)'.format(band))
        plt.ylabel('True Error')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"sdss_errs_scatp_%c.png" %(band), transparent=True)
        
    plt.figure(12)
    plt.xlabel( 'model #' )
    plt.ylabel( 'error' )
    plt.title( 'errors per stellar type (SDSS)' )
    
    return medians, mads, onesig_errs

    
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
            model_err.append( c.model_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['B','R']):
        ccc = filt.lower()
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip: err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.figure(1)
        plt.hist( err, bins=bins, linewidth=3, color=ccc, histtype='step', normed=True, label=filt+": %.2f -- %.2f" %(med, mad))
        
        plt.figure(2)
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err, np.abs(err), color=ccc, alpha=.5, marker='+')
        
        dbins = np.array( [0, .5, 1., 2., 5., 8.] )
        delt_bin = dbins[1:]-dbins[:-1]
        inbin = np.digitize( model_err, dbins )
        onesigs = []
        for iii in range(1,len(dbins)):
            try:
                mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
            except:
                mmm = onesigs[-1]
            onesigs.append( mmm )
        x = list(dbins[1:]-delt_bin/2)
        # extend the medians out to the edges for interpolation
        x = [dbins[0]] + x + [dbins[-1]]
        onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
        plt.scatter( x, onesigs, color=ccc, label=filt, s=50, alpha=.9 )
        f_med = interp1d(x, onesigs)
        x2 = np.linspace(min(x), max(x), 1000)
        plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
        
    plt.figure(1)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('USNOB1 - model')
    plt.ylabel('Normalized Count')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    plt.savefig(imdir+"usnob_errs_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure(2)
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.xlabel('Model Chi^2 (mags)')
    plt.ylabel('True error (mags)')
    plt.title('Errors over {} sources'.format(len(obs_matches)))
    plt.savefig(imdir+"usnob_errs_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.show()


def compare_to_UKIDSS( field_center, field_size, fn='data/ukidss_all.csv', clip=True, ignore_sdss=True ):
    '''
    # Note: UKIDSS is all in vegamag.
    #  yAB - yVega = 0.58693544 across all models.
    #  For now, just add that in to Y below.
    '''
    
    data = np.loadtxt( fn, delimiter=',')
    # get rid of any that don't have a yband
    dat = data[ data[:,2]>0 ]
    input_coords = dat[:,:2] # ra,dec
    Y = dat[:,2:4] + 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
    J = dat[:,4:6]           # mag, err
    H = dat[:,6:8]
    K = dat[:,8:10]
    
    c = gs.catalog( field_center, field_size, input_coords=input_coords, ignore_sdss=ignore_sdss )
    matches,tmp = gs.identify_matches( input_coords, c.coords, 3. )
    
    i_y = gs.FILTER_PARAMS['y'][-1]
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
            model_matches.append( c.SEDs[match,i_y] )
            obs_matches.append( Y[i][0] )
            model_err.append( c.full_errors[match,i_y] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    
    # the histogram
    plt.figure(1)
    bins = np.linspace(-1, 1, 20)
    err = obs_matches - model_matches
    if clip: err = gs.clip_me( err )
    med = np.median( err )
    mad = np.median( np.abs( err-np.median(err) ) )
    
    plt.hist( err, bins=bins, linewidth=3, histtype='step', normed=True )
    plt.xlabel('UKIDSS - model')
    plt.ylabel('Normalized Count')
    plt.title('Errors in y over {} sources'.format(len(obs_matches)))
    if ignore_sdss:
        plt.savefig(imdir+"ukidss_nosdss_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    else:
        plt.savefig(imdir+"ukidss_sdss_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    # the scatterplot
    plt.figure(2)
    err = obs_matches - model_matches
    plt.scatter( model_err, np.abs(err), alpha=.5, marker='+' )
    
    if ignore_sdss:
        dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
    else:
        dbins = np.array( [0, .5, 1., 2., 5., 8.] )
    delt_bin = dbins[1:]-dbins[:-1]
    inbin = np.digitize( model_err, dbins )
    onesigs = []
    for iii in range(1,len(dbins)):
        try:
            mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
        except:
            mmm = onesigs[-1]
        onesigs.append( mmm )
    x = list(dbins[1:]-delt_bin/2)
    # extend the medians out to the edges for interpolation
    x = [dbins[0]] + x + [dbins[-1]]
    onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
    plt.scatter( x, onesigs, s=50, alpha=.9 )
    f_med = interp1d(x, onesigs)
    x2 = np.linspace(min(x), max(x), 1000)
    plt.plot( x2, f_med(x2), '--', alpha=.9)
    
    plt.xlabel('Model Chi^2 (mags)')
    plt.ylabel('True error (mags)')
    plt.title('Errors in y over {} sources'.format(len(obs_matches)))
    if ignore_sdss:
        plt.savefig(imdir+"ukidss_nosdss_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    else:
        plt.savefig(imdir+"ukidss_sdss_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.show()


def web_plots_UKIDSS( size=1800., coords=COORDS, ignore='sdss' ):
    '''
    Produce web plots for UKIDSS y compared to model y
    '''
    
    data = np.loadtxt( 'data/ukidss_all.csv', delimiter=',')
    # get rid of any that don't have a yband
    dat = data[ data[:,2]>0 ]
    input_coords = dat[:,:2] # ra,dec
    Y = dat[:,2:4] + 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
    J = dat[:,4:6]           # mag, err
    H = dat[:,6:8]
    K = dat[:,8:10]
    
    medians, mads, onesig_errs = [],[],[]
    bins = np.linspace(-1, 1, 20)
    
    for i_field,field_center in enumerate(coords):
        print 'Starting field %.2f, %.2f\n\n' %(field_center[0], field_center[1])
        ccc = COLORS[i_field]
        
        c = gs.catalog( field_center, size, input_coords=input_coords, ignore=ignore )
        matches,tmp = gs.identify_matches( input_coords, c.coords )
        
        i_y = gs.FILTER_PARAMS['y'][-1]
        model_matches = []
        obs_matches = []
        model_err = []
        models = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_matches.append( c.SEDs[match,i_y] )
                obs_matches.append( Y[i][0] )
                model_err.append( c.model_errors[match] )
                models.append( c.models[match] )
        model_matches = np.array(model_matches)
        obs_matches = np.array(obs_matches)
        
        err = gs.clip_me( obs_matches - model_matches )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        medians.append(med)
        mads.append(mad)
        
        # the histogram
        plt.figure(1)
        plt.hist( err, bins=bins, alpha=.8, color=ccc, linewidth=3, histtype='step', normed=True,\
                       label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]))
                       
        # now the scatterplot
        err = obs_matches - model_matches
        
        plt.figure(2)
        plt.scatter( model_err, np.abs(err), c=ccc, alpha=.5, marker='+' )
        if 'sdss' in ignore:
            dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 8.] )
        delt_bin = dbins[1:]-dbins[:-1]
        inbin = np.digitize( model_err, dbins )
        onesigs = []
        for iii in range(1,len(dbins)):
            try:
                mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
            except:
                mmm = onesigs[-1]
            onesigs.append( mmm )
        x = list(dbins[1:]-delt_bin/2)
        # extend the medians out to the edges for interpolation
        x = [dbins[0]] + x + [dbins[-1]]
        onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
        plt.scatter( x, onesigs, c=ccc, label="(%.2f, %.2f)" %(field_center[0], field_center[1]), s=50, alpha=.9 )
        f_med = interp1d(x, onesigs)
        x2 = np.linspace(min(x), max(x), 1000)
        plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
        
        onesig_errs.append(onesigs)
        
        # error per model
        plt.figure(14)
        plt.scatter( models, err, alpha=.25, marker='.', c=ccc )
        
    plt.figure(1)
    plt.xlabel('Error in y-band (mag)')
    plt.ylabel('Normalized Count')
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"ukidss_errs_hist_y.png", transparent=True)
    
    plt.figure(2)
    plt.xlabel('Model Chi^2 (y-band)')
    plt.ylabel('True Error (mags)')
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"ukidss_errs_scatp_y.png", transparent=True)
    
    plt.figure(14)
    plt.xlabel( 'model #' )
    plt.ylabel( 'error' )
    plt.title( 'errors per stellar type (UKIDSS)' )
    
    return medians, mads, onesig_errs


def compare_to_APASS( field_center, field_size, fn='data/apass_all.csv', ignore_sdss=True, clip=True ):
    
    data = np.loadtxt( fn, delimiter=',', skiprows=1 )
    input_coords = data[:,:4:2]
    V = data[:,5:7]
    B = data[:,7:9]
    g = data[:,9:11]
    r = data[:,11:13]
    iband = data[:,13:15]
    obs = np.hstack( (g,r,iband,B,V) )
    
    c = gs.catalog( field_center, field_size, input_coords=input_coords, ignore_sdss=ignore_sdss )
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
    for i,filt in enumerate(['g','r','i','B','V']):
        ccc=COLORS[i]
        i_filt = gs.FILTER_PARAMS[filt][-1]
        
        # handle cases with a single bad band (value = -99.)
        obmatch = obs_matches[:,i][ obs_matches[:,i]>0]
        momatch = model_matches[:,i_filt][ obs_matches[:,i]>0]
        
        err = obmatch - momatch
        if clip: err = gs.clip_me(err)
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        # histogram
        plt.figure(1)
        plt.hist( err, label=filt+': %.2f -- %.2f' %(med, mad), color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
        
        # scatterplot
        plt.figure(2)
        err = obmatch - momatch
        merr = model_err[ obs_matches[:,i] > 0 ]
        plt.scatter( merr, np.abs(err), c=ccc, alpha=.5, marker='+')
        
        if ignore_sdss:
            dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 8.] )
        delt_bin = dbins[1:]-dbins[:-1]
        inbin = np.digitize( merr, dbins )
        medians = []
        for iii in range(1,len(dbins)):
            mmm = np.median(np.abs(err[inbin==iii]))
            if np.isnan(mmm):
                mmm = medians[-1]
            medians.append( mmm )
        x = list(dbins[1:]-delt_bin/2)
        # extend the medians out to the edges for interpolation
        x = [dbins[0]] + x + [dbins[-1]]
        medians = [medians[0]] + medians + [medians[-1]]
        plt.scatter( x, medians, color=ccc, s=50, alpha=.9,\
                    label=filt )
        f_med = interp1d(x, medians)
        x2 = np.linspace(min(x), max(x), 1000)
        plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
    
    plt.figure(1)
    plt.xlabel('APASS - model')
    plt.ylabel('Normalized Count')
    plt.title('Errors over {} sources'.format(len(err)))
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"apass_errs_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure(2)
    plt.xlabel('Model Chi^2 (mags)')
    plt.ylabel('True error (mags)')
    plt.title('Errors over {} sources'.format(len(err)))
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"apass_errs_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.show()


def web_plots_APASS( size=1800, ignore_sdss=True, coords=COORDS ):
    '''
    Produce plots and errors from comparisions with APASS data.
    '''
    data = np.loadtxt( 'data/apass_all.csv', delimiter=',', skiprows=1 )
    input_coords = data[:,:4:2]
    V = data[:,5:7]
    B = data[:,7:9]
    g = data[:,9:11]
    r = data[:,11:13]
    iband = data[:,13:15]
    obs = np.hstack( (g,r,iband,B,V) )
    
    medians = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    mads = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    onesig_errs = { 'g':[], 'r':[], 'i':[], 'B':[], 'V':[] }
    bins = np.linspace(-1, 1, 20)
    
    for i_field,field_center in enumerate(coords):
        
        c = gs.catalog( field_center, size, input_coords=input_coords, ignore_sdss=ignore_sdss )
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
        
        for i,filt in enumerate(['g','r','i','B','V']):
            ccc=COLORS[i_field]
            i_filt = gs.FILTER_PARAMS[filt][-1]
            
            # handle cases with a single bad band (value = -99.)
            obmatch = obs_matches[:,i][ obs_matches[:,i]>0]
            momatch = model_matches[:,i_filt][ obs_matches[:,i]>0]
            
            err = gs.clip_me( obmatch - momatch )
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[filt].append(med); mads[filt].append(mad)
            
            # histogram
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True,\
                           label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]))
            
            # scatterplot
            plt.figure(i+5)
            err = obmatch - momatch
            merr = model_err[ obs_matches[:,i] > 0 ]
            plt.scatter( merr, np.abs(err), c=ccc, alpha=.5, marker='+')
            
            if ignore_sdss:
                dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
            else:
                dbins = np.array( [0, .5, 1., 2., 5., 8.] )
            delt_bin = dbins[1:]-dbins[:-1]
            inbin = np.digitize( merr, dbins )
            onesigs = []
            for iii in range(1,len(dbins)):
                try:
                    mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
                except:
                    mmm = onesigs[-1]
                onesigs.append( mmm )
            x = list(dbins[1:]-delt_bin/2)
            # extend the medians out to the edges for interpolation
            x = [dbins[0]] + x + [dbins[-1]]
            onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
            plt.scatter( x, onesigs, color=ccc, s=50, alpha=.9,\
                        label="(%.2f, %.2f)" %(coords[i_field][0], coords[i_field][1]) )
            f_med = interp1d(x, onesigs)
            x2 = np.linspace(min(x), max(x), 1000)
            plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
            
            onesig_errs[filt].append(onesigs)
            
                        
    for i,filt in enumerate(['g','r','i','B','V']):
        plt.figure(i)
        plt.xlabel('Error in {}-band (mag)'.format(filt))
        plt.ylabel('Normalized Count')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"apass_errs_hist_%c.png" %(filt), transparent=True)
        
        plt.figure(i+5)
        plt.xlabel('Model X^2, {}-band (mag)'.format(filt))
        plt.ylabel('True Error')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"apass_errs_scatp_%c.png" %(filt), transparent=True)
        
    return medians, mads, onesig_errs

    
def check_APASS_phot(size=1800 ):
    '''
    To check how close APASS and SDSS photometry are.
    
    OBSOLETE.
    Result: they are quite close!
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


def compare_to_Stetson( field, ignore='sdss', clip=True, colors=['b','g','r','orange'] ):
    input_coords = np.loadtxt('data/'+field+'.pos.clean', skiprows=1)
    data = np.loadtxt('data/'+field+'.pho.clean', skiprows=1)
    obs = data[:,1:] #B,V,R,I
    
    field_center,field_size = gs.find_field( input_coords )
    c = gs.catalog( field_center, max(field_size), input_coords=input_coords, ignore=ignore )
    matches,tmp = gs.identify_matches( input_coords, c.coords )
    
    model_matches = []
    obs_matches = []
    model_err = []
    for i,match in enumerate(matches):
        if match >= 0:
            model_matches.append( c.SEDs[match] )
            obs_matches.append( obs[i] )
            model_err.append( c.model_errors[match] )
    model_matches = np.array(model_matches)
    obs_matches = np.array(obs_matches)
    model_err = np.array(model_err)
    
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['B','V','R','I']):
        ccc=colors[i]
        i_filt = gs.FILTER_PARAMS[filt][-1]
        
        # handle cases with a single bad band (value = -99.)
        obmatch = obs_matches[:,i][ obs_matches[:,i]!=99.999 ]
        momatch = model_matches[:,i_filt][ obs_matches[:,i]!=99.999 ]
        
        err = obmatch - momatch
        if clip: err = gs.clip_me(err)
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        # histogram
        plt.figure(1)
        plt.hist( err, label=filt+': %.2f -- %.2f' %(med, mad), color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
        
        # scatterplot
        plt.figure(2)
        err = obmatch - momatch
        merr = model_err[ obs_matches[:,i]!=99.999 ]
        plt.scatter( merr, np.abs(err), c=ccc, alpha=.5, marker='+')
        
        if ignore_sdss:
            dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 8.] )
        delt_bin = dbins[1:]-dbins[:-1]
        inbin = np.digitize( merr, dbins )
        onesigs = []
        for iii in range(1,len(dbins)):
            try:
                mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
            except:
                mmm = onesigs[-1]
            onesigs.append( mmm )
        x = list(dbins[1:]-delt_bin/2)
        # extend the medians out to the edges for interpolation
        x = [dbins[0]] + x + [dbins[-1]]
        onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
        plt.scatter( x, onesigs, color=ccc, s=50, alpha=.9,\
                    label=filt )
        f_med = interp1d(x, onesigs)
        x2 = np.linspace(min(x), max(x), 1000)
        plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
        
    
    plt.figure(1)
    plt.xlabel('Stetson - model')
    plt.ylabel('Normalized Count')
    plt.title('Errors over {} sources'.format(len(err)))
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"stetson_errs_hist_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.figure(2)
    plt.xlabel('Model Chi^2 (mags)')
    plt.ylabel('True error (mags)')
    plt.title('Errors over {} sources'.format(len(err)))
    leg = plt.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.0)
    plt.savefig(imdir+"stetson_errs_scatp_%.1f_%.1f.png" %(field_center[0], field_center[1]), transparent=True)
    
    plt.show()


FIELDS = ['L112', 'L113', 'L92', 'L95', 'PG0231'] #'IC1613', 
def web_plots_Stetson( fields=FIELDS, ignore='sdss', colors=['b','g','r','orange','grey','yellow'] ):
    
    medians = { 'B':[], 'V':[], 'R':[], 'I':[] }
    mads = { 'B':[], 'V':[], 'R':[], 'I':[] }
    onesig_errs = { 'B':[], 'V':[], 'R':[], 'I':[] }
    
    for ifield,field in enumerate(fields):
        print 'starting {}\n\n\n'.format(field)
        ccc = colors[ifield]
        
        input_coords = np.loadtxt('data/'+field+'.pos.clean', skiprows=1)
        data = np.loadtxt('data/'+field+'.pho.clean', skiprows=1)
        obs = data[:,1:] #B,V,R,I
    
        field_center,field_size = gs.find_field( input_coords )
        c = gs.catalog( field_center, max(field_size), input_coords=input_coords, ignore=ignore )
        matches,tmp = gs.identify_matches( input_coords, c.coords )
        
        model_matches = []
        obs_matches = []
        model_err = []
        models = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_matches.append( c.SEDs[match] )
                obs_matches.append( obs[i] )
                model_err.append( c.model_errors[match] )
                models.append( c.models[match] )
        model_matches = np.array(model_matches)
        obs_matches = np.array(obs_matches)
        model_err = np.array(model_err)
        
        bins = np.linspace(-1, 1, 20)
        for i,filt in enumerate(['B','V','R','I']):
            i_filt = gs.FILTER_PARAMS[filt][-1]
            
            # handle cases with a single bad band (value = -99.)
            obmatch = obs_matches[:,i][ obs_matches[:,i]!=99.999 ]
            momatch = model_matches[:,i_filt][ obs_matches[:,i]!=99.999 ]
            goodmodels = np.array(models)[obs_matches[:,i]!=99.999]
            
            err = gs.clip_me(obmatch - momatch)
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[filt].append(med); mads[filt].append(mad)
            
            # histogram
            plt.figure(i)
            plt.hist( err, label="(%.2f, %.2f)" %(field_center[0], field_center[1]),\
                      color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
                      
            # scatterplot
            plt.figure(i+4)
            err = obmatch - momatch
            merr = model_err[ obs_matches[:,i]!=99.999 ]
            plt.scatter( merr, np.abs(err), c=ccc, alpha=.5, marker='+')
            
            if 'sdss' in ignore:
                dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
            else:
                dbins = np.array( [0, .5, 1., 2., 5., 8.] )
            delt_bin = dbins[1:]-dbins[:-1]
            inbin = np.digitize( merr, dbins )
            onesigs = []
            for iii in range(1,len(dbins)):
                try:
                    mmm = scoreatpercentile( np.abs(err[inbin==iii]), 68.2 )
                except:
                    mmm = onesigs[-1]
                onesigs.append( mmm )
            x = list(dbins[1:]-delt_bin/2)
            # extend the medians out to the edges for interpolation
            x = [dbins[0]] + x + [dbins[-1]]
            onesigs = [onesigs[0]] + onesigs + [onesigs[-1]]
            plt.scatter( x, onesigs, color=ccc, s=50, alpha=.9,\
                        label="(%.2f, %.2f)" %(field_center[0], field_center[1]) )
            f_med = interp1d(x, onesigs)
            x2 = np.linspace(min(x), max(x), 1000)
            plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
            
            onesig_errs[filt].append(onesigs)
            
            # error per model
            plt.figure(13)
            plt.scatter( goodmodels, err, alpha=.25, marker='.', c=ccc )
            
            
    for i,band in enumerate(['B','V','R','I']):
        plt.figure(i)
        plt.xlabel('Error in {}-band (mag)'.format(band))
        plt.ylabel('Normalized Count')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"stetson_errs_hist_%c.png" %(band), transparent=True)
        
        plt.figure(i+4)
        plt.xlabel('Model X^2, {}-band (mag)'.format(band))
        plt.ylabel('True Error')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.0)
        plt.savefig(imdir+"stetson_errs_scatp_%c.png" %(band), transparent=True)
    
    plt.figure(13)
    plt.xlabel( 'model #' )
    plt.ylabel( 'error' )
    plt.title( 'errors per stellar type (Stetson)' )
    
    return medians, mads, onesig_errs



##########################################
# comparing to empirical power laws
##########################################

def PL_sdss2stet( fields=FIELDS, colors=['b','g','r','orange','grey','yellow'] ):
    '''
    Produce error diagrams akin to those above, but by running power law relations,
     not modeling the sources.
    Attempt to predict BVRI from SDSS ugriz.
    '''
    medians = { 'B':[], 'V':[], 'R':[], 'I':[] }
    mads = { 'B':[], 'V':[], 'R':[], 'I':[] }
    
    for ifield,field in enumerate(fields):
        ccc = colors[ifield]
        
        input_coords = np.loadtxt('data/'+field+'.pos.clean', skiprows=1)
        data = np.loadtxt('data/'+field+'.pho.clean', skiprows=1)
        obs = data[:,1:] #B,V,R,I
        
        field_center,field_size = gs.find_field( input_coords )
        q = gs.online_catalog_query( field_center[0], field_center[1], max(field_size) )
        matches,tmp = gs.identify_matches( input_coords, q.SDSS[:,:2] )
        
        matched_coords = []
        matched_sdss = []
        matched_stet = []
        for i,match in enumerate(matches):
            if match >= 0:
                matched_coords.append( input_coords[i] )
                matched_sdss.append( q.SDSS[match,2::2] ) #throw away sigmas
                matched_stet.append( obs[i] )
        matched_coords = np.array(matched_coords)
        matched_sdss = np.array(matched_sdss)
        matched_stet = np.array(matched_stet)
        
        # compare all predictions
        Bpred1 = matched_sdss[:,0] - 0.8116*(matched_sdss[:,0] - matched_sdss[:,1]) + 0.1313
        Bpred2 = matched_sdss[:,1] + 0.3130*(matched_sdss[:,1] - matched_sdss[:,2]) + 0.2271
        
        Vpred1 = matched_sdss[:,1] - 0.2906*(matched_sdss[:,0] - matched_sdss[:,1]) + 0.0885
        Vpred2 = matched_sdss[:,1] - 0.5784*(matched_sdss[:,1] - matched_sdss[:,2]) - 0.0038
        
        Rpred1 = matched_sdss[:,2] - 0.1837*(matched_sdss[:,1] - matched_sdss[:,2]) - 0.0971
        Rpred2 = matched_sdss[:,2] - 0.2936*(matched_sdss[:,2] - matched_sdss[:,3]) - 0.1439
        
        Ipred1 = matched_sdss[:,2] - 1.2444*(matched_sdss[:,2] - matched_sdss[:,3]) - 0.3820
        Ipred2 = matched_sdss[:,3] - 0.3780*(matched_sdss[:,3] - matched_sdss[:,4]) - 0.3974
        
        preds = [ [Bpred1, Bpred2], [Vpred1,Vpred2], [Rpred1,Rpred2], [Ipred1,Ipred2] ]
        '''
        Lupton equations from here:
        http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html
        
        B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
        B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107
        
        V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
        V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054
        
        R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
        R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
        
        I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
        I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063
        ''' 
            
        bins = np.linspace(-1, 1, 20)
        for i,band in enumerate(['B','V','R','I']):
            # handle cases with a single bad band (value = -99.)
            mask = matched_stet[:,i] != 99.999
            err1 = gs.clip_me( matched_stet[:,i][mask] - preds[i][0][mask] )
            err2 = gs.clip_me( matched_stet[:,i][mask] - preds[i][1][mask] )
            
            medians[band].append( [np.median(err1), np.median(err2)] )
            mads[band].append( [np.median( np.abs( err1-np.median(err1) ) ), np.median( np.abs( err2-np.median(err2) ) )] )
                        
            plt.figure(i)
            plt.hist( err1, label="{} - 1".format(field),\
                      color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
            plt.hist( err2, label="{} - 2".format(field),\
                      color=ccc, bins=bins, ls='dashed', linewidth=3, histtype='step', normed=True)
                      
    for i,band in enumerate(['B','V','R','I']):
        plt.figure(i)
        plt.title( band )
        plt.xlabel( 'Stetson - Power Law' )
        plt.ylabel( 'Normalized Count' )
        plt.legend()
    plt.show()
    return medians, mads


def PL_mass2izy( size=3600., coords=COORDS ):
    '''
    Produce error diagrams akin to those above, but by running power law relations,
     not modeling the sources.
    Attempt to predict zy from 2mass.
    '''
    medians = { 'i':[], 'z':[], 'y':[] }
    mads = { 'i':[], 'z':[], 'y':[] }
    
    data = np.loadtxt( 'data/ukidss_all.csv', delimiter=',')
    # get rid of any that don't have a yband
    dat = data[ data[:,2]>0 ]
    input_coords = dat[:,:2] # ra,dec
    y = dat[:,2] #+ 0.58693544 #+ 0.634   # correction to ABmag from hewett 2006
    
    bins = np.linspace(-1, 1, 20)
    for i_field,field_center in enumerate(coords):
        print 'Starting field %.2f, %.2f\n\n' %(field_center[0], field_center[1])
        ccc = COLORS[i_field]
    
        q = gs.online_catalog_query( field_center[0], field_center[1], size )
        matches_mass,tmp = gs.identify_matches( input_coords, q.MASS[:,:2] )
        matches_sdss,tmp = gs.identify_matches( input_coords, q.SDSS[:,:2] )
        
        y_match = []
        z_match = []
        i_match = []
        mass_match = []
        for i,match in enumerate(matches_mass):
            if match >= 0 and matches_sdss[i] >= 0:
                y_match.append( y[i] )
                z_match.append( q.SDSS[ matches_sdss[i] ][-2] )
                i_match.append( q.SDSS[ matches_sdss[i] ][-4] )
                mass_match.append( q.MASS[match][2::2] )
        y_match = np.array(y_match)
        z_match = np.array(z_match)
        i_match = np.array(i_match)
        mass_match = np.array(mass_match)
        
        """
        Equations from WFCam team:
        http://apm49.ast.cam.ac.uk/surveys-projects/wfcam/technical/photometry
        Z_wfcam  =  J +  0.95*(J - H)
        Y_wfcam  =  J +  0.50*(J - H)  + 0.08
        zAB = Z_wfcam + 0.533
        yAB = Y_wfcam + 0.634 # not necessary, since I'm using y_wfcam above (from UKIDSS)
        
        Equations from Vista Cam team:
        http://apm49.ast.cam.ac.uk/surveys-projects/vista/technical/the-vista-guide-camera-predicting-magnitudes-from-2mass-photometry
        i = J + (J-K)*1.175 + 0.459 (for H-K <= 0.15)
        i = J + (J-K)*1.175 + 0.459 + (H-K-0.15)*3.529   (for H-K >  0.15)
        """
        
        z_pred = mass_match[:,0] + .95*(mass_match[:,0] - mass_match[:,1]) + 0.533
        y_pred = mass_match[:,0] + .50*(mass_match[:,0] - mass_match[:,1]) + 0.08
        # apply two rules for i translation:
        HKup = mass_match[:,1] - mass_match[:,2] > .15
        i_pred = np.zeros_like( y_pred )
        for i,up in enumerate(HKup):
            J,H,K = mass_match[i]
            if up:
                i_pred[i] = J + (J-K)*1.175 + 0.459 + (H-K-0.15)*3.529
            else:
                i_pred[i] = J + (J-K)*1.175 + 0.459
        erry = gs.clip_me( y_match - y_pred )
        errz = gs.clip_me( z_match - z_pred )
        erri = gs.clip_me( i_match - i_pred )
        
        medians['i'].append( np.median(erri) )
        medians['z'].append( np.median(errz) )
        medians['y'].append( np.median(erry) )
        mads['i'].append( np.median( np.abs( erri-np.median(erri) ) ) )
        mads['z'].append( np.median( np.abs( errz-np.median(errz) ) ) )
        mads['y'].append( np.median( np.abs( erry-np.median(erry) ) ) )
        
        
        plt.figure(1)
        plt.hist( errz, label="{} - {}".format(field_center[0], field_center[1]),\
                  color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
        plt.figure(2)
        plt.hist( erry, label="{} - {}".format(field_center[0], field_center[1]),\
                  color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
        plt.figure(3)
        plt.hist( erri, label="{} - {}".format(field_center[0], field_center[1]),\
                color=ccc, bins=bins, linewidth=3, histtype='step', normed=True)
    plt.figure(1)
    plt.xlabel( 'sdss z - predicted z')
    plt.ylabel( 'normalized count' )
    plt.legend( loc='best' )
    plt.figure(2)
    plt.xlabel( 'ukidss y - predicted y')
    plt.ylabel( 'normalized count' )
    plt.legend( loc='best' )
    plt.figure(3)
    plt.xlabel( 'sdss i - predicted i')
    plt.ylabel( 'normalized count' )
    plt.legend( loc='best' )
    plt.show()
    return medians, mads
        
    


if __name__ == '__main__':
    # construct a dictionary of the mean errors,
    #  used to construct errors for catalogs.
    # saves to file a dictionary of form:
    # err_dict[mode][band] = y_array
    # err_dict[mode]['x'] = x_array
    #  then err = f(chi^2) = interp1d(x_array, y_array)
    
    import pickle
    
    def clear_all_figs(max=100):
        for i in range(max):
            if plt.fignum_exists(i):
                plt.figure(i)
                plt.clf()
        return
    
    err_dict = {}
    
    # first, mode 2 (USNOB)
    err_dict[2] = {}
    dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
    delt_bin = dbins[1:]-dbins[:-1]
    x = list(dbins[1:]-delt_bin/2)
    x = [dbins[0]] + x + [dbins[-1]]
    err_dict[2]['x'] = x
    
    med,mad,onesig = web_plots_SDSS(size=3600., ignore=['sdss','apass'])
    clear_all_figs()
    for band in onesig.keys():
        errs = np.mean( np.array(onesig[band]), axis=0 )
        err_dict[2][band] = errs
        
    med,mad,onesig = web_plots_UKIDSS(size=3600., ignore=['sdss','apass'])
    clear_all_figs()
    errs = np.mean( np.array(onesig), axis=0 )
    err_dict[2]['y'] = errs
    
    med,mad,onesig = web_plots_Stetson( ignore=['sdss','apass'] )
    clear_all_figs()
    for band in onesig.keys():
        if band in ['B','R']: continue
        errs = np.mean( np.array(onesig[band]), axis=0 )
        err_dict[2][band] = errs
    
    # now, mode 0 (SDSS)
    err_dict[0] = {}
    dbins = np.array( [0, .5, 1., 2., 5., 8.] )
    delt_bin = dbins[1:]-dbins[:-1]
    x = list(dbins[1:]-delt_bin/2)
    x = [dbins[0]] + x + [dbins[-1]]
    err_dict[0]['x'] = x
    
    med,mad,onesig = web_plots_UKIDSS(size=3600., ignore=['apass','usnob'])
    clear_all_figs()
    errs = np.mean( np.array(onesig), axis=0 )
    err_dict[0]['y'] = errs
    
    med,mad,onesig = web_plots_Stetson(ignore=['apass','usnob'])
    clear_all_figs()
    for band in onesig.keys():
        errs = np.mean( np.array(onesig[band]), axis=0 )
        err_dict[0][band] = errs
    
    # finally, mode 1 (apass)
    err_dict[1] = {}
    dbins = np.array( [0, .25, .5, 1., 1.5, 2.5] )
    delt_bin = dbins[1:]-dbins[:-1]
    x = list(dbins[1:]-delt_bin/2)
    x = [dbins[0]] + x + [dbins[-1]]
    err_dict[1]['x'] = x
    
    med,mad,onesig = web_plots_SDSS(size=3600., ignore=['sdss','usnob'])
    clear_all_figs()
    for band in onesig.keys():
        if band in ['g','r','i']: continue
        errs = np.mean( np.array(onesig[band]), axis=0 )
        err_dict[1][band] = errs
        
    med,mad,onesig = web_plots_UKIDSS(size=3600., ignore=['sdss','usnob'])
    clear_all_figs()
    errs = np.mean( np.array(onesig), axis=0 )
    err_dict[1]['y'] = errs
    
    med,mad,onesig = web_plots_Stetson( ignore=['sdss','usnob'] )
    clear_all_figs()
    for band in onesig.keys():
        if band in ['B','V']: continue
        errs = np.mean( np.array(onesig[band]), axis=0 )
        err_dict[1][band] = errs
        
    # now save our error dictionary to file
    pickle.dump( err_dict, open('err_dict.p','w') )
    
        
        
'''
ERROR ANALYSIS RESULTS

March 22, 2013


# One-sigma errors on the median value for a 1-degree field #

SDSS from USNOB
 u - 0.34
 g - 0.30
 r - 0.18
 i - 0.12
 z - 0.07
SDSS from APASS
 u - 0.13
 z - 0.05

UKIDSS from USNOB
 y - 0.03
UKIDSS from SDSS
 y - 0.03
UKIDSS from APASS
 y - 0.04  #real! probably due to fewer number of sources

Stetson from USNOB
 B - 0.48
 V - 0.24
 R - 0.14
 I - 0.10
Stetson from SDSS
 B - 0.07
 V - 0.07
 R - 0.05
 I - 0.03
Stetson from APASS
 R - 0.03
 I - 0.02

'''