"""
Suite of tests to compare my code's predictions to APASS, UKIDSS, SDSS, and USNOB.
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
    
    bins = np.linspace(-1, 1, 20)
    for i,filt in enumerate(['u','g','r','i','z']):
        ccc = colors[i]
        i_filt = gs.FILTER_PARAMS[filt][-1]
        err = obs_matches[:,i] - model_matches[:,i_filt]
        if clip: err = gs.clip_me( err )
        med = np.median( err )
        mad = np.median( np.abs( err-np.median(err) ) )
        
        plt.figure(1)
        plt.hist( err, bins=bins, alpha=.8, color=ccc, linewidth=3, histtype='step', normed=True, label=filt+": %.2f -- %.2f" %(med, mad))
        
        plt.figure(2)
        err = obs_matches[:,i] - model_matches[:,i_filt]
        plt.scatter( model_err, np.abs(err), color=ccc, alpha=.5, marker='+' )
        
        dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
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


def web_plots_SDSS( size=1800., coords=COORDS ):
    '''
    produce web plot and error report for SDSS predicted from USNOB
    '''
        
    medians = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    mads = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    onesig_errs = { 'u':[], 'g':[], 'r':[], 'i':[], 'z':[] }
    
    bins = np.linspace(-1, 1, 20)
    for ifield,field_center in enumerate(coords):
        print 'starting {}\n\n\n'.format(field_center)
        ccc = COLORS[ifield]
        
        q = gs.online_catalog_query( field_center[0], field_center[1], size )
        sdss = q.query_sdss()
        c = gs.catalog( field_center, size, ignore_sdss=True )
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
        
        for i,band in enumerate(['u','g','r','i','z']):
            # the histogram
            i_filt = gs.FILTER_PARAMS[band][-1]
            err = gs.clip_me(obs_matches[:,i] - model_matches[:,i_filt])
            med = np.median( err )
            mad = np.median( np.abs( err-np.median(err) ) )
            medians[band].append(med); mads[band].append(mad)
            
            plt.figure(i)
            plt.hist( err, bins=bins, alpha=.8, linewidth=3, histtype='step', normed=True, color=ccc,\
                           label="(%.2f, %.2f)" %(field_center[0], field_center[1]))
            
            # the scatterplot
            err = obs_matches[:,i] - model_matches[:,i_filt]
            
            plt.figure(i+5)
            plt.scatter( model_err, np.abs(err), alpha=.5, marker='+', c=ccc )
            dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
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
            plt.scatter( x, onesigs, label="(%.2f, %.2f)" %(field_center[0], field_center[1]),\
                        s=50, alpha=.9, c=ccc )
            f_med = interp1d(x, onesigs)
            x2 = np.linspace(min(x), max(x), 1000)
            plt.plot( x2, f_med(x2), '--', alpha=.9, c=ccc)
            
            onesig_errs[band].append(onesigs)
    
    
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
        
        dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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
        dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
    else:
        dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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


def web_plots_UKIDSS( size=1800., ignore_sdss=True, coords=COORDS ):
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
        
        c = gs.catalog( field_center, size, input_coords=input_coords, ignore_sdss=ignore_sdss )
        matches,tmp = gs.identify_matches( input_coords, c.coords )
        
        i_y = gs.FILTER_PARAMS['y'][-1]
        model_matches = []
        obs_matches = []
        model_err = []
        for i,match in enumerate(matches):
            if match >= 0:
                model_matches.append( c.SEDs[match,i_y] )
                obs_matches.append( Y[i][0] )
                model_err.append( c.model_errors[match] )
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
        if ignore_sdss:
            dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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
            dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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
                dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
            else:
                dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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


def compare_to_Stetson( field, ignore_sdss=True, clip=True, colors=['b','g','r','orange'] ):
    input_coords = np.loadtxt('data/'+field+'.pos.clean', skiprows=1)
    data = np.loadtxt('data/'+field+'.pho.clean', skiprows=1)
    obs = data[:,1:] #B,V,R,I
    
    field_center,field_size = gs.find_field( input_coords )
    c = gs.catalog( field_center, max(field_size), input_coords=input_coords, ignore_sdss=ignore_sdss )
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
            dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
        else:
            dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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


FIELDS = ['IC1613', 'L112', 'L113', 'L92', 'L95', 'PG0231']
def web_plots_Stetson( fields=FIELDS, ignore_sdss=True, colors=['b','g','r','orange','grey','yellow'] ):
    
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
        c = gs.catalog( field_center, max(field_size), input_coords=input_coords, ignore_sdss=ignore_sdss )
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
            i_filt = gs.FILTER_PARAMS[filt][-1]
            
            # handle cases with a single bad band (value = -99.)
            obmatch = obs_matches[:,i][ obs_matches[:,i]!=99.999 ]
            momatch = model_matches[:,i_filt][ obs_matches[:,i]!=99.999 ]
            
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
            
            if ignore_sdss:
                dbins = np.array( [0, .25, .5, 1., 1.5, 2., 5.] )
            else:
                dbins = np.array( [0, .5, 1., 2., 5., 7.5, 10.] )
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
    return medians, mads, onesig_errs

