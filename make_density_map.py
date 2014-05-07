'''
A script to produce a density map of our catalog coverage.
'''

import get_SEDs as gs
import matplotlib.pyplot as plt
imort matplotlib
from matplotlib import ticker
from astro.allskymap import AllSkyMap
import numpy as np
import pymongo as pm

ra_array = np.linspace( -180, 179, 360 )
dec_array = np.linspace( -89, 89, 179 )
size = 0.25 # means an area of 1/16 of a square degree per query

density_0 = np.empty( (len(ra_array), len(dec_array)) ) #sdss+2mass
density_1 = np.empty( (len(ra_array), len(dec_array)) ) #apass+2mass
density_2 = np.empty( (len(ra_array), len(dec_array)) ) #usnob+2mass
save = True

try:
    DB = pm.MongoClient().photometry
    assert DB.authenticate('phot','ometry')
except:
    raise IOError('cannot connect to database')

try:
    density_0 = np.load('density_0.npy')
    density_1 = np.load('density_1.npy')
    density_2 = np.load('density_2.npy')
except:
    print 'starting calculation'
    for j,dec in enumerate(dec_array):
        print 'total length:',len(ra_array),len(dec_array)
        raddec = np.deg2rad( dec )
        dr = (size/2)/np.cos(raddec) # in degrees
        dd = size/2
        for i,ra in enumerate(ra_array):
            print i,j
            ral = ra-dr
            rah = ra+dr
            decl = dec-dd
            dech = dec+dd
            if ral < -180: ral = ral+360
            if rah > 180: rah = rah-360
            box = { "type" : "Polygon", "coordinates" : [ [ [ral, decl], [rah, decl],
                                                            [rah, dech], [ral, dech],
                                                            [ral, decl] ] ] }
            curs = DB.mass.find( {"coords": {"$geoWithin": {"$geometry":box}} } )
            
            mass, sdss, usnob, apass = [], [], [], []
            count = 0
            while True:
                if not count%100: print '.',
                count +=1
                try:
                    obj = curs.next()
                except StopIteration:
                    break
                try:
                    ttt = [obj[thing] for thing in ['ra','dec','J','J_err','H','H_err','K','K_err']]
                    mass.append( True )
                except:
                    # hit here if this source doesn't have all passbands
                    continue
                if obj['ra'] > 180.0:
                    querycoords = [obj['ra']-360.0, obj['dec']]
                else:
                    querycoords = [obj['ra'], obj['dec']]
                query = {"coords": {"$near": {"$geometry":{"type": "Point", "coordinates": querycoords}, "$maxDistance":30.0}}}
                # now query other catalogs
                #first SDSS
                ocurs = DB['sdss'].find( query )
                try:
                    obj = ocurs.next()
                    ttt = [obj[thing] for thing in ['ra','dec','u','u_err','g','g_err','r','r_err','i','i_err','z','z_err']]
                    sdss.append( True )
                except:
                    sdss.append(None)
                #then APASS
                ocurs = DB['apass'].find( query )
                try:
                    obj = ocurs.next()
                    # allow sources with 3-5 APASS observations, marking the others as 0.0
                    keys = obj.keys()
                    row = [ obj['ra'], obj['dec'] ]
                    for band in ["g'","r'","i'",'B','V']:
                        if (band in keys) and (band+'_err' in keys):
                            # some apass errors are listed as < 0, so here we fix that
                            row += [obj[band], abs(obj[band+"_err"])]
                        elif (band in keys):
                            # for missing errors, adopt an error of 0.15mag
                            row += [obj[band], 0.15]
                        else:
                            # just put placeholders of 0.0 for missing bands (cut later)
                            row += [0.0, 0.0]
                    if len( [r for r in row[::2] if r!=0] ) < 3:
                        # don't bother with any sources that have fewer than 3 observations
                        apass.append(None)
                    else:
                        apass.append( True )
                except:
                    apass.append(None)
                #finally USNOB
                ocurs = DB['usnob'].find( query )
                try:
                    obj = ocurs.next()
                    # use an error of 0.3 for all USNOB observations, and ignore any I-band observations (I, I_sigma = 0.0)
                    ttt = [obj['ra'],obj['dec'],obj['B'],0.3,obj['R'],0.3,0.0,0.0]
                    usnob.append( True )
                except:
                    usnob.append(None)
            # now actually find the matches
            modes = []
            if len(mass) > 0:
                mass = np.array(mass)
                # Go through 2mass objects and assemble all objects present in 2mass and (sdss or apass or usnob)
                #  Preference ranking: 2MASS + (SDSS > APASS > USNOB)
                for iii,obj in enumerate(mass):
                    if (len(sdss) != 0) and (sdss[iii] != None):
                        mode = 0
                    elif (len(apass) != 0) and (apass[iii] != None):
                        mode = 1
                    elif (len(usnob) != 0) and (usnob[iii] != None):
                        mode = 2
                    else:
                        continue
                    modes.append( mode )
            
            # now, use the modes array to actually determine the densities here
            factor = size**-2 # puts the densities in sources per deg^2
            modes = np.array(modes)
            density_0[i,j] = len(modes[modes<1]) * factor
            density_1[i,j] = len(modes[modes<2]) * factor
            density_2[i,j] = len(modes) * factor

# save the results
np.save('density_0.npy',density_0)
np.save('density_1.npy',density_1)
np.save('density_2.npy',density_2)

# now make the plots
X,Y = np.meshgrid(ra_array, dec_array)
values = np.logspace(0, 5, 6)

font = {'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)

plt.figure(figsize=(12,7))
m = AllSkyMap(projection='hammer')
plt.title('SDSS Coverage', y=1.05)
matplotlib.rcParams['text.color'] = '#C21A00'
m.drawparallels(np.arange(-75,76,15), linewidth=0.5, dashes=[1,2],
                labels=[1,0,0,0], fontsize=12)
m.drawmeridians(np.arange(-150,151,30), linewidth=0.5, dashes=[1,2])
# Label a subset of meridians.
lons = np.arange(-150,151,30)
m.label_meridians(lons, fontsize=12, vnudge=1,
                  halign='left', hnudge=-1)  # hnudge<0 shifts to right
matplotlib.rcParams['text.color'] = 'k'
m.contourf(X.T, Y.T, density_0, values, latlon=True, cmap=plt.cm.bone_r, locator=ticker.LogLocator())
cb = m.colorbar()
cb.set_label(r'Sources per deg$^2$')
if save: plt.savefig('mode1.png', transparent=True)

plt.figure(figsize=(12,7))
m = AllSkyMap(projection='hammer')
plt.title('SDSS and APASS Coverage', y=1.05)
matplotlib.rcParams['text.color'] = '#C21A00'
m.drawparallels(np.arange(-75,76,15), linewidth=0.5, dashes=[1,2],
                labels=[1,0,0,0], fontsize=12)
m.drawmeridians(np.arange(-150,151,30), linewidth=0.5, dashes=[1,2])
# Label a subset of meridians.
lons = np.arange(-150,151,30)
m.label_meridians(lons, fontsize=12, vnudge=1,
                  halign='left', hnudge=-1)  # hnudge<0 shifts to right
matplotlib.rcParams['text.color'] = 'k'
m.contourf(X.T, Y.T, density_1, values, latlon=True, cmap=plt.cm.bone_r, locator=ticker.LogLocator())
cb = m.colorbar()
cb.set_label(r'Sources per deg$^2$')
if save: plt.savefig('mode12.png', transparent=True)

plt.figure(figsize=(12,7))
m = AllSkyMap(projection='hammer')
plt.title('SDSS, APASS, and USNOB Coverage', y=1.05)
matplotlib.rcParams['text.color'] = '#C21A00'
m.drawparallels(np.arange(-75,76,15), linewidth=0.5, dashes=[1,2],
                labels=[1,0,0,0], fontsize=12)
m.drawmeridians(np.arange(-150,151,30), linewidth=0.5, dashes=[1,2])
# Label a subset of meridians.
lons = np.arange(-150,151,30)
m.label_meridians(lons, fontsize=12, vnudge=1,
                  halign='left', hnudge=-1)  # hnudge<0 shifts to right
matplotlib.rcParams['text.color'] = 'k'
m.contourf(X.T, Y.T, density_2, values, latlon=True, cmap=plt.cm.bone_r, locator=ticker.LogLocator())
cb = m.colorbar()
cb.set_label(r'Sources per deg$^2$')
if save: plt.savefig('mode123.png', transparent=True)

plt.show()
