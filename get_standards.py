'''
find a set of stars in both NOMAD and SDSS
'''

from find_photometry import *
from subprocess import Popen, PIPE

center = (314.128304, -06.070742)
width = 600. #arcseconds
request = 'findnomad1 -c "{} {}" -bs {}'.format(center[0], center[1], width)
out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
o,e = out.communicate()
# parse the response
nomad_objects = parse_nomad1(o)

keepers = []
for obj in nomad_objects:
    if all(obj[2:]):
        # select objects with all bands included
        keepers.append(obj)

# now do sdss
request = 'findsdss8 -c "{} {}" -bs {} -e0,'.format(center[0], center[1], width)
out = Popen(request, shell=True, stdout=PIPE, stderr=PIPE)
o,e = out.communicate()
# parse the response
sdss_objects = parse_sdss8(o)

# match them
final_coords = []
for star in keepers:
    diffs = [((star[0]-other[0])**2 + (star[1]-other[1])**2)**.5 for other in sdss_objects[:,:2]]
    if min(diffs) < .5*2.778e-4:
        # a match is if they're within .5" of each other
        final_coords.append( (star[0], star[1]) )

# final_coords is now a list of tuples of coordinates, each
# for a star with a match in both NOMAD and SDSS
np.save( open('matches.npy','w'), final_coords )