from astropy.coordinates import SkyCoord
import astropy.units as u
import os

maindir = os.getcwd()
outdir = os.path.join(maindir,'results')
#catalog properties
query_cats = {'2mass':'II/246',
              'vvv'  :'II/348'}
#area used for similarsearches. VVV Catalogue is roughly from 295-10deg lon, +-2 deg lat
vvv_leftbot = SkyCoord(['297.0 -2.0'],unit=[u.deg,u.deg],frame='galactic')
vvv_topright= SkyCoord(['295.0  2.0'],unit=[u.deg,u.deg],frame='galactic')
vvvbands = ['Zmag3','Ymag3','Jmag3','Hmag3','Ksmag3']
twomassbands= ['Jmag','Hmag','Kmag']

#query properties
query_radius = 2.83 *u.deg
query_center = SkyCoord(l=295.,b=0.,unit=u.deg,frame='galactic')
comp_radii  = [45.,60.,80.,] *u.arcsec

srcname = 'HD101412'
