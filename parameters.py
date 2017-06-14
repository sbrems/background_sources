from astropy.coordinates import SkyCoord
import time
import astropy.units as u
import os

maindir = os.getcwd()
outdir = os.path.join(maindir,'results',time.strftime("%Y-%m-%d_%H_%M"))
vosadir = os.path.join(maindir,'results','vosa')
#catalog properties
query_cats = {'2mass':'II/246',
              'vvv'  :'II/348'}
#area used for similarsearches. VVV Catalogue is roughly from 295-10deg lon, +-2 deg lat
vvv_leftbot = SkyCoord(['297.0 -2.0'],unit=[u.deg,u.deg],frame='galactic')
vvv_topright= SkyCoord(['295.0  2.0'],unit=[u.deg,u.deg],frame='galactic')
nlat = 20
nlon = 20
vvv2vosa = {'Zmag3':'Paranal/VISTA.Z',
            'Ymag3':'Paranal/VISTA.Y',
            'Jmag3':'Paranal/VISTA.J',
            'Hmag3':'Paranal/VISTA.H',
            'Ksmag3':'Paranal/VISTA.Ks'}
twomassbands= ['Jmag','Hmag','Kmag']

#query properties
query_radius = 3. *u.deg
query_center = SkyCoord(l=295.,b=0.,unit=u.deg,frame='galactic')
comp_radii  = [90.,180.,] *u.arcsec

srcname = 'HD101412'


#if you want to use manual center points, enter them here and set keyword in get_background
grid = 'auto' #'auto' or 'manual'
manual_center_points = SkyCoord([296.85,295.95,296.35,295.95,295.15,295.45,
                                 295.45,296.15,295.85,296.55,#10 good centers
                                 296.85,296.85,296.75,296.55,296.75,],#5 bad centers
                                [1.744736,1.5394736,1.539473,1.95,1.33421,1.744736,
                                 1.5394736,1.7447368,1.5394736,1.3342105,
                                 -0.5131578,-0.102631578,-0.102631578,-0.102631578,-0.51315789],
                                unit=u.deg, frame='galactic')
best_centers = SkyCoord([295.950],
                        [1.539474],
                        unit=u.deg, frame='galactic')
best_radius = 180*u.arcsec
