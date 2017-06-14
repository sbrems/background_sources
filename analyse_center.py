from astropy.coordinates import SkyCoord
import time
import numpy as np
from astropy.table import Table,vstack
import astropy.units as u
import os
import fileinput

from parameters import *
from get_background import read_cat

def do():
    if not os.path.exists(vosadir):
        os.makedirs(vosadir)
    vvvcat = read_cat(twomass=False)

    idz1,idz2,sep2d,sep3d = vvvcat['skycoord'].search_around_sky(best_centers,best_radius)
    for icen,center in enumerate(best_centers):
        small_cat = vvvcat[idz2[np.where(idz1 == icen)]]
        make_vosa_input(small_cat,center)
        
def make_vosa_input(table,center):
    print('Making the VOSA compatible input table')
    vosatable = Table([['longdummystarname',],[-999.,],[-999.,],['---',],['---',],
                       ['longdummyfiltername',],
                       [-999.,],[-999.,],['nofit',],['longcoomentdummy',]],
                      names=['starname','RA','DEC','DIS','Av','Filter',
                             'Flux','Error','PntOpts','ObjOpts'])
    vosatable = vosatable[:0]
    nstars = len(table)
    for vvvband,vosabandname in vvv2vosa.items():
        temptable = Table(table['iauname','RAJ2000','DEJ2000',vvvband,'e_'+vvvband],
                          names=['starname','RA','DEC','Flux','Error'])
        temptable['Filter'] = vosabandname
        temptable['DIS'] = '---'
        temptable['Av'] = '---'
        temptable['PntOpts'] = '---'
        temptable['ObjOpts'] = 'Av:0.0/3.0'
        vosatable = vstack([vosatable,temptable])
    #now add the NaCo band. Use Ksband as proxy
    temptable = Table(table['iauname','RAJ2000','DEJ2000','Ksmag3','e_Ksmag3'],
                      names=['starname','RA','DEC','Flux','Error'])
    temptable['Filter'] = 'Paranal/NACO.Lp'
    temptable['DIS'] = '---'
    temptable['Av'] = '---'
    temptable['PntOpts'] = 'nofit'
    temptable['ObjOpts'] = 'Av:0.0/3.0'
    vosatable = vstack([vosatable,temptable])
    pntable = os.path.join(vosadir,'vosainputtable_cen{0:.4f}_{1:4f}.csv'.format(\
                                   center.galactic.l.value,center.galactic.b.value))
    vosatable['starname'] = [x.replace(" ","_") for x in vosatable['starname']]
    #comment what was done
    vosatable.meta['comments'] = ['Stars fromm VVV for center lon/lat:{}/{} (gal) and radius {}'.format(\
                                    center.galactic.l.value,center.galactic.b.value,best_radius)]
    #remove invalid/empty entries
    #vosatable.remove_rows(vosatable['Flux'].mask.nonzero())
    vosatable.write(pntable,
                    names=['starname','RA','DEC','DIS','Av','Filter','Flux',
                           'Error','PntOpts','ObjOpts'],
                    format='ascii.commented_header',delimiter=' ',
                    fill_values='---',overwrite=True)
    #replace "" with --- 

    with fileinput.FileInput(pntable, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace('""', '---'), end='')
    print('Saved Vosatable ({} rows) to {}'.format(len(vosatable),pntable))
    
    import ipdb;ipdb.set_trace()


#define some functions for the table savings
def process2vosa(val):
    val = val.strip()
    val.replace(" ","_")
    val.replace('"','')
