import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
from astropy.coordinates import Angle,SkyCoord
import astropy.units as u
from astropy.table import Table
from starclass import Star
from parameters import *
import os
import ipdb

star = Star(srcname)
star.coordinates = 'auto'



def do():
    '''Does everything except downloading the data. Use down_catalogues for that.'''
    vvvcat,twomasscat = read_cat()
    del vvvcat
    references = []
    results = []
    cencoords = make_comp_grid()
    print('Processing the 2Mass referencepoint')
    for radius in comp_radii:
        for wave in twomassbands:
            references.append(['2mass',star.coordinates,radius,wave] +
                              [x[0] for x in list(make_cumhists(twomasscat,star.coordinates,
                                                                radius,[wave,]))])
    print('Processing the 2Masssurvey')
    for ii,center in enumerate(cencoords.icrs):
        print('Processing center {} ({} of {})'.format(center,ii,len(cencoords)-1))
        for radius in comp_radii:
            for wave in twomassbands:
                results.append(['2mass',center,radius,wave] +
                               [x[0] for x in list(make_cumhists(twomasscat,SkyCoord([center,]),
                                                                 radius,[wave,]))])

    print('Done with cumhists. Plotting them now')
    wb2col = {'Jmag':'blue',
              'Hmag':'green',
              'Kmag':'red'}
    bins    = references[-1][-1][1:]
    cummax = np.max([x[-2] for x in references])
    fig = plt.figure(figsize=(5*len(comp_radii), 5*len(cencoords)))
    outergrid = gridspec.GridSpec(len(cencoords), len(comp_radii), wspace=0.2,hspace=0.2)
    for icen,center in enumerate(cencoords.icrs):
        for irad,radius in enumerate(comp_radii):
            iplot = irad*len(cencoords)+icen
            innergrid = gridspec.GridSpecFromSubplotSpec(3, 1,
                                subplot_spec=outergrid[icen,irad], wspace=0.1, hspace=0.1)
            ax1 = plt.Subplot(fig,innergrid[:-1])
            ax2 = plt.Subplot(fig,innergrid[-1])
            for wave in twomassbands:
                refvals = [x[-2] for x in references if ((x[2] == radius) and
                                                         (x[3] == wave))][0]
                scivals = [x[-2] for x in results if ((x[2] == radius) and 
                                                      (x[1].separation(center) < 0.01*u.arcsec) and
                                                      (x[3] == wave))][0]
                ax1.plot(bins,refvals,ls='--',c=wb2col[wave],label=wave)
                ax1.plot(bins,scivals,ls='-' ,c=wb2col[wave])
                ax1.set_ylim(0,cummax)
                ax2.set_ylim(-cummax/5,cummax/5)
                ax2.plot(bins,refvals-scivals,c=wb2col[wave])
            ax1.legend(loc='upper left')
            #ax1.set_title('2MASS, r={}, cen={}'.format(radius,center.galactic))
            #make the labels where needed
            import ipdb;ipdb.set_trace()
            if irad == 0:
                ax1.set_ylabel('Cumulative nr of sources')
                ax2.set_ylabel('Residuals')
                plt.subplot(outergrid[iplot]).set_ylabel('{}'.format(center.galactic))
            else:
                ax1.set_yticks([])
                ax2.set_yticks([])
            if icen == len(cencoords):
                plt.subplot(outergrid[iplot]).set_xlabel('mag')
            else:
                ax1.set_xticklabels([])
                ax2.set_xticklabels([])
            
            fig.add_subplot(ax1)
            fig.add_subplot(ax2,sharex=ax2)
    fp_plot = os.path.join(outdir,'sourcedensity_overview.pdf')
    fig.savefig(fp_plot)
    print('Saved plot to {}'.format(fp_plot))
    plt.close('all')
                        
    import ipdb;ipdb.set_trace()


            
                                                     


def make_cumhists(cat,center,radius,wavebands,magrange=[5,14.5]):
    '''Make a magnitude histogram around the center with a given radius for each waveband'''
    nbins = 200
    idz1,idz2,sep2d,sep3d = cat['skycoord'].search_around_sky(center,radius)
    
    small_cat = cat[idz2]
    hists = []
    for wave in wavebands:
        values, bins = np.histogram(small_cat[wave],bins=nbins,range=magrange)
        hists.append(np.cumsum(values))
    return hists,[bins,]
    
def read_cat():
    '''Read the catalogues. Mainly convert coordinates to the ones given in column
    skycoord.'''
    vvvcat   = Table.read(os.path.join(maindir,'vvv_catalogue.csv'),format='ascii.csv')
    twomasscat = Table.read(os.path.join(maindir,'2mass_catalogue.csv'),format='ascii.csv')
    vvvcat['skycoord']    = SkyCoord(vvvcat['RAJ2000'],vvvcat['DEJ2000'],frame='icrs',unit=u.deg)
    twomasscat['skycoord']= SkyCoord(twomasscat['RAJ2000'],twomasscat['DEJ2000'],
                                     frame='icrs',unit=u.deg)
    return vvvcat,twomasscat

def make_comp_grid():
    '''Makes the center points of the comparison. Returns a SkyCoord object with the
    centrag points'''
    max_rad = np.max(comp_radii)
    nlat = 2
    nlon = 2
    #make the longitude and latitude vectors in galactic coords
    veclon = np.linspace(vvv_topright.l + max_rad, vvv_leftbot.l - max_rad, nlon)
    veclat = np.linspace(vvv_leftbot.b + max_rad, vvv_topright.b - max_rad, nlat)
    
    vlon = np.meshgrid(veclon,veclat)[0].flatten()
    vlat = np.meshgrid(veclon,veclat)[1].flatten()
    
    cencoords = SkyCoord(vlon,vlat,unit=u.deg,frame='galactic')
    
    return cencoords
    

def down_catalogues():
    '''Download and save the catalogs'''
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    for cat in query_cats.keys():
        res = Vizier.query_region(query_center,catalog=query_cats[cat],
                               radius=query_radius)
        cat_pname = os.path.join(maindir,cat+'_catalogue.csv')
        res[0].write(cat_pname,
                       format='ascii.csv',
                       delimiter=',',overwrite=True)
        print('Saved cat {} ({} rows) to {}'.format(cat,len(res[0]),cat_pname))
        del res
    import ipdb;ipdb.set_trace()
        
