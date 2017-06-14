import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import gridspec
from astropy.coordinates import Angle,SkyCoord
import astropy.units as u
from astropy.table import Table
from starclass import Star
from parameters import *
import pickle
import os
import ipdb

star = Star(srcname)
star.coordinates = 'auto'



def do():
    '''Does everything except downloading the data. Use down_catalogues for that.'''
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    twomasscat = read_cat(vvv=False)
    references = []
    results = []#in order catalogue,coords,raudius,wavelength,[histvalues,],[histbins,]
    if grid == 'auto':
        print('Auto-making the grid')
        cencoords = make_comp_grid()
    elif grid == 'manual':
        print('Using manually center points')
        cencoords = load_comp_grid()
    else:
        raise ValueError('Unknown keyword {} for grid. Select auto or manual'.format(grid))
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
    pn_results    = os.path.join(outdir,'results.py')
    pn_references = os.path.join(outdir,'references.py')
    pickle.dump(results,open(pn_results,'wb'))
    pickle.dump(references,open(pn_references,'wb'))
    print('Done with cumhists. Saved them to {} , {}'.format(pn_results,pn_references))
    if len(cencoords) > 100:
        plot = False
    else:
        plot=True
    t_maxdiff = plot_bins(results,references,cencoords,plot=plot)
    #find the best result in a min of max diff sense
    t_diff_grouped = t_maxdiff.group_by(['cen_lon','cen_lat','radius'])
    t_diff_max = t_diff_grouped.groups.aggregate(np.max)
    t_diff_max.sort(['radius','maxdiff'])
    #mindiffglob = np.min(t_diff_max['maxdiff'])
    #t_bestcens  = t_diff_max.sort([np.where(t_diff_max['maxdiff'] == mindiffglob)]
    #save the results
    pn_tgrouped = os.path.join(outdir,'centerprops_grouped.csv')
    pn_tall     = os.path.join(outdir,'centerprops_all.csv')
    #pn_tbest= os.path.join('centerprops_best.csv')
    i_largest_radius = len(cencoords)*(len(comp_radii)-1)
    print('The 10 best centers seem {}.\n Saved them to {}'.format(\
          t_diff_max[i_largest_radius:i_largest_radius+10],pn_tgrouped))
    t_diff_grouped.write(pn_tall,format='ascii.csv',overwrite=True)
    t_diff_max.write(pn_tgrouped,format='ascii.csv',overwrite=True)
    #t_bestcens.write(pn_tbest,format='ascii.csv')

    print('DONE with get_background.do')


def plot_bins(results,references,cencoords,plot=True):
    '''make a grid plot in one file. Returns a table with simple statistics.'''
    import seaborn as sn
    wb2col = {'Jmag':'blue',
              'Hmag':'green',
              'Kmag':'red'}
    #store stuff in a table. Fill and delete dummy values for right column properties
    t_maxdiff = Table([[np.nan],[np.nan],[np.nan,],['Dummy_waveband',],[-999]],
                      names=['cen_lon','cen_lat','radius','waveband','maxdiff'])
    t_maxdiff = t_maxdiff[:0]
    bins    = references[-1][-1][1:]
    cummax = np.max([x[-2] for x in results])
    if plot: 
        fig = plt.figure(figsize=(5*len(comp_radii), 5*len(cencoords)))
        outergrid = gridspec.GridSpec(len(cencoords), len(comp_radii), wspace=0.12,hspace=0.12)
    for icen,center in enumerate(cencoords.icrs):
        if plot:
            print('Plotting and getting statistics for Position {} of {}'.format(icen,
                                                                        len(cencoords.icrs)-1))
        else:
            print('Not plotting. Getting statistics for Position {} of {}'.format(icen,
                                                                        len(cencoords.icrs)-1))
        for irad,radius in enumerate(comp_radii):
            if plot:
                innergrid = gridspec.GridSpecFromSubplotSpec(3, 1,
                                subplot_spec=outergrid[icen,irad], wspace=0.1, hspace=0.1)
                ax1 = plt.Subplot(fig,innergrid[:-1,0])
                ax2 = plt.Subplot(fig,innergrid[-1,0])
            for wave in twomassbands:
                refvals = [x[-2] for x in references if ((x[2] == radius) and
                                                         (x[3] == wave))][0]
                scivals = [x[-2] for x in results if ((x[2] == radius) and 
                                                      (x[1].separation(center) < 0.01*u.arcsec) and
                                                      (x[3] == wave))][0]
                if plot:
                    ax1.plot(bins,refvals,ls='--',c=wb2col[wave],label=wave)
                    ax1.plot(bins,scivals,ls='-' ,c=wb2col[wave])
                    ax1.set_ylim(0,cummax)
                    ax2.set_ylim(-cummax/8,cummax/8)
                    ax2.plot(bins,refvals-scivals,c=wb2col[wave])
                t_maxdiff.add_row([center.galactic.l,center.galactic.b,
                                   radius,wave,int(np.max(np.abs(refvals-scivals)))])
            if plot:
                ax1.legend(loc='upper left')
                #ax1.annotate('icen,irad: {},{}'.format(icen,irad),xy=(10,20))
                #ax1.set_title('2MASS, r={}, cen={}'.format(radius,center.galactic))
                #make the labels where needed
                if irad == 0:
                    ax1.set_ylabel('lon={0:.4f}, lat={1:4f} \n\nCumulative nr of sources'.format(center.galactic.l.value,center.galactic.b.value))
                    ax2.set_ylabel('Residuals')
                    #plt.subplot(outergrid[icen,irad]).set_ylabel('lon={10.6f}, lat={8.6f}\n\n\n'.format(center.galactic.l.value,center.galactic.b.value))
                    #plt.subplot(outergrid[icen,irad]).set_yticks([])
                    #plt.subplot(outergrid[icen,irad]).set_xticks([])
                else:
                    ax1.set_yticklabels([])
                    ax1.set_ylabel('')
                    ax2.set_ylabel('')
                    ax2.set_yticklabels([])
                if icen == len(cencoords)-1:
                    ax2.set_xlabel('mag')
                    ax1.set_xticklabels([])
                elif icen == 0:
                    ax1.set_title('search radius = {}'.format(radius))
                    ax1.set_xticklabels([])
                    #plt.subplot(outergrid[icen,irad]).set_ylabel('search radius = {}'.format(radius))
                    #plt.subplot(outergrid[icen,irad]).set_yticks([])
                    #plt.subplot(outergrid[icen,irad]).set_xticks([])
                else:
                    ax2.set_xlabel('')
                    ax1.set_xticklabels([])
                    ax2.set_xticklabels([])
            if plot:
                fig.add_subplot(ax1)
                fig.add_subplot(ax2,sharex=ax2)
    if plot:
        fp_plot = os.path.join(outdir,'sourcedensity_overview.pdf')
        try:
            fig.savefig(fp_plot)
        except:
            print('Somehow could not save the figure.')
        print('Saved plot to {}'.format(fp_plot))
        plt.close('all')

    return t_maxdiff

            
                                                     


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
    
def read_cat(vvv=True,twomass=True):
    '''Read the catalogues. Mainly convert coordinates to the ones given in column
    skycoord.'''
    if vvv:
        print('Reading vvvcat')
        vvvcat   = Table.read(os.path.join(maindir,'vvv_catalogue.csv'),format='ascii.csv')
        vvvcat['skycoord']    = SkyCoord(vvvcat['RAJ2000'],vvvcat['DEJ2000'],frame='icrs',unit=u.deg)
    if twomass:
        print('Reading 2masscat')
        twomasscat = Table.read(os.path.join(maindir,'2mass_catalogue.csv'),format='ascii.csv')
        twomasscat['skycoord']= SkyCoord(twomasscat['RAJ2000'],twomasscat['DEJ2000'],
                                     frame='icrs',unit=u.deg)
    if vvv and twomass:
        return vvvcat,twomasscat
    elif vvv:
        return vvvcat
    else:
        return towmasscat

def make_comp_grid():
    '''Makes the center points of the comparison. Returns a SkyCoord object with the
    centrag points'''
    max_rad = np.max(comp_radii)
    #make the longitude and latitude vectors in galactic coords
    veclon = np.linspace(vvv_topright.l + max_rad, vvv_leftbot.l - max_rad, nlon)
    veclat = np.linspace(vvv_leftbot.b + max_rad, vvv_topright.b - max_rad, nlat)
    
    vlon = np.meshgrid(veclon,veclat)[0].flatten()
    vlat = np.meshgrid(veclon,veclat)[1].flatten()
    
    cencoords = SkyCoord(vlon,vlat,unit=u.deg,frame='galactic')
    
    return cencoords
    
def load_comp_grid():
    '''Uses manually entered grid. Counterpart to make_comp_grid'''
    return manual_center_points
    
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
        
