#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.mlab import griddata
import numpy.ma as ma
import scipy.ndimage.filters 
# from skimage.filters import roberts, sobel, scharr, prewitt
from scipy.ndimage import convolve
import matplotlib
from lasif import colors
import os
import obspy.geodetics
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from pyproj import Geod
import random
import copy
import pycpt
import pyasdf
import math
from subprocess import call

lon_diff_weight_2 = np.array([[1., 0., -1.]])/2.
lat_diff_weight_2 = lon_diff_weight_2.T
lon_diff_weight_4 = np.array([[-1., 8., 0., -8., 1.]])/12.
lat_diff_weight_4 = lon_diff_weight_4.T
lon_diff_weight_6 = np.array([[1./60., 	-3./20.,  3./4.,  0., -3./4., 3./20.,  -1./60.]])
lat_diff_weight_6 = lon_diff_weight_6.T

lon_diff2_weight_2 = np.array([[1., -2., 1.]])
lat_diff2_weight_2 = lon_diff2_weight_2.T
lon_diff2_weight_4 = np.array([[-1., 16., -30., 16., -1.]])/12.
lat_diff2_weight_4 = lon_diff2_weight_4.T
lon_diff2_weight_6 = np.array([[1./90., 	-3./20.,  3./2.,  -49./18., 3./2., -3./20.,  1./90.]])
lat_diff2_weight_6 = lon_diff2_weight_6.T

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

class model2d(object):
    """
    An object to analyze/plot 2D spherical model on Earth
    ===========================================================================
    Parameters:
    dlon, dlat      - grid interval
    Nlon, Nlat      - grid number in longitude, latitude 
    lonArr, latArr  - arrays for grid location
    fieldtype       - field type (Tph, Tgr, Amp)
    ---------------------------------------------------------------------------
    Note: meshgrid's default indexing is 'xy', which means:
    lons, lats = np.meshgrid[lon, lat]
    in lons[i, j] or lats[i, j],  i->lat, j->lon
    ===========================================================================
    """
    def __init__(self, infname=None, minlon=None, maxlon=None, dlon=0.5, minlat=None, maxlat=None, dlat=0.5):
        if infname==None and (minlon==None or maxlon==None or minlat==None or maxlat==None):
            raise ValueError('Fail in initializing parameters!')
        if infname != None:
            try:
                Inarray=np.loadtxt(infname)
                self.lonArrIn=Inarray[:,0]
                self.latArrIn=Inarray[:,1]
                minlon = self.lonArrIn.min(); maxlon = self.lonArrIn.max()
                minlat = self.latArrIn.min(); maxlat = self.latArrIn.max()
            except:
                self.load(infname=infname)
                return
        self.Nlon   = int(round((maxlon-minlon)/dlon)+1)
        self.Nlat   = int(round((maxlat-minlat)/dlat)+1)
        self.dlon   = dlon
        self.dlat   = dlat
        self.lon    = np.arange(self.Nlon)*self.dlon+minlon
        self.lat    = np.arange(self.Nlat)*self.dlat+minlat
        self.lonArr, self.latArr = np.meshgrid(self.lon, self.lat)
        self.minlon = minlon
        self.maxlon = self.lon.max()
        self.minlat = minlat
        self.maxlat = self.lat.max()
        # self._get_dlon_dlat_km()
        # self.fieldtype=fieldtype
        self.Zarr   = np.zeros((self.Nlat, self.Nlon))
        self.mask   = np.zeros((self.Nlat, self.Nlon), dtype=bool)
        return
    
    
    def __add__(self, infield):
        if isinstance(infield, model2d):
            out = self.copy()
            self.ma2np(); infield.ma2np()
            out.Zarr = self.Zarr + infield.Zarr
            out.mask = ~((~self.mask)*(~infield.mask))
            return out
        else:
            raise TypeError
        
    def __sub__(self, infield):
        if isinstance(infield, model2d):
            out = self.copy()
            self.ma2np(); infield.ma2np()
            out.Zarr = self.Zarr - infield.Zarr
            out.mask = ~((~self.mask)*(~infield.mask))
            return out
        else:
            raise TypeError
        
    
    def copy(self): return copy.deepcopy(self)
    
    def _get_dlon_dlat_km(self):
        """Get longitude and latitude interval in km
        """
        self.dlon_km=np.array([])
        self.dlat_km=np.array([])
        for lat in self.lat:
            dist_lon, az, baz = obspy.geodetics.gps2dist_azimuth(lat, 0., lat, self.dlon)
            dist_lat, az, baz = obspy.geodetics.gps2dist_azimuth(lat, 0., lat+self.dlat, 0.)
            self.dlon_km=np.append(self.dlon_km, dist_lon/1000.)
            self.dlat_km=np.append(self.dlat_km, dist_lat/1000.)
        self.dlon_kmArr=(np.tile(self.dlon_km, self.Nlon).reshape(self.Nlon, self.Nlat)).T
        self.dlat_kmArr=(np.tile(self.dlat_km, self.Nlon).reshape(self.Nlon, self.Nlat)).T
        return
    
    def _get_basemap(self, projection='lambert', geopolygons=None, resolution='i'):
        """Plot data with contour
        """
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        lat_centre = (self.maxlat+self.minlat)/2.0
        lon_centre = (self.maxlon+self.minlon)/2.0
        if projection=='merc':
            m=Basemap(projection='merc', llcrnrlat=self.minlat-5., urcrnrlat=self.maxlat+5., llcrnrlon=self.minlon-5.,
                      urcrnrlon=self.maxlon+5., lat_ts=20, resolution=resolution)
            m.drawparallels(np.arange(-80.0,80.0,5.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,5.0), labels=[1,0,0,1])
            m.drawstates(color='g', linewidth=2.)
        elif projection=='global':
            m=Basemap(projection='ortho',lon_0=lon_centre, lat_0=lat_centre, resolution=resolution)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,10.0), labels=[1,0,0,1])
        elif projection=='regional_ortho':
            m1 = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution='l')
            m = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution=resolution,\
                llcrnrx=0., llcrnry=0., urcrnrx=m1.urcrnrx/mapfactor, urcrnry=m1.urcrnry/3.5)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,0],  linewidth=2,  fontsize=20)
            m.drawmeridians(np.arange(-170.0,170.0,10.0),  linewidth=2)
        elif projection=='lambert':
            distEW, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.minlat, self.maxlon) # distance is in m
            distNS, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.maxlat+2., self.minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS, rsphere=(6378137.00,6356752.3142), resolution='l', projection='lcc',\
                lat_1=self.minlat, lat_2=self.maxlat, lon_0=lon_centre, lat_0=lat_centre+1)
            m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=1, dashes=[2,2], labels=[1,1,0,0], fontsize=15)
            m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=1, dashes=[2,2], labels=[0,0,1,0], fontsize=15)
        m.drawcoastlines(linewidth=1.0)
        m.drawcountries(linewidth=1.)
        m.drawstates()
        m.fillcontinents(lake_color='#99ffff',zorder=0.2)
        # m.drawlsmask(land_color='0.8', ocean_color='#99ffff')
        m.drawmapboundary(fill_color="white")
        try: geopolygons.PlotPolygon(inbasemap=m)
        except: pass
        return m
    
    def read_irregular(self, infname, interpolate=True):
        """read 2d model file
        """
        try:
            Inarray=np.loadtxt(infname)
            with open(infname) as f:
                inline = f.readline()
                if inline.split()[0] =='#':
                    evlostr = inline.split()[1]
                    evlastr = inline.split()[2]
                    if evlostr.split('=')[0] =='enx':
                        self.evlo = float(evlostr.split('=')[1])
                    if evlastr.split('=')[0] =='eny':
                        self.evla = float(evlastr.split('=')[1])
        except:
            Inarray=np.load(infname)
        self.lonArrIn=Inarray[:,0]
        self.latArrIn=Inarray[:,1]
        self.ZarrIn=Inarray[:,2]
        if interpolate: self.interp_surface(workingdir='.', outfname='temp_model2d', tension=0.0)
        return
    
    def read(self, infname):
        """
        Read txt velocity model file
        """
        InArr   = np.loadtxt(infname)
        inlon   = InArr[:,0]
        inlat   = InArr[:,1]
        inZ     = InArr[:,2]
        self.mask = ~self.mask
        for i in xrange(inlon.size):
            lon=inlon[i]
            if lon < 0: lon+=360
            lat=inlat[i]
            index = np.where((self.lonArr==lon)*(self.latArr==lat))
            if inZ[i]==0 or math.isnan(inZ[i]): continue
            self.mask[index[0], index[1]]=False
            self.Zarr[index[0], index[1]]=inZ[i]
        return
    
    def save(self, outfname):
        self.ma2np()
        np.savez(outfname, self.Zarr, self.mask, self.lon, self.lat)
        return
    
    def load(self, infname):
        indata      = np.load(infname)
        self.Zarr   = indata['arr_0']
        self.mask   = indata['arr_1']
        self.lon    = indata['arr_2']
        self.lat    = indata['arr_3']
        self.Nlon   = self.lon.size
        self.Nlat   = self.lat.size
        self.dlon   = abs(self.lon[0]-self.lon[1])
        self.dlat   = abs(self.lat[0]-self.lat[1])
        self.lonArr, self.latArr = np.meshgrid(self.lon, self.lat)
        self.minlon = self.lon.min()
        self.maxlon = self.lon.max()
        self.minlat = self.lat.min()
        self.maxlat = self.lat.max()
        return
    
    def interp_surface(self, workingdir='.', outfname='temp_model2d', tension=0.0, resetmask=True):
        """Interpolate input data to grid point with gmt surface command
        =======================================================================================
        Input Parameters:
        workingdir  - working directory
        outfname    - output file name for interpolation
        tension     - input tension for gmt surface(0.0-1.0)
        =======================================================================================
        """
        if not os.path.isdir(workingdir):
            os.makedirs(workingdir)
        try:
            outlonArr   = self.lonArrIn
            outlatArr   = self.latArrIn
            outZArr     = self.ZarrIn
        except:
            outind      = self.mask == False
            outlonArr   = self.lonArr[outind]
            outlatArr   = self.latArr[outind]
            outZArr     = self.Zarr[outind]
        OutArr=np.append(outlonArr, outlatArr)
        OutArr=np.append(OutArr, outZArr)
        OutArr=OutArr.reshape(3, outlonArr.size)
        OutArr=OutArr.T
        np.savetxt(workingdir+'/'+outfname, OutArr, fmt='%g')
        fnameHD=workingdir+'/'+outfname+'.HD'
        tempGMT=workingdir+'/'+outfname+'_GMT.sh'
        grdfile=workingdir+'/'+outfname+'.grd'
        with open(tempGMT,'wb') as f:
            REG='-R'+str(self.minlon)+'/'+str(self.maxlon)+'/'+str(self.minlat)+'/'+str(self.maxlat)
            f.writelines('gmtset MAP_FRAME_TYPE fancy \n')
            f.writelines('surface %s -T%g -G%s -I%g %s \n' %( workingdir+'/'+outfname, tension, grdfile, self.dlon, REG ))
            f.writelines('grd2xyz %s %s > %s \n' %( grdfile, REG, fnameHD ))
        call(['bash', tempGMT])
        os.remove(grdfile)
        os.remove(tempGMT)
        Inarray=np.loadtxt(fnameHD)
        ZarrIn=Inarray[:,2]
        self.Zarr=(ZarrIn.reshape(self.Nlat, self.Nlon))[::-1, :]
        if resetmask: self.mask = np.zeros((self.Nlat, self.Nlon), dtype=bool)
        return
    
    def np2ma(self):
        """Convert data array to masked array according to self.mask array.
        """
        self.Zarr       = ma.masked_array(self.Zarr, mask=np.zeros(self.mask.shape) )
        self.Zarr.mask  = self.mask
        return
    
    def ma2np(self):
        """Convert the maksed data array to numpy data array
        """
        try:
            self.mask   = self.Zarr.mask
            self.Zarr   = ma.getdata(self.Zarr)
        except: print 'Data array is already numpy array'
        return
    
    def smooth(self, sigma):
        z_filtered=self.Zarr.copy()
        wZrr = self.Zarr * (~self.mask)
        for iteration in xrange(int(sigma)):
            for i in np.arange(1,self.Nlat-1):
                for j in np.arange(1,self.Nlon-1):
                    if self.mask[i, j]: continue
                    weight = 5 - (self.mask[i,j]+self.mask[i+1,j]+self.mask[i-1,j]+self.mask[i,j+1]+self.mask[i,j-1])
                    z_filtered[i,j]=(wZrr[i,j]+wZrr[i+1,j]+wZrr[i-1,j]+wZrr[i,j+1]+wZrr[i,j-1])/weight
        self.Zarr=z_filtered
        return
    
    def hist(self, vmin=-0.2, vmax=0.2, label='', unit='(km/sec)'):
        from matplotlib.ticker import FuncFormatter
        ax=plt.subplot()
        valid_data = self.Zarr[(~self.mask)]
        # the histogram of the data
        # n, bins, patches = plt.hist(valid_data, bins=50, normed=1, facecolor='red', alpha=0.75,
        #                             weights=np.zeros_like(valid_data) + 100. / valid_data.size)
        n, bins, patches = plt.hist(valid_data, bins=100, facecolor='red', alpha=0.75,
                                    weights=np.zeros_like(valid_data) + 100. / valid_data.size)
        plt.ylabel('percentage (%)', fontsize=15)
        plt.xlabel(label+unit, fontsize=15)
        # plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
        maxvalue = max(valid_data.max(), - valid_data.min())
        if vmin==None or vmax==None: vmin=-maxvalue; vmax=maxvalue
        plt.xlim([vmin-0.001, vmax+0.001])
        # plt.ylim([0, 50])
        plt.grid(True)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        plt.show()
        

    def plot(self, label='', unit='(km/sec)', projection='lambert', cmap='cv', contour=False, geopolygons=None, showfig=True, vmin=None, vmax=None):
        """
        """
        
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        x, y=m(self.lonArr, self.latArr)
        if cmap == 'lasif':
            cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        elif cmap == 'cv':
            try: cmap=pycpt.load.gmtColormap('cv.cpt')
            except: cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        elif cmap=='seismic':
            import matplotlib.cm
            cmap=matplotlib.cm.seismic
        elif cmap=='bwr':
            cmap='bwr'
            cmap =discrete_cmap(int(vmax-vmin)/2, cmap)
        cmap.set_bad('0.8',0.)
        im=m.pcolormesh(x, y, self.Zarr, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        #############################################################
        # stlaLst=np.arange(5)*0.25+33.25
        # stlo=110.
        # for stla in stlaLst:
        #     xs, ys=m(stlo, stla)
        #     plt.plot(xs,ys,'^r', markersize=8)
        # stlaLst=np.arange(5)*0.25+29
        # stlo=115.
        # for stla in stlaLst:
        #     xs, ys=m(stlo, stla)
        #     plt.plot(xs,ys,'^b', markersize=15)
        #############################################################
        # try:
        #     vrange=vmin+np.arange((vmax-vmin)/0.1+1)*0.1
        #     cb = m.colorbar(im, "bottom", size="3%", pad='2%', ticks=vrange)
        # except:
        try:
            # vrange=vmin+np.arange((vmax-vmin)/0.02+1)*0.02
            # vrange=vmin+np.arange((vmax-vmin)/0.05+1)*0.05
            vrange=vmin+np.arange((vmax-vmin)/0.1+1)*0.1
            cb = m.colorbar(im, "bottom", size="3%", pad='2%', ticks=vrange)
        except:
            cb = m.colorbar(im, "bottom", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=15)
        # cb.set_label(label+' '+unit, fontsize=12, rotation=0)
        cb.set_label(label+' '+unit, fontsize=15, rotation=0)
        if contour:
            # levels=np.linspace(ma.getdata(self.Zarr).min(), ma.getdata(self.Zarr).max(), 20)
            levels=np.linspace(ma.getdata(self.Zarr).min(), ma.getdata(self.Zarr).max(), 60)
            m.contour(x, y, self.Zarr, colors='k', levels=levels, linewidths=0.5)
        if showfig: plt.show()
        return
    
    
            
                
                    
    

    

