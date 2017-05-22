import model2d
import GeoPolygon
basins=GeoPolygon.GeoPolygonLst()
basins.ReadGeoPolygonLst('basin1')
minlat=23.
maxlat=51.
minlon=86.
maxlon=132.


# # 
# m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5)
# m2d.read(infname='EA_data/Tph_10.0_appV.lst')
# m2d.np2ma()
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eikonal_EA_0.0005.npz')
# 
# m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5)
# m2d.read(infname='EA_data/helm_10_phv.lst')
# m2d.np2ma()
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_helm_EA.npz')

# m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5)
# m2d.read(infname='EA_data/phV_10.lst')
# m2d.Zarr=m2d.Zarr+0.025
# # m2d.mask=m2d1.mask
# m2d.np2ma()
# m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eigen_EA.npz')

# 
# m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5)
# m2d.read(infname='EA_data/Tph_10.0_lplc.lst')
# # m2d.Zarr=m2d.Zarr+0.025
# # m2d.mask=m2d1.mask
# m2d.np2ma()
# # m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_T_lplc_EA.npz')
# 
# 
# #
# # 
m2d1=model2d.model2d(infname='./10_phv_eikonal_EA_0.0005.npz')
# m2d1=model2d.model2d(infname='./10_phv_helm_EA.npz')
m2d2=model2d.model2d(infname='./10_phv_eigen_EA.npz')
# m2d2=model2d.model2d(infname='./10_phv_helm_EA.npz')


# # 
# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Helmholtz phase velocity', vmin=2.8, vmax=3.3, geopolygons=basins)

# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Eikonal phase velocity', vmin=2.8, vmax=3.3, geopolygons=basins)
# # #
# m2d2.mask[(m2d2.lonArr>minlon)*(m2d2.lonArr<maxlon)*(m2d2.latArr>minlat)*(m2d2.latArr<maxlat) ]=0
# m2d2.np2ma()
# m2d2.minlat=minlat; m2d2.maxlat=maxlat; m2d2.minlon=minlon; m2d2.maxlon=maxlon
# m2d2.plot(label='Target phase velocity', vmin=2.8, vmax=3.4, geopolygons=basins, cmap='lasif')

diffm2d = m2d1-m2d2
diffm2d.minlat=minlat; diffm2d.maxlat=maxlat; diffm2d.minlon=minlon; diffm2d.maxlon=maxlon
# # # 
diffm2d.np2ma()
diffm2d.plot(label='phase velocity difference',vmin=-0.2, vmax=0.2, showfig=True)
diffm2d.hist(label='phase velocity difference ', vmin=-0.4, vmax=0.4)
# 
import numpy as np
print np.std(diffm2d.Zarr)
# # # m2d.smooth(10)
# # m2d.plot(vmin=2.8, vmax=3.4)
# # m2d.save('./10_phv_eikonal.npz')

# Quality control
# m2d1=model2d.model2d(infname='./10_phv_eikonal_EA.npz')
# m2dlplc=model2d.model2d(infname='./10_T_lplc_EA.npz')
# m2d2=model2d.model2d(infname='./10_phv_eigen_EA.npz')
# m2d1.mask[m2dlplc.Zarr > 0.001*1000]=1
# m2d1.mask[m2dlplc.Zarr < -0.001*1000]=1
# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Eikonal phase velocity', vmin=2.8, vmax=3.3, geopolygons=basins)
# # 
# diffm2d = m2d1-m2d2
# diffm2d.minlat=minlat; diffm2d.maxlat=maxlat; diffm2d.minlon=minlon; diffm2d.maxlon=maxlon
# 
# diffm2d.np2ma()
# diffm2d.plot(label='phase velocity difference',vmin=-0.2, vmax=0.2, showfig=True)
# diffm2d.hist(label='phase velocity difference ', vmin=-0.2, vmax=0.2)
# 
# import numpy as np
# print np.std(diffm2d.Zarr)
