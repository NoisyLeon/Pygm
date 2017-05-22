import model2d

per=40.
# vmin=2.8; vmax=3.4
# vmin=3.3; vmax=3.7
# vmin=3.5; vmax=4.0
vmin=3.7; vmax=4.1


minlat=24.
maxlat=50.
minlon=-120.0+360.
maxlon=-80.+360.

minlat=24.+3
maxlat=50.-3.
minlon=-120.0+360.+3.
maxlon=-80.+360.-3.

# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.np2ma()
# # m2d.plot(vmin=2.8, vmax=3.4, showfig=True)
# m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_smooth.npz')
# 
m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat)
m2d.read(infname='US_data_0.25/Tph_%.1f_appV.lst' %per)
m2d.np2ma()
# m2d.plot(vmin=vmin, vmax=vmax)
m2d.save('US_data_0.25/%d_phv_eikonal_US.npz' %per)
# # 
# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='US_data/helm_%d_phv.lst' %per)
# m2d.np2ma()
# m2d.plot(vmin=vmin, vmax=vmax)
# m2d.save('US_data/%d_phv_helm_US.npz' %per)
# # 
m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat)
m2d.read(infname='US_data_0.25/phV_%d.lst' %per)
m2d.Zarr=m2d.Zarr+0.05
# m2d.mask=m2d1.mask
m2d.np2ma()
m2d.smooth(10)
m2d.plot(vmin=vmin, vmax=vmax)
m2d.save('US_data_0.25/%d_phv_eigen_US.npz' %per)


# #
# minlat=24.
# maxlat=50.
# minlon=-120.0+360.
# maxlon=-80.+360.
m2d1=model2d.model2d(infname='US_data_0.25/%d_phv_eikonal_US.npz' %per)
# m2d1=model2d.model2d(infname='US_data/%d_phv_helm_US.npz' %per)
# m2d2=model2d.model2d(infname='US_data/%d_phv_helm_US.npz' %per)
m2d2=model2d.model2d(infname='US_data_0.25/%d_phv_eigen_US.npz' %per)
# # m2d2.mask[(m2d2.lonArr>minlon)*(m2d2.lonArr<maxlon)*(m2d2.latArr>minlat)*(m2d2.latArr<maxlat) ]=0
# 
# # m2d2=model2d.model2d(infname='./10_phv_helm_US.npz')
# 
# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Eikonal phase velocity', vmin=vmin, vmax=vmax, dv=0.1)
# # 
# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Helmholtz phase velocity', vmin=vmin, vmax=vmax)
# # # # #
# m2d2.mask=m2d1.mask
# m2d2.mask[(m2d2.lonArr>minlon)*(m2d2.lonArr<maxlon)*(m2d2.latArr>minlat)*(m2d2.latArr<maxlat) ]=0
# m2d2.np2ma()
# m2d2.minlat=minlat; m2d2.maxlat=maxlat; m2d2.minlon=minlon; m2d2.maxlon=maxlon
# m2d2.plot(label='Target phase velocity', vmin=vmin, vmax=vmax, dv=0.1)
# # # # # 
diffm2d = m2d1-m2d2
diffm2d.minlat=minlat; diffm2d.maxlat=maxlat; diffm2d.minlon=minlon; diffm2d.maxlon=maxlon
# # # # 
diffm2d.np2ma()
# diffm2d.plot(label='phase velocity difference',vmin=-0.02, vmax=0.02, dv=0.01, showfig=True)
# diffm2d.hist(label='phase velocity difference ', vmin=-0.02, vmax=0.02, dv=0.01)

diffm2d.plot(label='phase velocity difference',vmin=-0.1, vmax=0.1, dv=0.02, showfig=True)
diffm2d.hist(label='phase velocity difference ', vmin=-0.1, vmax=0.1, dv=0.02)
# 
# diffm2d.plot(label='phase velocity difference',vmin=-0.2, vmax=0.2, dv=0.05, showfig=True)
# diffm2d.hist(label='phase velocity difference ', vmin=-0.2, vmax=0.2, dv=0.05)

import numpy as np
print np.std(diffm2d.Zarr)

# # m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eikonal.npz')