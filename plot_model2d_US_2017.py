import model2d

# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.np2ma()
# m2d.plot(vmin=2.8, vmax=3.4, showfig=False)
# m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_smooth.npz')
# # 
# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='Tph_10.0_appV.lst')
# m2d.np2ma()
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eikonal.npz')

# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='helm_10_phv.lst')
# m2d.np2ma()
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_helm.npz')
# m2d1=model2d.model2d(infname='./10_phv_eikonal.npz')
# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
# m2d.read(infname='phV_10.lst')
# m2d.Zarr=m2d.Zarr+0.025
# m2d.mask=m2d1.mask
# m2d.np2ma()
# m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eigen.npz')


#
minlat=24.
maxlat=50.
minlon=-120.0+360.
maxlon=-80.+360.
m2d1=model2d.model2d(infname='./10_phv_eikonal_US.npz')
# m2d1=model2d.model2d(infname='./10_phv_helm_US.npz')
m2d2=model2d.model2d(infname='./10_phv_helm_US.npz')
# m2d2=model2d.model2d(infname='./10_phv_eigen_US.npz')
# m2d2.mask[(m2d2.lonArr>minlon)*(m2d2.lonArr<maxlon)*(m2d2.latArr>minlat)*(m2d2.latArr<maxlat) ]=0

# m2d2=model2d.model2d(infname='./10_phv_helm_US.npz')

# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Eikonal phase velocity', vmin=2.8, vmax=3.4)

# m2d1.np2ma()
# m2d1.minlat=minlat; m2d1.maxlat=maxlat; m2d1.minlon=minlon; m2d1.maxlon=maxlon
# m2d1.plot(label='Helmholtz phase velocity', vmin=2.8, vmax=3.4)
# # #
# m2d2.mask[(m2d2.lonArr>minlon)*(m2d2.lonArr<maxlon)*(m2d2.latArr>minlat)*(m2d2.latArr<maxlat) ]=0
# m2d2.np2ma()
# m2d2.minlat=minlat; m2d2.maxlat=maxlat; m2d2.minlon=minlon; m2d2.maxlon=maxlon
# m2d2.plot(label='Target phase velocity', vmin=2.8, vmax=3.4)
# # # 
diffm2d = m2d1-m2d2
minlat=24.
maxlat=50.
minlon=-120.0
maxlon=-80.
diffm2d.minlat=minlat; diffm2d.maxlat=maxlat; diffm2d.minlon=minlon; diffm2d.maxlon=maxlon
# # 
diffm2d.np2ma()
diffm2d.plot(label='phase velocity difference',vmin=-0.2, vmax=0.2, showfig=True)
diffm2d.hist(label='phase velocity difference ', vmin=-0.1, vmax=0.1)


import numpy as np
print np.std(diffm2d.Zarr)

# # # m2d.smooth(10)
# # m2d.plot(vmin=2.8, vmax=3.4)
# # m2d.save('./10_phv_eikonal.npz')