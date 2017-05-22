import model2d

# minlat=23.
# maxlat=51.
# minlon=86.
# maxlon=132.

minlat=24.
maxlat=50.
minlon=-120.0+360.
maxlon=-80.+360.

minlat=25
maxlat=50.
minlon=235
maxlon=295
# # 
m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.05, minlat=minlat, maxlat=maxlat, dlat=0.05)
# m2d.read(infname='/projects/life9360/code/SES3DPy/diffA.lst')
# m2d.read(infname='US_data/phV_100.lst')
#m2d.read(infname='/projects/life9360/code/SES3DPy/MAPS_US/smpkolya_phv_R_40')
# m2d.read(infname='/projects/life9360/US_phase_map/8_ANT.vel.HD_0.05_v1')
m2d.read_numpy_txt(infname='/projects/life9360/US_phase_map/8_ANT.vel.HD_0.05_v1')
m2d.np2ma()
m2d.plot(vmin=2.8, vmax=3.4)
# m2d.smooth(10)
#m2d.plot(vmin=3.6, vmax=4.1)
# m2d.plot(vmin=-20, vmax=20, cmap='bwr')
# m2d.save('./10_phv_beam.npz')
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
# m2d1=model2d.model2d(infname='./10_phv_eikonal.npz')
# m2d2=model2d.model2d(infname='./10_phv_eigen.npz')
# diffm2d = m2d1-m2d2
# 
# 
# diffm2d.np2ma()
# diffm2d.plot(vmin=-0.2, vmax=0.2, showfig=True)
# diffm2d.hist()
# # m2d.smooth(10)
# m2d.plot(vmin=2.8, vmax=3.4)
# m2d.save('./10_phv_eikonal.npz')
