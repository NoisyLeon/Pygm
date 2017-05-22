import model2d
import matplotlib.pyplot as plt
per=20.
# vmin=3.10; vmax=3.40
vmin=3.4; vmax=3.7
# vmin=3.65; vmax=4.0
# vmin=3.75; vmax=4.15
# minlat=24.
# maxlat=50.
# minlon=-120.0+360.
# maxlon=-80.+360.


minlat=24.+3
maxlat=50.-3.
minlon=-120.0+360.+3.
maxlon=-80.+360.-3.

m2d=model2d.model2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat)
# m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
m2d.read(infname='US_data_ak135/Tph_%.1f_appV.lst' %per)
m2d.save('US_data_ak135/%d_phv_eikonal_US.npz' %per)
m2d=model2d.model2d(infname='US_data_ak135/%d_phv_eikonal_US.npz' %per)
m2d.minlat=minlat; m2d.maxlat=maxlat; m2d.minlon=minlon; m2d.maxlon=maxlon
m2d.np2ma()
m2d.plot(label='Eikonal phase velocity', cmap='seismic_r', vmin=vmin, vmax=vmax, dv=0.05)
m2d.Zarr = m2d.Zarr - m2d.Zarr.mean()
# m2d.hist(label='phase velocity perturbation ', vmin=-0.05, vmax=0.05, dv=0.01, showfig=False)
m2d.hist(label='phase velocity perturbation ', vmin=-0.05, vmax=0.05, dv=0.01, showfig=False)
plt.show()