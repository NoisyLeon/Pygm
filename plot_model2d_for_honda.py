import model2d


# 
m2d=model2d.model2d(minlon=235, maxlon=295, dlon=0.5, minlat=25, maxlat=50)
# m2d.read(infname='/projects/life9360/US_phase_map/20_ANT.vel.HD_0.05_v1')
m2d.load('US_20sec.npz')
m2d.np2ma()
m2d.plot(vmin=3.2, vmax=3.7)
# m2d.save('US_20sec.npz')


