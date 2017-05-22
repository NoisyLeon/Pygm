import model2d
import matplotlib.pyplot as plt
import numpy as np
dxArr = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
perArr= np.array([10., 20., 30., 40.])

# vmin=2.8; vmax=3.4
# vmin=3.3; vmax=3.7
# vmin=3.5; vmax=4.0
vmin=3.7; vmax=4.1
vArr = [[2.8, 3.4], [3.3, 3.7], [3.5, 4.0], [3.7, 4.1]]

outdir='/home/life9360/figures_grid_test'
# #
minlat=24.
maxlat=50.
minlon=-120.0+360.
maxlon=-80.+360.
for dx in dxArr:
    print 'dx =',dx
    for per, v in zip(perArr, vArr):
        vmin=v[0]; vmax= v[1]
        m2d=model2d.model2d(infname='US_data_%g/%d_phv_eikonal_US.npz' %(dx, per) )
        # # # # #
        plt.figure()
        m2d.np2ma()
        m2d.minlat=minlat; m2d.maxlat=maxlat; m2d.minlon=minlon; m2d.maxlon=maxlon
        m2d.plot(label='Eikonal phase velocity', vmin=vmin, vmax=vmax, dv=0.1, showfig=False)
        phVname = outdir +'/eikonalV_%dsec_%g.pdf' %(per, dx)
        plt.savefig(phVname)