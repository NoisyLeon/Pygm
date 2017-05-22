import model2d
import matplotlib.pyplot as plt
import numpy as np
dxArr = np.array([0.25, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
perArr= np.array([10., 20., 30., 40.])

outdir='/home/life9360/figures_grid_test'
# #
minlat=24.
maxlat=50.
minlon=-120.0+360.
maxlon=-80.+360.
for dx in dxArr:
    print 'dx =',dx
    for per in perArr:
        m2d1=model2d.model2d(infname='US_data_%g/%d_phv_eikonal_US.npz' %(dx, per) )
        m2d2=model2d.model2d(infname='US_data_0.25/%d_phv_eigen_US.npz' %(per) )
        # # # # # 
        diffm2d = m2d1-m2d2
        minlat=24.
        maxlat=50.
        minlon=-120.0
        maxlon=-80.
        diffm2d.minlat=minlat; diffm2d.maxlat=maxlat; diffm2d.minlon=minlon; diffm2d.maxlon=maxlon
        # # # 
        diffm2d.np2ma()
        import numpy as np
        mystd= np.std(diffm2d.Zarr)
        
        diffname = outdir +'/diff_eikonal_eigen_%dsec_%g.pdf' %(per, dx)
        plt.figure()
        diffm2d.plot(label='phase velocity difference',vmin=-0.1, vmax=0.1, dv=0.02, showfig=False)
        plt.savefig(diffname)
        
        histname = outdir +'/hist_eikonal_eigen_%dsec_%g.pdf' %(per, dx)
        plt.figure()
        diffm2d.hist(label='phase velocity difference ', vmin=-0.1, vmax=0.1, dv=0.02, showfig=False)
        plt.savefig(histname)