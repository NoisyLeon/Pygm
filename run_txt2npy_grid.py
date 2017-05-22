import model2d
import numpy as np
dxArr = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
perArr= np.array([10., 20., 30., 40.])
for dx in dxArr:
    print 'dx =',dx
    for per in perArr:
        infname='US_data_%g/Tph_%.1f_appV.lst' %(dx, per)
        m2d=model2d.model2d(infname='/projects/life9360/US_phase_map/10_ANT.vel.HD_0.05_v1')
        m2d.read(infname='US_data_%g/Tph_%.1f_appV.lst' %(dx, per))
        m2d.np2ma()
        m2d.save('US_data_%g/%d_phv_eikonal_US.npz' %(dx, per))
