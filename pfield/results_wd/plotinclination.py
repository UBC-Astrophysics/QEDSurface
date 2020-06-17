import numpy as np
import matplotlib.pyplot as plt

angle, fraction = np.loadtxt("2008MNRAS_386_2167M_3T.dat",unpack=True)

mu, nu, m, r, inclination, sumq, nq = np.loadtxt("cone_summary_6.1",unpack=True)

for rr in [9.2e8, 5.5e8]:
    for nn in [0,0.5e18,2e18]:
        ii=np.logical_and(nu==nn,np.logical_and(mu==1e35,r==rr))
        x1=inclination[ii]
        y1=-np.interp(inclination[ii],angle,fraction)*sumq[ii]/nq[ii]
        y1=y1[np.argsort(x1)]
        x1=x1[np.argsort(x1)]
        if (nn==0):
            plt.plot(x1,y1,label='QED-off')
        else:
            plt.plot(x1,y1,label='E = %3.1f keV' % (nn*4.14e-18) )
    plt.legend(bbox_to_anchor=(0.35, 1))
    plt.title(r"Radius = %4.0f km, $\mu = 10^{35}$ G cm$^{3}$" % (rr/1e5))
    plt.xlabel("Inclination [degrees]")
    plt.ylabel("Polarization Fraction")
    plt.savefig("wd_61_%5.1e.png" % rr)
    plt.close()
