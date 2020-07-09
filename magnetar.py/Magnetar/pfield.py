import numpy as np
class pfield:
    def __init__(self):
        self.data = []
        self.imean = 0
        self.qmean = 0
        self.mu = 0
        self.nu = 0
        self.mass = 0
        self.radius = 0
        self.theta = 0
        self.angkb = []
        self.ndotb = []
        self.pt = 1
        self.filename = ""
        self.ebins=[]
        self.iint=[]
        self.qint=[]        
        self.oneSided=True

    def loaddata(self, pfield_file):
        with open(pfield_file, 'r') as f:
            line = f.readline()
            a = line.split()
            try:
                self.mu = float(a[2])
                self.nu = float(a[3])
                self.mass = float(a[4])
                self.radius = float(a[5])
                self.theta = float(a[6])
            except:
                pass
        self.data = np.genfromtxt(
            pfield_file, unpack=True, names=True, skip_header=25)
        self.data = self.data[~np.isnan(self.data['s1'])]
        self.imean = np.mean(self.data['X'] + self.data['O'])
        self.qmean = np.mean(self.data['Q'] *
                             (self.data['X'] - self.data['O']))
        hld = np.cos(np.radians(self.data['mag_colat']))**2
        self.ndotb = (4 * hld / (3 * hld + 1))**0.5
        self.angkb = np.degrees(
            np.arccos((self.ndotb * np.cos(np.radians(self.data['theta'])) +
                       (1 - self.ndotb**2)**0.5 * np.sin(
                           np.radians(self.data['theta'])) * np.cos(
                               np.radians(self.data['phi'])))))
        self.filename=pfield_file
        return self

    def calcIQ(self, energies, atmo_map):
        ones = np.ones(len(energies))
        dataarray = np.array([
            np.kron(ones, self.data['theta']),
            np.kron(ones, self.data['phi']),
            np.kron(energies, np.ones(len(self.data['phi'])))
        ])
        ii, qq = atmo_map.calcIQ(dataarray)
        ii = np.reshape(ii, (len(energies), len(self.data['phi'])))
        qq = np.reshape(qq, (len(energies), len(self.data['phi'])))
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        return np.mean(ii, axis=1), np.mean(-qq * dotprod, axis=1)

    def calcmeanIQ(self, energies, atmo_map):
        ones = np.ones(len(energies))
        dataarray = np.array([
            np.kron(ones, self.data['mag_colat']),
            np.kron(ones, self.angkb),
            np.kron(energies, np.ones(len(self.angkb)))
        ])
        ii, qq = atmo_map.calcmeanIQ(dataarray)
        ii = np.reshape(ii, (len(energies), len(self.angkb)))
        qq = np.reshape(qq, (len(energies), len(self.angkb)))
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        return np.mean(ii, axis=1), -np.mean(qq * dotprod, axis=1)

    def recalculate(self, energy, atmo_map, gtt=1):
        dataarray = np.array([
            self.data['mag_colat'], self.data['theta'], self.data['phi'],
            self.pt * np.full(len(self.data['mag_colat']), energy / gtt)
        ])
        ii, qq = atmo_map.calcIQ(dataarray)
        self.data['X'] = 0.5 * (ii - qq)
        self.data['O'] = 0.5 * (ii + qq)
        self.imean = np.mean(ii) * gtt**3
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        self.qmean = np.mean(-dotprod * qq) * gtt**3
        #self.qmean=np.mean(-self.data['Q']*qq)*gtt**3
        return self.imean, self.qmean

    def recalculatephiaverage(self, energy, atmo_map, gtt=1):
        dataarray = [
            self.data['mag_colat'], self.angkb,
            self.data['mag_colat'] * 0 + energy / gtt
        ]
        ii, qq = atmo_map.calcmeanIQ(dataarray)
        self.data['X'] = (ii - qq) * 0.5
        self.data['O'] = (ii + qq) * 0.5
        self.imean = np.mean(ii) * gtt**3
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        self.qmean = np.mean(-dotprod * qq) * gtt**3
        return self.imean, self.qmean

    def calcvalues(self,surfacemodel,ebins=None,gtt=1):
        if ebins is None:
            ebins=np.logspace(-0.5,1.5,100)
        self.ebins=np.array(ebins)
        self.iint=[]
        self.qint=[]
        for en in ebins:
            self.recalculate(en,surfacemodel,gtt=gtt)
            self.iint.append(self.imean)
            self.qint.append(self.qmean)
        self.iint=np.array(self.iint)
        self.qint=np.array(self.qint)
    def __str__(self):
        outstring='#   Phi[rad]  Energy[keV]            I          Q/I\n'
        for ee,ii,rr in zip(self.ebins,self.iint,self.qint/self.iint):
            outstring='%s%12g %12g %12g %12g\n' % (outstring,self.theta*np.pi/180.0,ee,ii,rr)
        return outstring
        
    def plot(self,
             maxrings=15,
             size=0.08,
             plotpoints=False,
             drawedge=True,
             drawaxes=False,
             contours=[15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180],
             lev=100,
             datamap=None,
             cmap="viridis",
             ellipsecolor=[0.9, 0.0, 0.0],
             contourcolor='k',
             oneSided=True):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse

        def _doellipse(bb, be, s1a, s2a, s3a):
            e = Ellipse(
                xy=[bb * np.cos(be), bb * np.sin(be)],
                width=size,
                height=size * s3a / (s1a * s1a + s2a * s2a + s3a * s3a)**0.5,
                angle=(np.arctan2(s2a, s1a) * 0.5 + be) * 180 / 3.1415)
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.6)
            e.set_facecolor(ellipsecolor)
            e.set_edgecolor(ellipsecolor)
            e.set_linewidth(0.5)
            return e

        fig = plt.figure(0, figsize=(8, 8))
        ax = fig.add_subplot(111, aspect='equal')
        b, beta, s1, s2, s3, mag_colat = self.data['b'], self.data[
            'beta'], self.data['s1'], self.data['s2'], self.data[
                's3'], self.data['mag_colat']
        bmax = max(b)
        x = b * np.cos(beta) / bmax
        y = b * np.sin(beta) / bmax
        if contours is not None:
            ax.tricontour(
                np.concatenate((x, x)),
                np.concatenate((y, -y)),
                np.concatenate((mag_colat, mag_colat)),
                levels=contours,
                linewidths=0.5,
                colors=contourcolor)
        if datamap is None:
            if drawedge:
                e = Ellipse(xy=[0, 0], width=2, height=2)
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_alpha(1.0)
                e.set_facecolor([1.0, 1.0, 1.0])
                e.set_edgecolor([0.0, 0.0, 0.0])
                e.set_linewidth(1.0)
        else:
            if oneSided:
                cntr2 = ax.tricontourf(
                    np.concatenate((x, x)),
                    np.concatenate((y, -y)),
                    np.concatenate((datamap, datamap)),
                    levels=lev,
                    cmap=cmap)
            else:
                cntr2 = ax.tricontourf(x, y, datamap, levels=lev, cmap=cmap)
            fig.colorbar(cntr2, ax=ax)

        bmax = max(b)
        nb = (int)(len(b) / 2)**0.5
        if nb > maxrings:
            skip = (int)(nb / maxrings)
            nbf = nb / skip
            nstart = np.arange(0, nbf + 1)
            istart = nstart * nstart * skip * skip * 2
            iend = (nstart * skip + 1) * (nstart * skip + 1) * 2
            ival = []
            for ii, jj in zip(istart, iend):
                ival = np.concatenate((ival, np.arange(ii + skip / 2, jj,
                                                       skip)))
            ival = ival[ival < len(b)].astype(int)
            b, beta, s1, s2, s3 = b[ival], beta[ival], s1[ival], s2[ival], s3[
                ival]

        if plotpoints:
            ax.plot(b / bmax * np.cos(beta), b / bmax * np.sin(beta), 'b.')
            if oneSided:
                ax.plot(b / bmax * np.cos(beta), -b / bmax * np.sin(beta),
                        'b.')
        for bb, be, s1a, s2a, s3a in zip(b / bmax, beta, s1, s2, np.abs(s3)):
            _doellipse(bb, be, s1a, s2a, s3a)
            if oneSided:
                _doellipse(bb, -be, s1a, -s2a, s3a)

        if drawaxes == False:
            plt.axis('off')

        return fig
    
        def display(self):
            return self.plot(datamap=np.log10(self.data['X']+self.data['O']),cmap='inferno',ellipsecolor=[0,0.9,0])

class pfield_array:
    def __init__(self):
        self.mufil,self.pfi = [],[]
    def loaddata(self, files, modeltype=None):
        self.mufil,self.pfi = [],[]
        if type(files) is str:
            files = [
                files,
            ]
        for i in files:
            pfii= pfield()
            pfii.loaddata(i)
            self.pfi.append(pfii)
            self.mufil.append(pfii.theta)
        ind = np.argsort(self.mufil)
        self.mufil = [self.mufil[i] for i in ind]
        self.pfi = [self.pfi[i] for i in ind]
        return self
    def calcvalues(self,surfacemodel,ebins=None,gtt=1):
        for pf in self.pfi:
            pf.calcvalues(surfacemodel,ebins,gtt=1)
    def __str__(self):
        outstring=''
        for pf in self.pfi:
            outstring=outstring+str(pf)
        return outstring
    def recalculate(self,energy,surfacemodel,gtt=1):
        for pf in self.pfi:
            pf.recalculate(energy,surfacemodel,gtt=gtt)
    def loadfiles(self,dirpath):
        import os
        import fnmatch
        from tkinter import Tcl
        import glob
        import matplotlib.pylab as plt

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)
        #define energy array and asociated NH absorption
        self.e_pfield   = np.logspace(-0.5,2.0,30)
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(self.e_pfield,eabs,sabs)/(self.e_pfield)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.6e22*ssabs)
        # load surface emission 2d-array (energy and theta)
        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,surfaceaux,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
            self.i_pfield2d.append(ivec)
            self.iabs_pfield2d.append(totabs*ivec/1.7e4)
            self.q_pfield2d.append(qvec)
            #plt.loglog(self.e_pfield,ivec/self.e_pfield**2)
            if(pfieldaux.theta<=180.0):
                ax1.plot(self.e_pfield,np.array(totabs*ivec/1.7e4))
                ax2.plot(self.e_pfield,np.array(qvec)/np.array(ivec))
        #plt.ylim(1,1e4)
        ax1.set_yscale("log")
        plt.show()
        self.i_pfield2d = np.array(self.i_pfield2d)
        self.iabs_pfield2d = np.array(self.iabs_pfield2d)
        self.q_pfield2d = np.array(self.q_pfield2d)
        return self.e_pfield

    def loadfullcompton(self,dirpath):
        import os
        import fnmatch
        from tkinter import Tcl
        import glob
        import matplotlib.pylab as plt

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)
        #define energy array and asociated NH absorption
        self.e_pfield   = np.logspace(-0.5,2.0,30)
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(self.e_pfield,eabs,sabs)/(self.e_pfield)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.6e22*ssabs)
        # load surface emission 2d-array (energy and theta)
        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #repeat
            allsurface150=surface_model().loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=allsurface150.mcolat
            # add the angles
            allsurface150.mcolat=aa+[180-i for i in aa[::-1]]
            # add the patches
            allsurface150.patches=allsurface150.patches+allsurface150.patches[::-1]
            #apply complete comptonization on every slab
            ii150,qq150 = [],[]
            for ii, pp in enumerate(surfaceaux.patches):
                allsurface150.patches[ii] = Complete_Comptonization(pp,2.1)
                ii150int,qq150int=allsurface150.patches[ii].fluxIQ(self.e_pfield)
                ii150.append(ii150int)
                qq150.append(qq150int)

            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,allsurface150,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
            self.i_pfield2d.append(ivec)
            self.iabs_pfield2d.append(totabs*ivec/1.7e4)
            self.q_pfield2d.append(qvec)
            #plt.loglog(self.e_pfield,ivec/self.e_pfield**2)
            if(pfieldaux.theta<=180.0):
                ax1.plot(self.e_pfield,np.array(totabs*ivec/1.7e4))
                ax2.plot(self.e_pfield,np.array(qvec)/np.array(ivec))
        #plt.ylim(1,1e4)
        ax1.set_yscale("log")
        plt.show()
        self.i_pfield2d = np.array(self.i_pfield2d)
        self.iabs_pfield2d = np.array(self.iabs_pfield2d)
        self.q_pfield2d = np.array(self.q_pfield2d)
        return self.e_pfield
    
    def loadRCS(self,dirpath):
        import numpy as np
        import matplotlib.pyplot as plt
        import glob


        # Load a surface model
        ee=np.logspace(-2,2,401)
        allsurface=surface_model().loaddata(glob.glob('../QEDSurface/B14.11T6.48/*.int'))

        # Calculate the spectra
        iilist,qqlist, mcolat = [],[],[]
        for ii, pp in enumerate(allsurface.patches):
            iiint,qqint = pp.fluxIQ(ee)
            iilist.append(iiint)
            qqlist.append(qqint)

        eo = ee/1.3
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(eo,eabs,sabs)/(eo)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.3e22*ssabs)

        #from Magnetar import partial_res_comptonization
        new_patch150=partial_res_comptonization(allsurface.patches[3],-0.151,150)
        #comp_patch3=allsurface.patches[3]**0.151
        ii150,qq150=new_patch150.fluxIQ(ee)
        new_patch150u=partial_res_comptonization(allsurface.patches[3],0.151,150)
        #comp_patch3=allsurface.patches[3]**0.151
        ii150u,qq150u=new_patch150u.fluxIQ(ee)
        factor=4.5e-6 # 0.346 keV

        plt.semilogx(eo,qq150u/ii150u,label='H  Compt 150 u')
        plt.show()
        plt.loglog(eo,factor*totabs*(ii150u*ee),'k-.',label='H atmo + RCS')
        plt.legend()
        plt.xlabel("Energy [keV]")
        plt.ylabel(r'$E f_E$ [keV$^{-2}$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
        #plt.loglog(ee,nout*ee**3)
        plt.xlim(0.5,1e2/1.3)
        plt.ylim(1e-4,1e-1)
        plt.show()
        
        self.e_pfield = ee
        self.th_pfield = np.linspace(0,360,4)
        self.i_pfield2d = np.concatenate([[ii150u]] * 4, axis=0)
        self.iabs_pfield2d = np.concatenate([[factor*totabs*(ii150u*ee)+1e-20]]* 4, axis=0)
        self.q_pfield2d = np.concatenate([[qq150u]] * 4, axis=0)

    def loadRCS2(self,filepath):
        import os
        import fnmatch
        from tkinter import Tcl

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)

        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/B14.11T6.48/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,surfaceaux,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
      

    def lightcurve(self, inclination, field_angle, energy_lc):
        from scipy.interpolate import RectBivariateSpline
        import matplotlib.pylab as plt

        self.energy_lc = energy_lc
        self.qoveri =  self.q_pfield2d/self.i_pfield2d

        i_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, self.iabs_pfield2d)
        log_i_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, np.log10(self.iabs_pfield2d))
        q_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, self.q_pfield2d)
        qoveri_spline   = RectBivariateSpline(self.th_pfield, self.e_pfield, self.qoveri)

        self.inclination = np.radians(inclination)
        self.field_angle = np.radians(field_angle)
        self.phase_array = np.radians(np.linspace(0.36, 360.05730682, 10))

        self.theta_array = np.arccos(
            np.cos(self.field_angle) * np.cos(self.inclination) +
            np.sin(self.field_angle) * np.sin(self.inclination
                                              ) * np.cos(self.phase_array))

        self.eta_array = np.arcsin(np.sin(self.phase_array)*np.sin(self.field_angle)/
                                   np.sin(self.theta_array))

        outputfile = open('ixpesw-ixpeobssim-572fa0092131/ixpeobssim/config/ascii/mag_4U_2hotspots.dat', 'w')
        outputfile.write("\t%.4f\t%.4f\n" % (inclination, field_angle))
        #outputfile.write("%.4f   %.4f\n" % (inclination, field_angle))

        #fig, ax1 = plt.subplots()
        #fig, ax2 = plt.subplots()
        i_array = []
        q_array = []
        i_array2 = []
        qoveri_array =[]
        pa_array = []
        for i in range(len(self.phase_array)-1):
            i_aux   = 10.0**log_i_interp_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            q_aux   =           q_interp_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            qoveri_aux  =         qoveri_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            pa_aux = 0.5*np.arctan2(-np.sin(2.0*self.eta_array[i])*qoveri_aux[0],np.cos(2.0*self.eta_array[i])*qoveri_aux[0])
            pa_aux = np.where(pa_aux<0,np.degrees(np.pi+pa_aux),np.degrees(pa_aux))
            i_array.append(i_aux[0])
            q_array.append(q_aux[0])
            qoveri_array.append(qoveri_aux[0])
            pa_array.append(pa_aux)
            #ax1.plot(self.energy_lc,i_aux[0])
            #ax2.plot(self.energy_lc,qoveri_aux[0])
            #save model to file
            outputfile.write("\t%d\t%.8f\t%.8f\n" % (i+1,self.phase_array[i], self.phase_array[i+1]))
            #outputfile.write("%d %.8f %.8f\n" % (i+1,self.phase_array[i], self.phase_array[i+1]))

            for j in range(len(self.energy_lc)-1):
                e0 = np.log10(self.energy_lc[j])
                e1 = np.log10(self.energy_lc[j+1])
                i0 = i_aux[0][j]/(self.energy_lc[j])**2.0
                pd = np.abs(qoveri_aux[0][j])
                pd = np.where(pd>1.0, 0.99, pd)
                pa = pa_aux[j]
                #print("\t%.6f\t%.6f\t%.7f\t%.6f\t%.4f\n" % (e0,e1,i0,pd,pa))
                outputfile.write("\t%.6f\t%.6f\t%.7e\t%.6f\t%.4f\n" % (e0,e1,i0,pd,pa))

        outputfile.close()
        #ax1.set_yscale("log")
        self.i_array = np.array(i_array)
        self.q_array = np.array(q_array)
        self.qoveri_array = np.array(qoveri_array)
        self.pa_array = np.array(pa_array)

        return 10
