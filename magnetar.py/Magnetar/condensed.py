from numpy import sin, cos, arccos, exp, log, sqrt, clip, minimum, maximum, where
from Magnetar.utils import atmosphere

class condensed_surface(atmosphere):
    def __init__(self,effective_temperature,mag_strength,mag_inclination,density,fixed_ions=True):
        self.effective_temperature=effective_temperature,mag_strength,mag_inclination
        self.mag_strength=mag_strength
        self.mag_inclination=mag_inclination
        self.blackbody_temperature=self.effective_temperature
        self.fixed_ions=fixed_ions
        self.dens=density
        self.Z=26.
        self.A=56.
    def _bbfunk(self, ee):  # per mode
        return 208452.792 * ee**3 / np.expm1(ee / self.blackbody_temperature) / 2
    #
    # This computes the emissivity from a condensed surface using the approximate treatment
    # by Potekhin et al. (2012, A&A, 546, A121) https://arxiv.org/pdf/1208.6582.pdf
    #
    def emissivity_xo(self,dataarray):
        # get dataarray values
        thetak=dataarray[-3]
        phik=dataarray[-2]
        ene=dataarray[-1]
        # magnetic field
        B13=self.mag_strength/1e13
        # useful energies
        epe=0.0288*sqrt(self.dens*self.Z/self.A)
        eci=0.0635*(self.Z/self.A)*B13
        ece=115.77*B13
        n0=sqrt(1.+epe**2/(2.*ece*eci))
        ec=eci+epe**2/ece
        ecfx=epe**2/ece
        #
        # geometric stuff
        #
        thetab=self.mag_inclination
        stb=sin(thetab)
        ctb=cos(thetab)
        stk=sin(thetak)
        ctk=cos(thetak)
        ctk=clip(ctk,0.0,0.9999)
        stk=clip(stk,(1-0.9999**2)**0.5,1.0)
        cali = (stb*stk*cos(phik)-ctb*ctk)
        calr = (stb*stk*cos(phik)+ctb*ctk)
        cali=clip(cali,-0.9999,0.9999)
        calr=clip(calr,-0.9999,0.9999)
        alphai=arccos(cali)
        alphar=arccos(calr)
        alpha=minimum(alphai,alphar)
        #
        # some energies that depend on angle
        #
        epet=epe*sqrt(3.-2.*ctk)
        ectfx=epet**2/ece
        ect=eci+epet**2/ece
        #
        # translation matrix between surface modes and X and O
        #
        ypie1i = (cos(thetab)*sin(thetak)+sin(thetab)*cos(thetak)*cos(phik))/sin(alphai)
        xpie1i = (sin(thetab)*sin(phik))/sin(alphai)
        ypie2i = (sin(thetab)*sin(phik))/sin(alphai)
        xpie2i = (-cos(thetab)*sin(thetak)-sin(thetab)*cos(thetak)*cos(phik))/sin(alphai)
        #
        ypre1r = (cos(thetab)*sin(thetak)-sin(thetab)*cos(thetak)*cos(phik))/sin(alphar)
        xpre1r = (sin(thetab)*sin(phik))/sin(alphar)
        ypre2r = (-sin(thetab)*sin(phik))/sin(alphar)
        xpre2r = (cos(thetab)*sin(thetak)-sin(thetab)*cos(thetak)*cos(phik))/sin(alphar)
        #
        if (self.fixed_ions):
            ject=1./2.+0.05*(1.+ctb*stk)/(1.+B13)-0.15*(1.-ctb)*sin(alpha)
            jeci=2.*n0*(1.+(ctb-ctk)/2./(1.+B13))/(1.+n0)**2
            pfx=0.1*(1.+stb)/(1.+B13)
            jbfx=ject/(1-pfx+pfx*(ectfx/ene)**0.6)
            ject1= 1./2. + 0.05/(1. + B13) + stb/4.
            jb1fx=ject1/(0.1+0.9*(ectfx/ene)**0.4)
            #
            # if (ene lt ectfx) then begin
            emis_case1=jbfx
            emis1_case1=jb1fx
            # endif
            # if (ene gt ectfx) then begin
            el=epe*(1.+1.2*(1.-ctk)**1.5)*(1.-stb**2/3.)
            ntfx=sqrt(1.-epet**2/(ece*(ene)))
            jcfx=4.*ntfx/(1.+ntfx)**2
            wlfx=0.8*(ectfx/epe)**0.2*sqrt(sin(alpha/2.))*(1.+stb**2)
            xfx=(ene-el)/(1.-ctk)/2./epe/wlfx
            ellfx=stk**2*wlfx*(0.17*epe/ecfx/(1.+xfx**4)+0.21*exp(-(ene/epe)**2))
            rellfx=stb**0.25*(2.-(sin(alpha))**4)*ellfx/(1.+ellfx)
            fLfx = 1./(1.+exp(5*((el-ene)/(el-ectfx))))
            emis_case2=jbfx*(1.-jcfx)+jcfx/(1.+ellfx)
            emis1_case2=jb1fx*(1.-jcfx)+jcfx*(1.-rellfx)
            # choose the right case
            emis=where(ene<ectfx,emis_case1,emis_case2)
            emis1=where(ene<ectfx,emis1_case1,emis1_case2)
        else:
            #
            aa1=1.-ctb**2*ctk-stb**2*(cos(al))
            aa1=maximum(aa1,0.99999)
            #
            # There are three cases:
            #
            #   1) energy of photon < energy of ion cyclotron
            #   2) energy of ion cyclotron < energy of photon < ect
            #   3) energy of photon > ect
            #
            # To take advantage of the numpy vectorization, we have to calculate all three possibilities
            #
            #
            #  if (ene lt eci) then begin
            n0p=sqrt(1.+epe**2/(ece*(ene+eci)))
            n0m=sqrt(1.-epe**2/(ece*(ene-eci)))
            r0p=(n0p-1.)**2/(n0p+1.)**2
            r0m=(n0m-1.)**2/(n0m+1.)**2
            j0=1.-(r0p+r0m)/2.
            j00=4./((sqrt(ec/eci)+1.)*(sqrt(eci/ec)+1.))
            a=(1.-ctb)/2./sqrt(1.+B13)+(0.7-0.45/j00)*stk**4*(1.-cos(al))
            ja=(1.-a)*j0
            ject=1./2.+0.05*(1.+ctb*stk)/(1.+B13)-0.15*(1.-ctb)*sin(al)
            jeci=2.*n0*(1.+(ctb-ctk)/2./(1.+B13))/(1.+n0)**2
            p=log(ject/jeci)/log(ect/eci)
            jb=(ene/ect)**p*ject
            ject1=0.5+0.05/(1.+B13)+0.25*stb
            jeci1=(1.-aa1)*jeci
            p1=log(ject1/jeci1)/log(ect/eci)
            jb1=(ene/ect)**p1*ject1
            a1=aa1/(1.+0.6*B13*ctb**2)
            ja1=(1.-a1)*ja
            emis_case1=maximum(ja,jb)
            emis1_case1=maximum(ja1,jb1)
            # endif   
            ject=1./2.+0.05*(1.+ctb*stk)/(1.+B13)-0.15*(1.-ctb)*sin(al)
            jeci=2.*n0*(1.+(ctb-ctk)/2./(1.+B13))/(1.+n0)**2
            p=alog(ject/jeci)/alog(ect/eci)
            jb=(ene/ect)**p*ject
            ject1=0.5+0.05/(1.+B13)+0.25*stb
            jeci1=(1.-aa1)*jeci
            p1=log(ject1/jeci1)/log(ect/eci)
            jb1=(ene/ect)**p1*ject1
            # if (ene ge eci and ene le ect) then begin
            emis_case2=jb
            emis1_case2=jb1
            # endif
            # if (ene gt ect) then begin
            nt=sqrt(1.-epet**2/(ece*(ene-eci)))
            jc=4.*nt/(1.+nt)**2
            el=epe*(1.+1.2*(1.-ctk)**1.5)*(1.-stb**2/3.)
            wl=0.8*(ect/epe)**0.2*sqrt(sin(al/2.))*(1.+stb**2)
            x=(ene-el)/(1.-ctk)/2./epe/wl
            ell=stk**2*wl*(0.17*epe/ec/(1.+x**4)+0.21*exp(-(ene/epe)**2))
            rell=stb**0.25*(2.-(sin(al))**4)*ell/(1.+ell)
            emis_case3=jb*(1.-jc)+jc/(1.+ell)
            emis1_case3=jb1*(1.-jc)+jc*(1.-rell)
            # endif
            # choose the right case
            emis=where(ene<eci,emis_case1,where(ene<ect,emis_case2,emis_case3))
            emis1=where(ene<eci,emis1_case1,where(ene<ect,emis1_case2,emis1_case3))    
        # fix the range of the values
        emis1=clip(emis1,2*emis-1,2*emis)
        emis1=clip(emis1,0,1)
        # calculate mode 2
        emis2=clip(2.*emis-emis1,0,1)
        # perform geometric transformations and clip
        emisX = clip(xpre1r**2*emis2+xpre2r**2*emis1,2*emis-1,2*emis)
        emisX = clip(emisX,0,1)
        emisO = clip(ypre1r**2*emis2+ypre2r**2*emis1,0,1)
        return emisX,emisO
    def calcIQ(self, dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return (eo+ex)*bbintens,(eo-ex)*bbintens
    def xintensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return ex*bbintens
    def ointensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return eo*bbintens
