try:
    # can be set purposefully fail if numba code is not ready yet 
    from numba import jit
    NumbaLoaded=True
except:
    # define a dummy decorator factory so the decorators don't fail
    def jit(*args,**kwargs):
        def jitd(func):
            def helper(*args, **kwargs):
                return func(*args, **kwargs)
            return helper
        return jitd
    NumbaLoaded=False
    
from numpy import sin, cos, arccos, radians, exp, expm1, log, sqrt, minimum, maximum, where, pi, logspace, array
from Magnetar.utils import atmosphere
from scipy.integrate import simps
from Magnetar.simple_atmospheres import bbfunk

if NumbaLoaded:
    from numpy import empty_like
    @jit(nopython=True,parallel=True)
    def clip(a, a_min, a_max, out=None):
        if out is None:
            out = empty_like(a)
        for i in range(len(a)):
            if a[i] < a_min:
                out[i] = a_min
            elif a[i] > a_max:
                out[i] = a_max
            else:
                out[i] = a[i]
        return out
    @jit(nopython=True,parallel=True)
    def clipa(a, a_min, a_max, out=None):
        if out is None:
            out = empty_like(a)
        for i in range(len(a)):
            if a[i] < a_min[i]:
                out[i] = a_min[i]
            elif a[i] > a_max[i]:
                out[i] = a_max[i]
            else:
                out[i] = a[i]
        return out
    @jit(nopython=True,parallel=True)
    def choose2a(en,ethresh,choice1a,choice2a,choice1b,choice2b):
        outa=choice1a
        outb=choice1b
        for i in range(len(outa)):
            if (en[i]>ethresh[i]):
                outa[i]=choice2a[i]
                outb[i]=choice2b[i]
        return outa, outb
else:
    from numpy import clip
    def clipa(*args,**kwargs):
        return clip(*args,**kwargs)
    def choose2a(en,ethresh,choice1a,choice2a,choice1b,choice2b):
        return where(en<ethresh,choice1a,choice2a),where(en<ethresh,choice1b,choice2b)
        
class condensed_surface(atmosphere):
    def __init__(self,effective_temperature,mag_strength,mag_inclination,density,fixed_ions=True):
        self.effective_temperature=effective_temperature
        self.mag_strength=mag_strength
        self.mag_inclination=mag_inclination
        self.fixed_ions=fixed_ions
        self.dens=density
        self.Z=26.
        self.A=56.
        self.surface_temperature=self.effective_temperature
        self.adjust_surface_temperature()
    def __str__(self):
        outstring='''#
# class condensed_surface
#
# effective_temperature %12g keV
# surface_temperature   %12g keV
# mag_strength          %12g Gauss
# mag_inclination       %12g radians
# fixed_ions            %12s
# density               %12g g/cc
# Z                     %12g
# A                     %12g
#       
''' % (self.effective_temperature, self.surface_temperature, self.mag_strength, self.mag_inclination, self.fixed_ions, self.dens, self.Z, self.A )
        return outstring+atmosphere.__str__(self)
    #
    # This computes the emissivity from a condensed surface using the approximate treatment
    # by Potekhin et al. (2012, A&A, 546, A121) https://arxiv.org/pdf/1208.6582.pdf
    #
    
    @jit(nopython=True,parallel=True)
    def _emissivity_xo(dataarray,mag_strength,mag_inclination,dens,fixed_ions,Z,A):
        # get dataarray values
        thetak=radians(dataarray[-3])
        phik=radians(dataarray[-2])
        ene=dataarray[-1]
        # magnetic field
        B13=mag_strength/1e13
        # useful energies
        epe=0.0288*sqrt(dens*Z/A)
        eci=0.0635*(Z/A)*B13
        ece=115.77*B13
        n0=sqrt(1.+epe**2/(2.*ece*eci))
        ec=eci+epe**2/ece
        ecfx=epe**2/ece
        #
        # geometric stuff
        #
        thetab=radians(mag_inclination)
        stb=sin(thetab)
        ctb=cos(thetab)
        stk=sin(thetak)
        ctk=cos(thetak)
        clip(ctk,0.0,0.9999,out=ctk)
        clip(stk,(1-0.9999**2)**0.5,1.0,out=stk)
        cali = (stb*stk*cos(phik)-ctb*ctk)
        calr = (stb*stk*cos(phik)+ctb*ctk)
        clip(cali,-0.9999,0.9999,out=cali)
        clip(calr,-0.9999,0.9999,out=calr)
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
        if (fixed_ions):
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
            # emis,emis1=choose2a(ene,ectfx,emis_case1,emis_case2,emis1_case1,emis1_case2)
        else:
            #
            aa1=1.-ctb**2*ctk-stb**2*(cos(alpha))
            aa1=minimum(aa1,0.99999)
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
            a=(1.-ctb)/2./sqrt(1.+B13)+(0.7-0.45/j00)*stk**4*(1.-cos(alpha))
            ja=(1.-a)*j0
            ject=1./2.+0.05*(1.+ctb*stk)/(1.+B13)-0.15*(1.-ctb)*sin(alpha)
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
            ject=1./2.+0.05*(1.+ctb*stk)/(1.+B13)-0.15*(1.-ctb)*sin(alpha)
            jeci=2.*n0*(1.+(ctb-ctk)/2./(1.+B13))/(1.+n0)**2
            p=log(ject/jeci)/log(ect/eci)
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
            wl=0.8*(ect/epe)**0.2*sqrt(sin(alpha/2.))*(1.+stb**2)
            x=(ene-el)/(1.-ctk)/2./epe/wl
            ell=stk**2*wl*(0.17*epe/ec/(1.+x**4)+0.21*exp(-(ene/epe)**2))
            rell=stb**0.25*(2.-(sin(alpha))**4)*ell/(1.+ell)
            emis_case3=jb*(1.-jc)+jc/(1.+ell)
            emis1_case3=jb1*(1.-jc)+jc*(1.-rell)
            # endif
            # choose the right case
        # fix the range of the values
        #emis=emis_case1
        #emis1=emis1_case1
        #for i in range(len(emis)):
        #    if ene[i]>ectfx[i]:
        #        emis[i]=emis_case2[i]
        #        emis1[i]=emis1_case2[i]
        if fixed_ions:        
            ethresh1=ectfx
        else:
            ethresh1=eci+0*ene
        emis,emis1=choose2a(ene,ethresh1,emis_case1,emis_case2,emis1_case1,emis1_case2)
        if not fixed_ions:        
            emis,emis1=choose2a(ene,ect,emis,emis_case3,emis1,emis1_case3)
        clipa(emis1,2*emis-1,2*emis,out=emis1)
        clip(emis1,0.0,1.0,out=emis1)
        # calculate mode 2
        emis2=2.0*emis-emis1
        clip(emis2,0.0,1.0,out=emis2)
        # perform geometric transformations and clip
        emisX=xpre1r**2*emis2+xpre2r**2*emis1
        #clipa(emisX,2*emis-1,2*emis,out=emisX)
        clip(emisX,0.0,1.0,out=emisX)
        emisO=ypre1r**2*emis2+ypre2r**2*emis1
        clip(emisO,0.0,1.0,out=emisO)
        return emisX,emisO
    def emissivity_xo(self,dataarray):
        #a,b = condensed_surface._emissivity_xo(array(dataarray),self.mag_strength,self.mag_inclination,self.dens,self.fixed_ions,self.Z,self.A)
        #print(a,b)
        #return a,b
        return condensed_surface._emissivity_xo(array(dataarray),self.mag_strength,self.mag_inclination,self.dens,self.fixed_ions,self.Z,self.A)
    def calcIQ(self, dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return (eo+ex)*bbintens,(eo-ex)*bbintens
    def xintensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return ex*bbintens
    def ointensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return eo*bbintens
