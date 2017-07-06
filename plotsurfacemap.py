import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

import sys

size=0.1
skip=4

for n,f in enumerate(sys.argv[1:]):
   fig = plt.figure(n)
   ax = fig.add_subplot(111, aspect='equal')
   b,beta,s1,s2,s3,mago,o1,o2,o3,mag_colat,theta,phi,X,O,Q,IOmdL =  np.loadtxt(f,unpack=True)
   bmax=max(b)
   nb=(int) (len(b)/2)**0.5
   if nb>10:
      skip=(int) (nb/10)
      nbf=nb/skip
      nstart=np.arange(nbf)
      istart=nstart*nstart*skip*skip*2
      ival=[]
      for ii,jj in zip(istart[:-1],istart[1:]):
         ival=np.concatenate((ival,np.arange(ii,jj,skip)))
      ival=ival.astype(int)
      b=b[ival]
      beta=beta[ival]
      s1=s1[ival]
      s2=s2[ival]
      s3=s3[ival]
   for bb,be,s1a,s2a,s3a in zip(b/bmax,beta,s1,s2,np.abs(s3)):
       e = Ellipse(xy=[bb*np.cos(be),bb*np.sin(be)],width=size,height=size*s3a/(s1a*s1a+s2a*s2a+s3a*s3a)**0.5,angle=(np.arctan2(s2a,s1a)*0.5+be)*180/3.1415)
       ax.add_artist(e)
       e.set_clip_box(ax.bbox)
       e.set_alpha(0.6)
       e.set_facecolor([0.9,0.0,0.0])
   for bb,be,s1a,s2a,s3a in zip(b/bmax,-beta,s1,-s2,np.abs(s3)):
       e = Ellipse(xy=[bb*np.cos(be),bb*np.sin(be)],width=size,height=size*s3a/(s1a*s1a+s2a*s2a+s3a*s3a)**0.5,angle=(np.arctan2(s2a,s1a)*0.5+be)*180/3.1415)
       ax.add_artist(e)
       e.set_clip_box(ax.bbox)
       e.set_alpha(0.6)
       e.set_facecolor([0.9,0.0,0.0])

   plt.xlim(-1.1,1.1)
   plt.ylim(-1.1,1.1)
   plt.show()
    
