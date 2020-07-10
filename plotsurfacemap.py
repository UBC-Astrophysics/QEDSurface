import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

import sys

def _doellipse(bb,be,s1a,s2a,s3a,size,ax):
   e = Ellipse(xy=[bb*np.cos(be),bb*np.sin(be)],width=size,height=size*s3a/(s1a*s1a+s2a*s2a+s3a*s3a)**0.5,angle=(np.arctan2(s2a,s1a)*0.5+be)*180/3.1415)
   ax.add_artist(e)
   e.set_clip_box(ax.bbox)
   e.set_alpha(0.6)
   e.set_facecolor([0.9,0.0,0.0])
   e.set_edgecolor([0.9,0.0,0.0])
   e.set_linewidth(0.5)

def plotfiles(filelist,
              maxrings=15,size=0.08,plotpoints=False,drawedge=True,drawaxes=False,
              contours=[15,30,45,60,75,90,105,120,135,150,165,180]):

   for n,f in enumerate(filelist):
      fig = plt.figure(n,figsize=(8,8))
      ax = fig.add_subplot(111, aspect='equal')
      b,beta,s1,s2,s3,mago,o1,o2,o3,mag_colat,theta,phi,X,O,Q,IOmdL=np.loadtxt(f,unpack=True,usecols=range(16))
      bmax=max(b)
      if contours is not None:
         x=b*np.cos(beta)/bmax
         y=b*np.sin(beta)/bmax
         ax.tricontour(np.concatenate((x,x)),np.concatenate((y,-y)),
                       np.concatenate( (mag_colat, mag_colat)),
                       levels=contours, linewidths=0.5, colors='k')
      if drawedge:
         e = Ellipse(xy=[0,0],width=2,height=2)
         ax.add_artist(e)
         e.set_clip_box(ax.bbox)
         e.set_alpha(1.0)
         e.set_facecolor([1.0,1.0,1.0])
         e.set_edgecolor([0.0,0.0,0.0])
         e.set_linewidth(1.0)

      bmax=max(b)
      nb=(int) (len(b)/2)**0.5
      if nb>maxrings:
         skip=(int) (nb/maxrings)
         nbf=nb/skip
         nstart=np.arange(0,nbf+1)
         istart=nstart*nstart*skip*skip*2
         iend=(nstart*skip+1)*(nstart*skip+1)*2
         ival=[]
         for ii,jj in zip(istart,iend):
            ival=np.concatenate((ival,np.arange(ii+skip/2,jj,skip)))
         ival=ival[ival<len(b)].astype(int)
         b,beta,s1,s2,s3=b[ival],beta[ival],s1[ival],s2[ival],s3[ival]
      if plotpoints:
         plt.plot(b/bmax*np.cos(beta),
                  b/bmax*np.sin(beta),'b.')
         plt.plot(b/bmax*np.cos(beta),
                  -b/bmax*np.sin(beta),'b.')
      for bb,be,s1a,s2a,s3a in zip(b/bmax,beta,s1,s2,np.abs(s3)):
         _doellipse(bb,be,s1a,s2a,s3a,size,ax)
         _doellipse(bb,-be,s1a,-s2a,s3a,size,ax)

      plt.xlim(-1.1,1.1)
      plt.ylim(-1.1,1.1)
      if drawaxes==False:
         plt.axis('off')
 
   plt.savefig(f.rsplit('/', 1)[-1]+'.pdf',bbox_inches='tight')
   plt.show()
   return fig

def _main():
   flist=[]
#
# set the defaults
#
   plotpoints=False
   skip=False
   maxrings=15
   size=0.08
   contours=[15,30,45,60,75,90,105,120,135,150,165,180]
   outfile=None
   drawedge=True
   drawaxes=False
#
# parse the command line
#   
   for i,s in enumerate(sys.argv[1:]):
      if s=='--points':
         plotpoints=True
      elif s=='--nopoints':
         plotpoints=False
      elif s=='--axes':
         drawaxes=True
      elif s=='--noaxes':
         drawaxes=False
      elif s=='--edge':
         drawedge=True
      elif s=='--noedge':
         drawedge=False
      elif s=='--nocontours':
         contours=None
      elif s=='--contours':
         contours=[float(i) for i in sys.argv[i+2].split(',')]
         skip=True
      elif s=='--maxrings':
         maxrings=int(sys.argv[i+2])
         skip=True
      elif s=='--size':
         size=float(sys.argv[i+2])
         skip=True
      elif s=='--outfile':
         outfile=sys.argv[i+2]
         skip=True
      else:
         if skip:
            skip=False
         else:
            flist.append(s)
#
# generate the plots
#

   f=plotfiles(flist,contours=contours,plotpoints=plotpoints,
               maxrings=maxrings,size=size,drawedge=drawedge,drawaxes=drawaxes)

#
# output the final plot 
#
   if outfile is not None:
      f.savefig(outfile)


if __name__ == '__main__':
    _main()

