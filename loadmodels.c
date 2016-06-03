#include "models.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define ALLOCSTEP 20


void
initializenode(data_node *nodeptr) {
  nodeptr->description=NULL;
  nodeptr->data_description=NULL;
  nodeptr->ndata=0;
  nodeptr->data=NULL;

  nodeptr->nchildren=0;
  nodeptr->children=NULL;
  nodeptr->parent=NULL;

}

int
parse_trj_line(char *buffer, trj_line *trjptr) {
  int res;
  int tr_idx, mu_idx, phi_idx;
  double mu, theta, cosphi, phi, angbk, angbmkkmz;

  res=sscanf(buffer,"%d %d %d %lg %lg %lg %lg %lg %lg",&tr_idx,&mu_idx,&phi_idx,&mu,&theta,&cosphi,&phi,&angbk,&angbmkkmz);

  if ( res == 9 ) {
    trjptr->tr_idx=tr_idx;
    trjptr->mu_idx=mu_idx;
    trjptr->phi_idx=phi_idx;
    trjptr->mu=mu;
    trjptr->theta=theta;
    trjptr->phi=phi;
    trjptr->angbk=angbk;
    trjptr->angbmkkmz=angbmkkmz;
  }
  return(res);
}


data_node *malloc_data_node(int ndata_node) {
  int i;
  data_node *res;
  res = (data_node *) malloc(ndata_node*sizeof(data_node));
  if (res!=NULL) {
    for (i=0;i<ndata_node;i++) {
      initializenode(res+i);
    }
  }

  return(res);
}

data_node *realloc_data_node(data_node *nodeptr, int ndata_node) {
  return( (data_node *) realloc((void *) nodeptr, ndata_node));
}

int
loadtrjfile(char *filename, data_node *nodeptr) {
  char buffer[200], *stripped_filename;
  FILE *in;
  int i, j;
  trj_line line_hld;
  data_node *childptr, *grandchildptr;

  nodeptr->description="TRJ Node";
  nodeptr->data_description="Magnetic Colatitude [degrees]";

  nodeptr->ndata=1;
  if ( (nodeptr->data = (double *)malloc(sizeof(double))) == NULL ) {
    return(-1);
  }
  stripped_filename=strrchr(filename,'/');

  if (stripped_filename==NULL) {
    sscanf(filename,"%lg",nodeptr->data);
  } else {
    sscanf(stripped_filename,"/%lg",nodeptr->data);
  }
  
  if ((in=fopen(filename,"r"))==NULL) {
    return(-1);
  }

  // go to the end of the file to get the total number of mu and phi points
  while ( fgets(buffer,200,in) ) {
    if (buffer[0]=='#') continue;
    if ( parse_trj_line(buffer,&line_hld) != 9) break;
  }
  fclose(in);
  
  // allocate the children!
  nodeptr->nchildren=line_hld.mu_idx+1;
  
  if ((nodeptr->children = malloc_data_node(nodeptr->nchildren))==NULL) {
    return(-1);
  }

  for (i=0,childptr=nodeptr->children;i<nodeptr->nchildren;i++,childptr++) {
    childptr->description="Altitude Node";
    childptr->data_description="Mu";
    childptr->nchildren=line_hld.phi_idx+1;
    if ((childptr->children = malloc_data_node(childptr->nchildren))==NULL) {
      return(-1);
    }
    childptr->ndata=1;
    if ( (childptr->data=(double *) malloc(sizeof(double))) == NULL ) {
      return(-1);
    }
    childptr->parent=nodeptr;
    for (j=0,grandchildptr=childptr->children;j<childptr->nchildren;j++,grandchildptr++) {
      grandchildptr->description="Azimuth Node";
      grandchildptr->data_description="Phi [degrees]";
      grandchildptr->ndata=1;
      if ( (grandchildptr->data=(double *) malloc(sizeof(double))) == NULL ) {
	return(-1);
      }
      grandchildptr->parent=childptr;
    }
  }

  // fill up the children

  if ((in=fopen(filename,"r"))==NULL) {
    return(-1);
  }

  // go to the end of the file to get the total number of mu and phi points
  while ( fgets(buffer,200,in) ) {
    if (buffer[0]=='#') continue;
    if ( parse_trj_line(buffer,&line_hld) != 9) break;
    nodeptr->children[line_hld.mu_idx].data[0]=line_hld.mu;
    nodeptr->children[line_hld.mu_idx].children[line_hld.phi_idx].data[0]=line_hld.phi;
  }

  fclose(in);
  
}


int 
loadgnufile(char *filename, data_node *nodeptr) {
  char buffer[200], gnufilename[200];
  char *stripped_filename;  
  int nalloc, ntotal;
  double *energy_array;
  data_node *childptr, *grandchildptr, *ggchildptr;
  int i, id, j, jd, k, kd;
  FILE *in;

  strcpy(gnufilename,filename);
  stripped_filename=strrchr(gnufilename,'.');

  if (stripped_filename==NULL) {
    return(-1);
  } else {
    strcpy(stripped_filename,".gnu");
  }

  ntotal=0;
  nalloc=ALLOCSTEP;
  if ( (energy_array = (double *) malloc(sizeof(double)*nalloc)) == NULL ) {
    return(-1);
  }
  if ((in=fopen(gnufilename,"r"))==NULL) {
    return(-1);
  }
  while ( fgets(buffer,200,in) ) {
    if (buffer[0]=='#') continue;
    sscanf(buffer,"%*d %*g %lg",energy_array+ntotal);
    energy_array[ntotal]=log(energy_array[ntotal]);
    ntotal++;
    if (ntotal>=nalloc) {
      nalloc+=ALLOCSTEP;
      if ( (energy_array = (double *) realloc((void *) energy_array, sizeof(double)*nalloc)) == NULL ) {
	return(-1);
      }
    }
  }
  fclose(in);
  if (ntotal<nalloc) {
    if ( (energy_array = (double *) realloc((void *) energy_array, sizeof(double)*ntotal)) == NULL ) {
      return(-1);
    }
  }

  // build an extra layer onto the bottom of the tree!
  id=nodeptr->nchildren;
  for (i=0,childptr=nodeptr->children;i<id;i++,childptr++) {
    jd=childptr->nchildren;
    for (j=0,grandchildptr=childptr->children;j<jd;j++,grandchildptr++) {
      grandchildptr->nchildren=ntotal;
      if ((grandchildptr->children = malloc_data_node(ntotal))==NULL) {
	return(-1);
      }
      kd=ntotal;
      for (k=0,ggchildptr=grandchildptr->children;k<kd;k++,ggchildptr++) {
	ggchildptr->description="Frequency, Intensities";
	ggchildptr->data_description=" nu [Hz], total, X, O [natural logs]";
	ggchildptr->ndata=NDATA;
	if (( ggchildptr->data=(double *) malloc(NDATA*sizeof(double))) == NULL) {
	  return(-1);
	}
	ggchildptr->data[0]=energy_array[k];
	ggchildptr->data[1]=ggchildptr->data[2]=ggchildptr->data[3]=0.0;
	ggchildptr->parent=grandchildptr;
     }
    }
  }

  free((void *) energy_array);
  return(0);
}

int 
loadintfile(char *filename, data_node *nodeptr) {
  char buffer[200], intfilename[200];
  char *stripped_filename;  
  int nalloc, ntotal;
  double *dataptr;
  data_node *childptr, *grandchildptr, *ggchildptr;
  int i, id, j, jd, k, kd, trdx;
  double t, x, o;
  FILE *in;

  strcpy(intfilename,filename);
  stripped_filename=strrchr(intfilename,'.');

  if (stripped_filename==NULL) {
    return(-1);
  } else {
    strcpy(stripped_filename,".int");
  }

  if ((in=fopen(intfilename,"r"))==NULL) {
    return(-1);
  }

  id=nodeptr->nchildren;
  jd=nodeptr->children->nchildren;
  kd=nodeptr->children->children->nchildren;
  while ( fgets(buffer,200,in) ) {
    if (buffer[0]=='#') continue;
    sscanf(buffer,"%d %d %lg %lg %lg",&trdx,&k,&t,&x,&o);
    i=trdx/jd;
    j=trdx%jd;
    if ((i<id) && (k<kd)) {
      dataptr=nodeptr->children[i].children[j].children[k].data;
      dataptr[1]=log(t);
      dataptr[2]=log(x);
      dataptr[3]=log(o);
    }
  }
  fclose(in);
  return(0);
}



void
printtree_reset(data_node *nodeptr, int reset) {
  int i;
  static int depth, ndescendents;

  if (reset) {
    depth=0; ndescendents=0;
  } else {
    depth++; ndescendents++;
  }
  
  printf("Node description (%d): %s\n",depth,nodeptr->description);
  if (nodeptr->data_description!=NULL) {
    printf("Data description: %s\n",nodeptr->data_description);
  }

  for (i=0;i<nodeptr->ndata;i++) {
    printf("data[%d] = %g\n",i,nodeptr->data[i]);
  }

  printf("The node has %d children\n",nodeptr->nchildren);

  if (nodeptr->nchildren>0) {
    printf("{\n");
    for (i=0;i<nodeptr->nchildren;i++) {
      printtree_reset(nodeptr->children+i,0);
    }
    printf("}\n");
  }
  depth--;
  if (reset) {
    printf("Total descendents %d\n",ndescendents);
  }
}

void
printtree(data_node *nodeptr) {
  printtree_reset(nodeptr,1);
}

void
evaltree(data_node *nodeptr, double *args, int nargs, double *res) {
  double res1hld[NDATA], res2hld[NDATA], d1, d2, w1, w2;
  int i;
  double ascend;
  data_node *p1, *p2, *mid;
  double value=*args;

  if ((nargs==0) || (nodeptr->nchildren==0)) {
    for (i=0;i<nodeptr->ndata;i++) {
      res[i]=nodeptr->data[i];
      // printf("resout[%d]=%g\n",i,res[i]);
    }
    return;
  }
  if (nodeptr->nchildren==1) {
    evaltree(nodeptr->children,args+1,nargs-1,res);
    return;
  }
  if (nodeptr->nchildren==2) {
    p1=nodeptr->children; p2=nodeptr->children+1;
  } else {
    // do binary search
    p1=nodeptr->children;
    p2=nodeptr->children+nodeptr->nchildren-1;
    ascend=p2->data[0]-p1->data[0];
    //    printf("ascend=%g\n",ascend);
    if ( (value-p2->data[0])*ascend > 0 ) {
      p2;
      p1=p2-1;
    } else if ( (value-p1->data[0])*ascend < 0 ) {
      p2=p1+1;
    } else {
      while (p1<p2) {
	mid=p1+(p2-p1)/2;
	if ( (value-mid->data[0])*ascend < 0 ) {
	  p2=mid;
	} else {
	  p1=mid+1;
	}
        // printf("value=%g %g %g\n",value,p1->data[0],p2->data[0]);
      }
      p1--;
      // printf("Done!\n");
    }
  }

  evaltree(p1,args+1,nargs-1,res1hld);
  evaltree(p2,args+1,nargs-1,res2hld);
  d1=p1->data[0];
  d2=p2->data[0];

  w1=(d2-value)/(d2-d1);
  w2=(value-d1)/(d2-d1);
  //  printf("value %s=%g %g %g %g %g\n",p1->data_description,value,d1,d2,w1,w2);
    
  for (i=0;i<NDATA;i++) {
    res[i]=w1*res1hld[i]+w2*res2hld[i];
  }

  return;



}
#if 0
int
main(int argc, char *argv[]) {
  data_node parent_node;
  int i;
  double res[NDATA];
  double args[]={85,0.560472,46.2245, 40.0};

  if (argc<2) {
    printf("Format:\n\n   loadmodels [ TRJ_FILES ]\n");
    return(-1);
  }

  initializenode(&parent_node);

  parent_node.description="Parent Node";

  parent_node.parent=NULL;
  parent_node.nchildren=argc-1;
  
  if ((parent_node.children = malloc_data_node(parent_node.nchildren))==NULL) {
    return(-1);
  }
  
  for (i=0;i<parent_node.nchildren;i++) {
    loadtrjfile(argv[i+1],parent_node.children+i);
    loadgnufile(argv[i+1],parent_node.children+i);
    loadintfile(argv[i+1],parent_node.children+i);
    parent_node.children[i].parent=&parent_node;
  }

  evaltree(&parent_node,args,4, res);
  for (i=0;i<NDATA;i++) {
    printf("res[%d]=%g\n",i,res[i]);
  }
  //  printtree(&parent_node);
    
    return(0);
}
#endif



