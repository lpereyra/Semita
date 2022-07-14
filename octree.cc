#define _POSIX_C_SOURCE 200112L 
#define _XOPEN_SOURCE 600

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>

#include "variables.hh"
#include "properties.hh"
#include "octree.hh"

#define KERNEL_LENGTH 10000

static struct NODE *nodes;

static struct NODE *last;

static TYPE_INT numnodestotal;        /* total number of nodes */
static TYPE_INT MaxNodes;

static TYPE_REAL xmin[3],xmax[3];

static TYPE_INT N;

static float knlrad [KERNEL_LENGTH+1],
             knlpot [KERNEL_LENGTH+1];

/*************************************************************************/

static void force_setkernel(void)
{
  int i;
  double u;
  double u2;
  double u3;
  double u4;
  double u5;

  for(i = 0; i <= KERNEL_LENGTH; i++){
    u = (double)i/(double)KERNEL_LENGTH;

    knlrad[i] = u;

    if(u <= 0.5){
      u2 = u*u;
      u4 = u2*u2;
      u5 = u4*u;
      knlpot[i]  = 16.0/3.0*u2;
      knlpot[i] -= 48.0/5.0*u4;
      knlpot[i] += 32.0/5.0*u5;
      knlpot[i] -= 14.0/5.0;
    }else{
      u2 = u*u;
      u3 = u2*u;
      u4 = u2*u2;
      u5 = u4*u;
      knlpot[i]  = 1.0/15.0/u;
      knlpot[i] += 32.0/3.0*u2;
      knlpot[i] -= 16.0*u3;
      knlpot[i] += 48.0/5.0*u4;
      knlpot[i] -= 32.0/15.0*u5;
      knlpot[i] -= 16.0/5.0;
    }
  }
}

static void add_particle_props_to_node(struct NODE *no, TYPE_REAL *pos, TYPE_INT p)
{
  int i;
  for( i = 0 ; i < 3 ; i++)
    no->s[i] += (TYPE_REAL)cp.Mpart * (pos[3*p+i] - no->center[i]);

  no->mass += (TYPE_REAL)cp.Mpart;
}

static void force_setupnonrecursive(struct NODE *no)
{
  int i;
  struct NODE *nn;
  
  if(last)
    last->next = no;

  last = no;
  
  for(i = 0; i < 8; i++){
    nn = no->suns[i];
    if(nn != NULL)
      force_setupnonrecursive(nn);
  }
}

/****************************************/

/* packs the particles of group 'gr' into into BH-trees */
extern int force_treebuild(struct grup_data *Prop, float thetamax)
{
  TYPE_INT i;
  int j;
  int    subp,subi,p,subnode,fak;
  TYPE_REAL  length;
  TYPE_REAL  dx,dy,dz;
  struct NODE *nfree,*th,*nn,*ff;

  N = Prop->npart;

  nfree = nodes;
  numnodestotal = 0;

  for(j = 0 ; j < 3 ; j++)                        /* find enclosing rectangle */
    xmin[j] = xmax[j] = Prop->pos[j];

  for(i = 0 ; i < N ; i++)
    for(j = 0 ; j < 3 ; j++)
    {
      if(Prop->pos[3*i+j] > xmax[j]) 
        xmax[j] = Prop->pos[3*i+j];
      if(Prop->pos[3*i+j] < xmin[j]) 
        xmin[j] = Prop->pos[3*i+j];
    }
  
  for(j = 1 , length = xmax[0]-xmin[0] ; j < 3 ; j++)  /* determine maxmimum extension */
    if((xmax[j]-xmin[j]) > length)
      length = xmax[j]-xmin[j];
  
  length *= 1.001;

  /* insert first particle in root node */
  for(j = 0 ; j < 3 ; j++)
    nfree->center[j] = (xmax[j]+xmin[j])/2;
  nfree->len = length;
  
  /*inicializa variables*/
  nfree->father = 0;
  for(j = 0 ; j < 8 ; j++)
    nfree->suns[j] = 0;
  nfree->partind = 0;

  nfree->mass = (TYPE_REAL)cp.Mpart;

  for(j = 0 ; j < 3 ; j++)
    nfree->s[j] = (TYPE_REAL)cp.Mpart*(Prop->pos[j] - nfree->center[j]);
  
	/* initializes the variable 
	 * that points to the sibling */ 
  nfree->sibling = 0;

  numnodestotal++; nfree++;
  
  if(numnodestotal >= MaxNodes)
  {
    fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
            numnodestotal,__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* insert all other particles */
  for(i = 1 ; i < N ; i++)
  {
    th = nodes;
  
    while(1)
    {
      add_particle_props_to_node(th,Prop->pos,i);

      if(th->partind >= 0)
        break;
    
      for(j = 0 , subnode = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
        if(Prop->pos[3*i+j] > th->center[j])
          subnode += fak;

      nn = th->suns[subnode];
      if(nn != NULL)
        th = nn;
      else
        break;
    }
      
    if(th->partind >= 0)  /* cell is occcupied with one particle */
    {
      while(1)
      {
        p = th->partind;

        for( j = 0 , subp = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(Prop->pos[3*p+j] > th->center[j])
            subp += fak;

        nfree->father = th;
        
        for( j = 0 ; j < 8 ; j++)
          nfree->suns[j] = 0;

        nfree->sibling = 0;
        
        nfree->len = th->len/2;
    
          for(j = 0 ; j < 3 ; j++)
          nfree->center[j] = th->center[j];

        for(j = 0 ; j < 3 ; j++)
          if(Prop->pos[3*p+j] > nfree->center[j])
            nfree->center[j] += nfree->len/2;
          else
            nfree->center[j] -= nfree->len/2;

        nfree->partind = p;

        nfree->mass = (TYPE_REAL)cp.Mpart;

        for(j = 0 ; j < 3 ; j++)
          nfree->s[j] = (TYPE_REAL)cp.Mpart*(Prop->pos[3*p+j] - nfree->center[j]);
        
        th->partind = -1;
        th->suns[subp] = nfree;
      
        numnodestotal++; nfree++;
        if(numnodestotal >= MaxNodes)
        {
          fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
                  numnodestotal,__FILE__,__LINE__);
          fprintf(stderr,"i=%d\n", i);
          exit(EXIT_FAILURE);
        }

        for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(Prop->pos[3*i+j] > th->center[j])
            subi += fak;

        if(subi == subp)   /* the new particle lies in the same sub-cube */
        {
          th = nfree-1;
          add_particle_props_to_node(th,Prop->pos,i);
        }
        else
          break;
      }
    }
      
    for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
      if(Prop->pos[3*i+j] > th->center[j])
        subi += fak;
      
    nfree->father = th;
      
    for(j = 0 ; j < 8 ; j++)
      nfree->suns[j] = 0;

    nfree->sibling = 0;

    nfree->len = th->len/2;

    for(j = 0 ; j < 3 ; j++)
      nfree->center[j] = th->center[j];

    for(j = 0 ; j < 3 ; j++)
      if(Prop->pos[3*i+j] > nfree->center[j])
        nfree->center[j] += nfree->len/2;
      else
        nfree->center[j] -= nfree->len/2;

    nfree->mass = (TYPE_REAL)cp.Mpart;

    for(j = 0 ; j < 3 ; j++)
      nfree->s[j] = (TYPE_REAL)cp.Mpart*(Prop->pos[3*i+j] - nfree->center[j]);

    nfree->partind = i;
    th->suns[subi] = nfree;
      
    numnodestotal++; nfree++;

    if(numnodestotal >= MaxNodes)
    {
      fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
              numnodestotal,__FILE__,__LINE__);
      fprintf(stderr,"i=%d\n", i);
      exit(EXIT_FAILURE);
    }
  }

  /* now finish-up center-of-mass and quadrupol computation */
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
  {
    for(j = 0; j < 3; j++)
      th->s[j] /= th->mass;
      
    if(th->partind < 0)   /* cell contains more than one particle */
    {
      dx = th->s[0];
      dy = th->s[1];
      dz = th->s[2];
    
      th->oc  = (TYPE_REAL)sqrt(dx*dx + dy*dy + dz*dz);
      th->oc += th->len/(thetamax); 
      th->oc *= th->oc;     /* used in cell-opening criterion */
    }

    th->s[0] += th->center[0];
    th->s[1] += th->center[1];
    th->s[2] += th->center[2];
 
    for(j = 7, nn = 0; j >= 0; j--)    /* preparations for non-recursive walk */
    {
      if(th->suns[j])
      {
        th->suns[j]->sibling = nn;
        nn = th->suns[j];
      }
    }
  }

  last = 0;
  force_setupnonrecursive(nodes);    /* set up non-recursive walk */
  last->next = 0;
  
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
    if(!(th->sibling))
    {
      ff = th;
      nn = ff->sibling;

      while(!nn)
      {
        ff = ff->father;
        if(!ff)
        break;
        nn = ff->sibling;
      }
  
      th->sibling = nn;
    }

  return numnodestotal;
}

/****************************************************************/
extern TYPE_REAL force_treeevaluate_potential(const TYPE_REAL pos [])
{
  struct NODE *no;
  TYPE_REAL r2,dx,dy,dz,r,u,h,ff;
  TYPE_REAL wp;
  TYPE_REAL h_inv;
  TYPE_REAL local_pot;
  int ii;

  h = 2.8*(TYPE_REAL)cp.soft;
  h_inv = 1.0 / h;

  local_pot = 0.;

  no = nodes;
  while(no){
    dx = no->s[0];     
    dy = no->s[1];     
    dz = no->s[2];
    dx -= pos[0];     
    dy -= pos[1];     
    dz -= pos[2];

    r2 = dx*dx + dy*dy + dz*dz;

    if(no->partind >= 0){   /* single particle */
      r = (TYPE_REAL)sqrt(r2);  
     
      u = r * h_inv;

      if(u >= 1){
        local_pot -= no->mass/r;
      }else{
        ff = u*KERNEL_LENGTH;
        ii = (int)ff;
        ff -= knlrad[ii]*KERNEL_LENGTH;
        wp = knlpot[ii] + (knlpot[ii+1] - knlpot[ii])*ff;
        
        local_pot += no->mass*h_inv*wp;
      }
      no = no->sibling;
    }else{
      if(r2 < no->oc){
        no = no->next;  /* open cell */
      }else{
        r = (TYPE_REAL)sqrt(r2);  
        u = r*h_inv;
    
        if(u >= 1){  /* ordinary quadrupol moment */
          local_pot += -no->mass/r;
        }else{    /* softened monopole moment */
          ff = u*KERNEL_LENGTH;
          ii = (int)ff;
          ff -= knlrad[ii]*KERNEL_LENGTH;
          wp = knlpot[ii] + (knlpot[ii+1]-knlpot[ii])*ff;

          local_pot += no->mass*h_inv*wp;
        }
        no = no->sibling;
      }
    }
  }

  return local_pot;
}

/***ALLOCATACION Y FREE***/
extern void force_treeallocate(int maxnodes) 
{
  MaxNodes = maxnodes;
  nodes = (struct NODE *) malloc(MaxNodes*sizeof(struct NODE));
  assert(nodes!=NULL);
  force_setkernel();
}

extern void force_treefree(void)
{
  free(nodes);
}
/************************************************************************/


