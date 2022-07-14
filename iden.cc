#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "variables.hh"
#include "properties.hh"
#include "allocate.hh"
#include "leesnap.hh"
#include "grid.hh"
#include "iden.hh"
#include "io.hh"

#define DIV_CEIL(x,y) (x+y-1)/y

static TYPE_INT nstep;
static struct iden_st iden;
static struct temporary Temp;
static TYPE_INT **gr;
#ifdef LOCK
  static omp_lock_t *lock;
#endif
 
extern inline TYPE_INT Find(TYPE_INT i, TYPE_INT * __restrict__ ar)
{
  if(i != ar[i])
    ar[i] = Find(ar[i],ar);
 
  return ar[i];
}

static inline void Union(TYPE_INT u, TYPE_INT v, TYPE_INT * __restrict__ ar)
{
 
  TYPE_INT z;

  while(ar[u] != ar[v])
  { 
      if(ar[u] < ar[v])
      {
#ifdef LOCK
          if(u == ar[u])
          {
            omp_set_lock(&(lock[u]));
            z = 0;

            if(u == ar[u])
            {
              ar[u] = ar[v];  
              z = 1;
            } 

            omp_unset_lock(&(lock[u]));
            if(z==1) break;             

          }
#else
          if(u == ar[u])
          {
            ar[u] = ar[v];  
            break;             
          }
#endif
          
          z = ar[u];   
          ar[u] = ar[v];
          u = z;

      }else{
#ifdef LOCK
          if(v == ar[v])
          {
            omp_set_lock(&(lock[v]));
            z = 0;

            if(v == ar[v])
            {
              ar[v] = ar[u];   
              z = 1;
            }

            omp_unset_lock(&(lock[v]));
            if(z == 1) break;            
          }
#else
          if(v == ar[v])
          {
            ar[v] = ar[u];   
            break;             
          }
#endif

          z = ar[v];   
          ar[v] = ar[u];   
          v = z;

      }
  }

}

static void busv(const TYPE_INT ic)
{
  long ixc, iyc, izc, ibox;
  TYPE_INT i, niv;
  TYPE_REAL xx, yy, zz;

#ifdef COLUMN

  ixc  = (long)(P.x[ic]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(P.y[ic]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(P.z[ic]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));

#else

  ixc  = (long)(P[ic].pos[0]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(P[ic].pos[1]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(P[ic].pos[2]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));

#endif

  #ifndef PERIODIC
    for(long ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
  #else
    for(long ixx = ixc-1; ixx <= ixc+1; ixx++)
  #endif
  {
    #ifndef PERIODIC
      for(long iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
    #else
      for(long iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
    #endif
    {
      #ifndef PERIODIC
        for(long izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
      #else
        for(long izz = izc-1 ; izz <= izc+1 ; izz++)
      #endif
      {

      	#ifdef PERIODIC
          ibox = igrid(( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) ),\
                       ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ),\
                       ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) ),\
                       (long)grid.ngrid);
        #else
          ibox = igrid(ixx,iyy,izz,(long)grid.ngrid);
        #endif

        for( i = (grid.icell[ibox]   > ic ? grid.icell[ibox] : ic+1) ; \
             i < (grid.icell[ibox+1] > ic ? grid.icell[ibox+1] : ic) ; \
             i++)
        {

#ifdef COLUMN           
          xx = P.x[i] - P.x[ic];
          yy = P.y[i] - P.y[ic];
          zz = P.z[i] - P.z[ic];
#else
          xx = P[i].pos[0] - P[ic].pos[0];
          yy = P[i].pos[1] - P[ic].pos[1];
          zz = P[i].pos[2] - P[ic].pos[2];
#endif
          #ifdef PERIODIC
          xx = ( xx >  cp.lbox*0.5f ) ? xx - cp.lbox : xx ;
          yy = ( yy >  cp.lbox*0.5f ) ? yy - cp.lbox : yy ;
          zz = ( zz >  cp.lbox*0.5f ) ? zz - cp.lbox : zz ;
          xx = ( xx < -cp.lbox*0.5f ) ? xx + cp.lbox : xx ;
          yy = ( yy < -cp.lbox*0.5f ) ? yy + cp.lbox : yy ;
          zz = ( zz < -cp.lbox*0.5f ) ? zz + cp.lbox : zz ;
          #endif

          for(niv=0;niv<nstep;niv++)
          {
      	    if(xx*xx + yy*yy + zz*zz < iden.r0[niv])
            {
              Union(ic,i,gr[niv]);
            }
          }
        } /*end particles loop*/
      } /*end izz*/
    } /*end iyy*/
  } /*end ixx*/

}

static void linkedlist(TYPE_INT * __restrict__ ar)
{
  TYPE_INT i, g;
  
  Temp.ll = (TYPE_INT *) calloc(iden.nobj,sizeof(TYPE_INT));

  iden.ngrupos = 0;
  for(i=0;i<iden.nobj;i++)
  {
    ar[i] = Find(i,ar);
    assert(ar[i]>=i);
    if(Temp.ll[ar[i]] < NPARTMIN)
    {
      Temp.ll[ar[i]]++;
      if(Temp.ll[ar[i]]==NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[ar[i]] = NPARTMIN + iden.ngrupos;
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (TYPE_INT *) malloc(iden.ngrupos*sizeof(TYPE_INT));
  Temp.npgrup = (TYPE_INT *) malloc(iden.ngrupos*sizeof(TYPE_INT));

  for(i=0;i<iden.ngrupos;i++)
  {
    Temp.head[i]   = iden.nobj;
    Temp.npgrup[i] =  0;
  }

  for(i=0;i<iden.nobj;i++)
  {
    if(Temp.ll[ar[i]]>NPARTMIN)
    {
      ar[i] = Temp.ll[ar[i]] - NPARTMIN;
    }else{
#ifdef COLUMN
      P.sub[i] = ar[i] = 0;
#else
      P[i].sub = ar[i] = 0;
#endif
    }

    g = ar[i];

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }

  return;
}

static void Write_Groups(const TYPE_INT niv)
{
  TYPE_INT i,j,k,ilink,dim,npar,gn,save_sub;
  TYPE_REAL dx[3];
  char filename[200];
  FILE *pfout, *pfcentros;
  #ifdef FILE_ASCII
    FILE *pfcentros_ascii;
  #endif

  j = iden.ngrupos-1; // LE RESTO UNO POR EL GRUPO 0 PARA ESCRIBIR EN EL ARCHIVO
	npar = gn = 0;

	///////////////////////////////////////////////////////
	sprintf(filename,"%.2d_level_%.2d_fof_%.2f.bin",snap.num,niv,fof[niv]);
	pfout=fopen(filename,"w");
	fwrite(&j,sizeof(TYPE_INT),1,pfout);
	//////////////////////////////////////////////////////

	i = 0;
	if(niv == 0)
	{
	  for(ilink=1;ilink<iden.ngrupos;ilink++)
	  {
	    j = 0;
	    k = Temp.head[ilink];
	    save_sub = ilink;
	    fwrite(&save_sub,sizeof(TYPE_INT),1,pfout);
	    fwrite(&ilink,sizeof(TYPE_INT),1,pfout);
	    fwrite(&Temp.npgrup[ilink],sizeof(TYPE_INT),1,pfout);    
	
	    while(k != iden.nobj)
	    {
	#ifdef COLUMN
	      fwrite(&P.id[k],sizeof(TYPE_INT),1,pfout);
	      P.sub[k] = ilink;
	#else
	      fwrite(&P[k].id,sizeof(TYPE_INT),1,pfout);
	      P[k].sub = ilink;
	#endif
	      k = Temp.ll[k];
	      j++;
	    }
	    npar+=j;
	    gn++;
			i++;
	  }

	}else{

		cp.ngrup = iden.ngrupos-1;
	  sprintf(filename,"%.2d_level_%.2d_halos_%.2f.bin",snap.num,niv,fof[niv]);
	  pfcentros=fopen(filename,"w");
	  fwrite(&cp.ngrup,sizeof(TYPE_INT),1,pfcentros);
	  //////////////////////////////////////////////////////
	  #ifdef FILE_ASCII
	    sprintf(filename,"%.2d_level_%.2d_halos_%.2f.dat",snap.num,niv,fof[niv]);
	    pfcentros_ascii=fopen(filename,"w");
	  #endif  
	  //////////////////////////////////////////////////////
		
		// Allocate Group Properties
		Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

	  for(ilink=1;ilink<iden.ngrupos;ilink++)
	  {
	    j = 0;
	    k = Temp.head[ilink];
	#ifdef COLUMN           
	    save_sub = P.sub[k];
	#else
	    save_sub = P[k].sub;
	#endif

			Gr[i].save  = save_sub;
			Gr[i].id    = ilink;
	    Gr[i].npart = Temp.npgrup[ilink];

	    fwrite(&Gr[i].save,  sizeof(TYPE_INT), 1, pfout);
	    fwrite(&Gr[i].id,    sizeof(TYPE_INT), 1, pfout);
	    fwrite(&Gr[i].npart, sizeof(TYPE_INT), 1, pfout);    
	
	    for(dim = 0; dim < 3; dim++)
		  {
			  Gr[i].pcm[dim] = 0.;
	#ifdef STORE_VELOCITIES        
			  Gr[i].vcm[dim] = 0.;
			  Gr[i].L[dim]   = 0.;
			  Gr[i].sig[dim] = 0.;
	#endif
		  }
	    Gr[i].pos = (TYPE_REAL *) malloc(3*Gr[i].npart*sizeof(TYPE_REAL));
	#ifdef STORE_VELOCITIES        
	    Gr[i].vel = (TYPE_REAL *) malloc(3*Gr[i].npart*sizeof(TYPE_REAL));
	#endif
	
	    while(k != iden.nobj)
	    {
	
	      // cuidado con el orden {pos[i]-centro} en este caso
	#ifdef COLUMN
	      Gr[i].pos[3*j+0] = P.x[k];
	      Gr[i].pos[3*j+1] = P.y[k];
	      Gr[i].pos[3*j+2] = P.z[k];
	#ifdef STORE_VELOCITIES        
	      Gr[i].vel[3*j+0] = P.vx[k];
	      Gr[i].vel[3*j+1] = P.vy[k];
	      Gr[i].vel[3*j+2] = P.vz[k];
	#endif
	
	      dx[0] = P.x[k] - P.x[Temp.head[ilink]];
	      dx[1] = P.y[k] - P.y[Temp.head[ilink]];
	      dx[2] = P.z[k] - P.z[Temp.head[ilink]];
	#else
	      Gr[i].pos[3*j+0] = P[k].pos[0];
	      Gr[i].pos[3*j+1] = P[k].pos[1];
	      Gr[i].pos[3*j+2] = P[k].pos[2];
	#ifdef STORE_VELOCITIES        
	      Gr[i].vel[3*j+0] = P[k].vel[0];
	      Gr[i].vel[3*j+1] = P[k].vel[1];
	      Gr[i].vel[3*j+2] = P[k].vel[2];
	#endif
	
	      dx[0] = P[k].pos[0] - P[Temp.head[ilink]].pos[0];
	      dx[1] = P[k].pos[1] - P[Temp.head[ilink]].pos[1];
	      dx[2] = P[k].pos[2] - P[Temp.head[ilink]].pos[2];
	#endif
	
	      for(dim=0; dim<3; dim++)
	      {
	        #ifdef PERIODIC
	        if(dx[dim] >  cp.lbox*0.5)
	        {
	          dx[dim] -= cp.lbox;
	          Gr[i].pos[3*j+dim] -= cp.lbox;
	        }
	
	        if(dx[dim] < -cp.lbox*0.5)
	        {
	          dx[dim] += cp.lbox;
	          Gr[i].pos[3*j+dim] += cp.lbox;
	        }
	        #endif
	
	        Gr[i].pcm[dim] += dx[dim];
	#ifdef STORE_VELOCITIES        
	        Gr[i].vcm[dim] += Gr[i].vel[3*j + dim];
	#endif
	      }
	
	#ifdef COLUMN
	      fwrite(&P.id[k],sizeof(TYPE_INT),1,pfout);
	#else
	      fwrite(&P[k].id,sizeof(TYPE_INT),1,pfout);
	#endif
	      k = Temp.ll[k];
	      j++;
	    }
	    
	    assert(j == Temp.npgrup[ilink]);
	
	    for(dim = 0; dim < 3; dim++)
		  {
		  	Gr[i].pcm[dim] /= (TYPE_REAL)Gr[i].npart;
	#ifdef STORE_VELOCITIES        
		  	Gr[i].vcm[dim] /= (TYPE_REAL)Gr[i].npart;
	#endif
	
	      //	Recenter	
	      Gr[i].pcm[dim] += Gr[i].pos[dim];
	
	#ifdef PERIODIC
	      Gr[i].pcm[dim] = Gr[i].pcm[dim]<0.0f     ? cp.lbox+Gr[i].pcm[dim] : Gr[i].pcm[dim];
	      Gr[i].pcm[dim] = Gr[i].pcm[dim]>=cp.lbox ? Gr[i].pcm[dim]-cp.lbox : Gr[i].pcm[dim];
	#endif
	
	    }
	
	    properties(&Gr[i]);
	
			#ifdef CHANGE_POSITION
  			for(dim=0;dim<3;dim++)
      		Gr[i].pcm[dim] += pmin[dim];
			#endif

	    write_properties(pfcentros, Gr[i]);
	
	    #ifdef FILE_ASCII
	      fprintf(pfcentros_ascii,"%u %u ",save_sub,i);
	      write_properties_ascii(pfcentros_ascii, Gr[i]);
	    #endif
	     
	    free(Gr[i].pos);
	#ifdef STORE_VELOCITIES        
	    free(Gr[i].vel);
	#endif
	
	    npar+=j;
	    gn++;
			i++;
	  }
	
	  assert(gn == (iden.ngrupos-1));
	  fclose(pfcentros);
	  #ifdef FILE_ASCII
	    fclose(pfcentros_ascii);
	  #endif
	}
	
	fclose(pfout);
	fprintf(stdout,"num de grupos %u num de particulas en grupos %u\n",gn,npar);
	fflush(stdout);
	
  return;
}

extern void identification(void)
{
  TYPE_INT i, j, tid;

  nstep = nfrac;
  iden.nobj = cp.npart;
  iden.r0 = (double *) malloc(nstep*sizeof(double));
  gr      = (TYPE_INT **) malloc(nstep*sizeof(TYPE_INT *));

  for(j=0;j<nstep;j++)
  {
    gr[j] = (TYPE_INT *) malloc(iden.nobj*sizeof(TYPE_INT));
    iden.r0[j]  = cbrt(1./(1.+fof[j]));
    iden.r0[j] *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0f; //In Kpc

    if(iden.r0[j] <= cp.soft)
    {
      fprintf(stdout,"cambia Linking length = %f \n",iden.r0[j]);
      iden.r0[j] = cp.soft;
    }

    fprintf(stdout,"Linking length %d = %f \n",j,iden.r0[j]);
  }

  #ifdef LOCK
    lock = (omp_lock_t *) malloc(iden.nobj*sizeof(omp_lock_t));
  #endif
  for(i=0;i<iden.nobj;i++)
  {
#ifdef COLUMN           
      P.sub[i] = 0;
#else
      P[i].sub = 0;
#endif
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
    for(j=0; j<nstep; j++)
    {
      gr[j][i] = i;
    }
  }

  grid.ngrid = (long)(cp.lbox/iden.r0[0]);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  
  grid.nobj = iden.nobj;
  grid_init();
  grid_build();

  for(j=0;j<nstep;j++)
    iden.r0[j] *= iden.r0[j];

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif
    
  fprintf(stdout,"Comienza identificacion.....\n");
  fprintf(stdout,"Running on %d threads\n",NTHREADS);
  fflush(stdout);

  #ifdef LOCK
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,lock,stdout)   
  #else
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,stdout)   
  #endif
  {
    tid = omp_get_thread_num(); 

    for(i = tid*DIV_CEIL(iden.nobj,NTHREADS);
    i<(tid==NTHREADS-1 ? iden.nobj : (tid+1)*DIV_CEIL(iden.nobj,NTHREADS));
    i++)
    {
     
      if(i%1000000==0) fprintf(stdout,"%u %u %u %.4f\n",tid,i,iden.nobj,(float)i/(float)iden.nobj);

      #pragma omp task
      {
        busv(i);
      }
    }

  }  /****** TERMINA SECCION PARALELA ****************************/

  fprintf(stdout,"Sale del paralelo\n"); fflush(stdout);

  #ifdef LOCK
  for(i=0;i<iden.nobj;i++) 
    omp_destroy_lock(&(lock[i]));
  free(lock);
  #endif

  for(j=0;j<nstep;j++)
  {
    linkedlist(gr[j]);
    Write_Groups(j);

    free(gr[j]);
    free(Temp.ll);
    free(Temp.head);
    free(Temp.npgrup);
  }

  j = 0;
  for(i=0;i<cp.npart;i++)
  {
#ifdef COLUMN           
    if(P.sub[i] != 0)
    {
      TYPE_REAL Px_save  = P.x[j];
      TYPE_REAL Py_save  = P.y[j];
      TYPE_REAL Pz_save  = P.z[j];
      #ifdef STORE_VELOCITIES
        TYPE_REAL Pvx_save = P.vx[j];
        TYPE_REAL Pvy_save = P.vy[j];
        TYPE_REAL Pvz_save = P.vz[j];
      #endif
      #ifdef STORE_IDS
        TYPE_INT Pid_save  = P.id[j];
      #endif
      TYPE_INT Psub_save = P.sub[j];

      P.x[j] = P.x[i];
      P.y[j] = P.y[i];
      P.z[j] = P.z[i];
      #ifdef STORE_VELOCITIES
        P.vx[j] = P.vx[i];
        P.vy[j] = P.vy[i];
        P.vz[j] = P.vz[i];
      #endif
      #ifdef STORE_IDS
        P.id[j] = P.id[i];
      #endif
      P.sub[j] = P.sub[i];

      P.x[i]  = Px_save;
      P.y[i]  = Py_save;
      P.z[i]  = Pz_save;
      #ifdef STORE_VELOCITIES
        P.vx[i] = Pvx_save;
        P.vy[i] = Pvy_save;
        P.vz[i] = Pvz_save;
      #endif
      #ifdef STORE_IDS
        P.id[i] = Pid_save;
      #endif
      P.sub[i] = Psub_save;

      j++;
    }
#else
    if(P[i].sub != 0)
    {
      struct particle_data P_save = P[j];
      P[j] = P[i];
      P[i] = P_save;
      j++;
    }
#endif
  }
 
  cp.npart     = j;
  //if(!reallocate_particles(&P, cp.npart))  exit(1);

  free(iden.r0);
  grid_free();

  return;
}
