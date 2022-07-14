#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.hh"
#include "grid.hh"

extern void grid_init(void)
{
  unsigned long nalloc = (grid.ngrid*grid.ngrid*grid.ngrid) + 1;

  fprintf(stdout,"Allocating %.5f gb\n",
          (double)(nalloc*sizeof(int))/1024.0/1024.0/1024.0);
  fflush(stdout);

  grid.icell	= (TYPE_INT *) malloc(nalloc*sizeof(TYPE_INT));
  assert(grid.icell != NULL);
  for(unsigned long i = 0; i<=nalloc; i++)
    grid.icell[i] = cp.npart;

  return;
}

/*  This function brings the particles 
 *  into the same order as
 *  the sorted auxiliary. 
 */
static void reorder_particles(TYPE_INT *tmp_Id)
{ 
  TYPE_INT i;

  for(i=0;i<cp.npart;i++)
  {
    if(tmp_Id[i] != i)
    {

#ifdef COLUMN

      TYPE_REAL Px_source  = P.x[i];
      TYPE_REAL Py_source  = P.y[i];
      TYPE_REAL Pz_source  = P.z[i];
      #ifdef STORE_VELOCITIES
        TYPE_REAL Pvx_source = P.vx[i];
        TYPE_REAL Pvy_source = P.vy[i];
        TYPE_REAL Pvz_source = P.vz[i];
      #endif
      #ifdef STORE_IDS
        TYPE_INT Pid_source  = P.id[i];
      #endif

#else

      struct particle_data P_source = P[i];

#endif

      TYPE_INT  idsource = tmp_Id[i];
      TYPE_INT  dest     = tmp_Id[i];

      while(1)
      {
#ifdef COLUMN
        TYPE_REAL Px_save = P.x[dest];
        TYPE_REAL Py_save = P.y[dest];
        TYPE_REAL Pz_save = P.z[dest];
        #ifdef STORE_VELOCITIES
          TYPE_REAL Pvx_save = P.vx[dest];
          TYPE_REAL Pvy_save = P.vy[dest];
          TYPE_REAL Pvz_save = P.vz[dest];
        #endif
        #ifdef STORE_IDS
          TYPE_INT  Pid_save = P.id[dest];
        #endif
#else
        struct particle_data P_save = P[dest];
#endif

	      TYPE_INT idsave = tmp_Id[dest];

#ifdef COLUMN
        P.x[dest]  = Px_source;
        P.y[dest]  = Py_source;
        P.z[dest]  = Pz_source;
        #ifdef STORE_VELOCITIES
          P.vx[dest] = Pvx_source;
          P.vy[dest] = Pvy_source;
          P.vz[dest] = Pvz_source;
        #endif
        #ifdef STORE_IDS
          P.id[dest] = Pid_source;
        #endif
#else
        P[dest] = P_source;
#endif

        tmp_Id[dest] = idsource;

        if(dest == i)  break;

#ifdef COLUMN
        Px_source  = Px_save;
        Py_source  = Py_save;
        Pz_source  = Pz_save;
  	#ifdef STORE_VELOCITIES
          Pvx_source = Pvx_save;
          Pvy_source = Pvy_save;
          Pvz_source = Pvz_save;
	#endif
	#ifdef STORE_IDS
          Pid_source = Pid_save;
	#endif

#else
        P_source = P_save;
#endif

        idsource = idsave;
  	    dest = idsource;
      } // close while
    }  // close if

  } // close for

  return;
}

extern void grid_build(void)
{
  const long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;
  TYPE_INT i, j;
  TYPE_INT *tmp_Id; // Auxiliary array
  long ix, iy, iz, ibox;
  double fac;

  fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);
  tmp_Id = (TYPE_INT *) malloc(grid.nobj*sizeof(TYPE_INT));
  assert(tmp_Id != NULL);

  for(i = 0; i < grid.nobj; i++)
  {
#ifdef COLUMN
    ix = (long)((double)P.x[i]*fac);
    iy = (long)((double)P.y[i]*fac);
    iz = (long)((double)P.z[i]*fac);
#else
    ix = (long)((double)P[i].pos[0]*fac);
    iy = (long)((double)P[i].pos[1]*fac);
    iz = (long)((double)P[i].pos[2]*fac);
#endif

    #ifdef PERIODIC
      ibox = igrid(( (ix >= (long)grid.ngrid) ? ix-(long)grid.ngrid : ( (ix<0) ? ix + (long)grid.ngrid : ix ) ),\
                   ( (iy >= (long)grid.ngrid) ? iy-(long)grid.ngrid : ( (iy<0) ? iy + (long)grid.ngrid : iy ) ),\
                   ( (iz >= (long)grid.ngrid) ? iz-(long)grid.ngrid : ( (iz<0) ? iz + (long)grid.ngrid : iz ) ),\
                   (long)grid.ngrid);
    #else
      ibox = igrid(( (ix >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (ix<0) ? 0 : ix ) ),\
                   ( (iy >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (iy<0) ? 0 : iy ) ),\
                   ( (iz >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (iz<0) ? 0 : iz ) ),\
                   (long)grid.ngrid);
    #endif

    assert(ibox>=0 && ibox<nalloc);

    tmp_Id[i] = grid.icell[ibox];
    grid.icell[ibox] = i;
  }

  // Sorted tmp_Id

  j = 0;
  for(ibox = 0; ibox < nalloc; ibox++)
  {      
    i = grid.icell[ibox];
    grid.icell[ibox] = j;

    while(i != cp.npart)
    {
      unsigned long tmp = tmp_Id[i];
      tmp_Id[i] = j;
      j++;	
      i = tmp;
    }
  }

  //Ghost cell
  grid.icell[nalloc] = cp.npart;

  reorder_particles(tmp_Id);

  free(tmp_Id);

  return;
}


extern void grid_free(void)
{
  if(grid.icell!=NULL) free(grid.icell);

  return;
}
