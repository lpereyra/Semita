/*  file allocate.c
 *  Routines for allocating particle
*/

#include <stdlib.h>
#include <stdio.h>
#include "variables.hh"
#include "allocate.hh"

/*
 *  Allocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
#ifdef COLUMN
  extern int allocate_particles(struct particle_data  *Q, const TYPE_INT size)
#else
  extern int allocate_particles(struct particle_data **Q, const TYPE_INT size)
#endif
{

#ifdef COLUMN
  Q->x  = NULL;
  Q->y  = NULL;
  Q->z  = NULL;

  Q->x  = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));
  Q->y  = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));
  Q->z  = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));

  if(!Q->x || !Q->y || !Q->z) 
  {
    fprintf(stderr, "cannot allocate pos particles\n" );
    return(0);
  }

  #ifdef STORE_VELOCITIES
  Q->vx = NULL;
  Q->vy = NULL;
  Q->vz = NULL;

  Q->vx = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));
  Q->vy = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));
  Q->vz = (TYPE_REAL *) malloc(size*sizeof(TYPE_REAL));

  if(!Q->vx || !Q->vy || !Q->vz) 
  {
    fprintf(stderr, "cannot allocate vel particles\n" );
    return(0);
  }
  #endif

  #ifdef STORE_IDS
  Q->id = NULL;

  Q->id = (TYPE_INT  *) malloc(size*sizeof(TYPE_INT));
  
  if(!Q->id) 
  {
    fprintf(stderr, "cannot allocate ids particles\n" );
    return(0);
  }
  #endif

  Q->sub = NULL;

  Q->sub = (TYPE_INT  *) malloc(size*sizeof(TYPE_INT));
  
  if(!Q->sub) 
  {
    fprintf(stderr, "cannot allocate sub particles\n" );
    return(0);
  }

#else

  *Q = NULL;

  *Q = (struct particle_data *) malloc(size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot allocate particles\n" );
    return(0);
  }    

#endif

  return ( 1 );
}

/*
 *  Deallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
#ifdef COLUMN
  extern void free_particles(struct particle_data  *Q)
#else
  extern void free_particles(struct particle_data **Q)
#endif
{

#ifdef COLUMN
  if(Q->x)  free(Q->x);
  if(Q->y)  free(Q->y);
  if(Q->z)  free(Q->z);
  #ifdef STORE_VELOCITIES
  if(Q->vx) free(Q->vx);
  if(Q->vy) free(Q->vy);
  if(Q->vz) free(Q->vz);
  #endif   
  #ifdef STORE_IDS
  if(Q->id) free(Q->id);
  #endif   
  if(Q->sub) free(Q->sub);
#else
  if(*Q) free(*Q);
#endif  

}

/*
 *  Reallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
#ifdef COLUMN
  extern int reallocate_particles(struct particle_data  *Q, const TYPE_INT size)
#else
  extern int reallocate_particles(struct particle_data **Q, const TYPE_INT size)
#endif
{

#ifdef COLUMN
  Q->x  = (TYPE_REAL *) realloc(Q->x,size*sizeof(TYPE_REAL));
  Q->y  = (TYPE_REAL *) realloc(Q->y,size*sizeof(TYPE_REAL));
  Q->z  = (TYPE_REAL *) realloc(Q->z,size*sizeof(TYPE_REAL));

  if(!Q->x || !Q->y || !Q->z) 
  {
    fprintf(stderr, "cannot reallocate pos particles\n" );
    return(0);
  }

  #ifdef STORE_VELOCITIES
  Q->vx = (TYPE_REAL *) realloc(Q->vx,size*sizeof(TYPE_REAL));
  Q->vy = (TYPE_REAL *) realloc(Q->vy,size*sizeof(TYPE_REAL));
  Q->vz = (TYPE_REAL *) realloc(Q->vz,size*sizeof(TYPE_REAL));

  if(!Q->vx || !Q->vy || !Q->vz) 
  {
    fprintf(stderr, "cannot reallocate vel particles\n" );
    return(0);
  }
  #endif

  #ifdef STORE_IDS
  Q->id = (TYPE_INT  *) realloc(Q->id,size*sizeof(TYPE_INT));
  
  if(!Q->id) 
  {
    fprintf(stderr, "cannot reallocate ids particles\n" );
    return(0);
  }
  #endif

  Q->sub = (TYPE_INT  *) realloc(Q->sub,size*sizeof(TYPE_INT));
  
  if(!Q->sub) 
  {
    fprintf(stderr, "cannot reallocate sub particles\n" );
    return(0);
  }

#else

  *Q = (struct particle_data *) realloc(*Q,size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot reallocate particles\n" );
    return(0);
  }    

#endif

  return ( 1 );
}


