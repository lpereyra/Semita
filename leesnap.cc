/* file leesnap.c
 * Routine for read 
 * GADGET snapshot 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.hh"
#include "allocate.hh"
#include "leesnap.hh"
#include "colors.hh"

/* Header for the standard 
 * file GADGET format.
 */
static struct io_header header;

/* Read header GADGET format
*/
static void read_header(const char *filename)
{
  FILE *pf;
  TYPE_INT d1, d2;
#ifdef TYPE_TWO_GADGET
  TYPE_INT blocksize;
  char label[4];
#endif

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  assert(d1==256);

  // set structure cosmoparam
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[TYPE_PART];
  cp.npart_tot = cp.npart;
  cp.Mpart     = header.mass[TYPE_PART];
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

#ifdef TYPE_TWO_GADGET

  TYPE_INT k;
	int n;
  TYPE_REAL mass;
  int flag = 0;
  
  while(1)
  {
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label, sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
  
    flag = strncmp(label, "MASS", 4);
    if(flag == 0)
    {
  
      fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
      fflush(stdout);
  
      fread(&d1, sizeof(d1), 1, pf);
      for(k = 0; k < N_part_types; k++){
        for(n = 0; n < header.npart[k]; n++){
          fread(&mass, sizeof(mass), 1, pf);
          if(k == TYPE_PART){ /*ONLY KEEP DARK MATTER PARTICLES*/
            cp.Mpart = (double)mass;
          }
        }
      }
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
      break;
  
    }else{
  
      fread(&d1, sizeof(d1), 1, pf);
      fseek(pf, d1,SEEK_CUR);
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
  
    }
  
  }

  sprintf(message,"change Masa por particula = %f\n",cp.Mpart);
  RED(message);

#endif

  fclose(pf);

  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"*   Parametros de la simulacion   * \n");
  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"  Particles Number = %u \n", cp.npart);
  fprintf(stdout,"  Boxsize  = %g \n", cp.lbox);
  fprintf(stdout,"  Redshift = %g \n", cp.redshift);
  fprintf(stdout,"  Omega Matter = %g \n", cp.omegam);
  fprintf(stdout,"  Omega Lambda = %g \n", cp.omegal);
  fprintf(stdout,"  Hubble parameter = %g \n",cp.hparam);
  fprintf(stdout,"  Particle Mass = %g \n",cp.Mpart);
  fprintf(stdout,"  Softening = %g\n",cp.soft);
  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"*********************************** \n");
  fflush(stdout);

  return;
}

static void lee(const char *filename, TYPE_INT * __restrict__ ind)
{
  FILE *pf;
  TYPE_INT d1, d2;
  TYPE_INT k, pc; 
	int n;
#ifdef TYPE_TWO_GADGET
  TYPE_INT blocksize;
  char label[4];
#endif

  TYPE_REAL r[3];
  #ifdef STORE_VELOCITIES
  TYPE_REAL v[3];
  #endif
  #ifdef STORE_IDS
  TYPE_INT id;
  #endif

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label, sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif
 
  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  assert(d1==256);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], sizeof(TYPE_REAL), 3, pf);
      if(k == TYPE_PART){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN
        P.x[pc] = r[0]*POSFACTOR;
        P.y[pc] = r[1]*POSFACTOR;
        P.z[pc] = r[2]*POSFACTOR;

#ifdef CHANGE_POSITION
        if(P.x[pc] > pmax[0]) pmax[0] = P.x[pc];
        if(P.x[pc] < pmin[0]) pmin[0] = P.x[pc];
        if(P.y[pc] > pmax[1]) pmax[1] = P.y[pc];
        if(P.y[pc] < pmin[1]) pmin[1] = P.y[pc];
        if(P.z[pc] > pmax[2]) pmax[2] = P.z[pc];
        if(P.z[pc] < pmin[2]) pmin[2] = P.z[pc];
#endif
      #else
        P[pc].pos[0] = r[0]*POSFACTOR;
        P[pc].pos[1] = r[1]*POSFACTOR;
        P[pc].pos[2] = r[2]*POSFACTOR;

#ifdef CHANGE_POSITION
        if(P[pc].pos[0] > pmax[0]) pmax[0] = P[pc].pos[0];
        if(P[pc].pos[0] < pmin[0]) pmin[0] = P[pc].pos[0];
        if(P[pc].pos[1] > pmax[1]) pmax[1] = P[pc].pos[1];
        if(P[pc].pos[1] < pmin[1]) pmin[1] = P[pc].pos[1];
        if(P[pc].pos[2] > pmax[2]) pmax[2] = P[pc].pos[2];
        if(P[pc].pos[2] < pmin[2]) pmin[2] = P[pc].pos[2];
#endif
      #endif

        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&v[0], sizeof(TYPE_REAL), 3, pf);
      if(k == TYPE_PART){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN        
        P.vx[pc] = v[0]*VELFACTOR;
        P.vy[pc] = v[1]*VELFACTOR;
        P.vz[pc] = v[2]*VELFACTOR;
      #else
        P[pc].vel[0] = v[0]*VELFACTOR;
        P[pc].vel[1] = v[1]*VELFACTOR;
        P[pc].vel[2] = v[2]*VELFACTOR;
      #endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_IDS
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, sizeof(TYPE_INT), 1, pf);
      if(k == TYPE_PART){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN
        P.id[pc] = id;
      #else
        P[pc].id = id;
      #endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind = pc;
  
  fclose(pf);
}

extern void read_gadget(void)
{
  char filename[200];
  int ifile;
  TYPE_INT ind;
  size_t total_memory;

#ifdef CHANGE_POSITION
  for(ind = 0; ind < 3; ind++)
  {
    pmin[ind] =  1.E26; 
    pmax[ind] = -1.E26;
  }
#endif

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  read_header(filename);

  /****** Allocate total particles ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  fprintf(stdout,"Allocating %.5zu Gb for %u particles\n", total_memory, cp.npart);
  fflush(stdout);
  if(!allocate_particles(&P, cp.npart))  exit(1);

  /****** Read POS, VEL and ID ********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++)
  {
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,&ind);
  }

  cp.lbox *= POSFACTOR;

  fprintf(stdout,"cp.lbox %f ...\n",cp.lbox);
  fprintf(stdout,"End reading snapshot file(s)...\n"); 
  fflush(stdout);
}

#ifdef CHANGE_POSITION

extern void change_positions(TYPE_INT n)
{
  TYPE_INT ip;

  fprintf(stdout,"xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  fprintf(stdout,"ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  fprintf(stdout,"zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip=0;ip<n;ip++)
  {
#ifdef COLUMN
    P.x[ip] -= pmin[0];
    P.y[ip] -= pmin[1];
    P.z[ip] -= pmin[2];
#else
    P[ip].pos[0] -= pmin[0];
    P[ip].pos[1] -= pmin[1];
    P[ip].pos[2] -= pmin[2];
#endif
  }

  cp.lbox = 0.0f;
  for(ip = 0; ip < 3; ip++)
    if(cp.lbox < (pmax[ip] - pmin[ip])) 
      cp.lbox = (pmax[ip] - pmin[ip]);

  cp.lbox *= 1.001;

  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  fflush(stdout);
}

#endif
