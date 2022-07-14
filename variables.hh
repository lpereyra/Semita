/* file variables.h
 * brief declares global variables.
 *
 * This file declares all global variables. 
 * Further variables should be added here, and declared as 'extern'. 
 * The actual existence of these variables is provided by the file 'variables.c'. 
 * To produce 'variables.c' from 'variables.h', do the following:
 *
 *    - Erase all #define's, typedef's, and enum's
 *    - add #include "variables.h", delete the #ifndef VARIABLES_H conditional
 *    - delete all keywords 'extern'
 *    - delete all struct definitions enclosed in {...}, e.g.
 *      "extern struct cosmoparam {....} cp;"
 *      becomes "struct cosmoparam cp;"
 */

#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NPARTMIN
  #define NPARTMIN 20
#endif

#ifndef TYPE_PART
  #define TYPE_PART 1
#endif

#ifndef LEVEL_PRUNED
  #define LEVEL_PRUNED 4
#endif

#ifdef ITERA
  #ifndef NITERA
    #define NITERA 1
  #endif
#endif

#ifdef FIX_NSMOOTH
  #ifndef NSMOOTH
    #define NSMOOTH 30
  #endif
#else
  #ifndef RSPACE
    #define RSPACE 500.0
  #endif
#endif

#ifndef RADIUS_TUBE
  #define RADIUS_TUBE 2000.0
#endif

#define N_part_types 6    /* Number of particle types */

/* If defined, the type variable 
 * is set to "double", otherwise to "float" 
 */
#ifdef PRECDOUBLE
typedef double TYPE_REAL;
#else
typedef float TYPE_REAL;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long long TYPE_INT;
#else
typedef unsigned int TYPE_INT;
#endif

#define RHOCRIT 2.77525E11   /* Critical density of the Universe [Msol h² / Mpc³] */
#define GCONS 6.67300E-20    /* Gravitational Constant [km³ / kg / seg²]          */
#define Msol 1.9891E30       /* Sun Mass [Kg]                                     */
#define Kpc 3.08568025E16    /* Kiloparsec -> Kilometers                          */  

extern struct cosmoparam
{
  double   omegam    ;  /* Omega Matter                             */
  double   omegal    ;  /* Omega Lambda                             */
  double   omegak    ;  /* Omega Curvature                          */
  double   hparam    ;  /* Hubble constant in units of 100 km/s/Mpc */
  double   lbox      ;  /* Boxsize [Kpc / h]                        */
  double   Mpart     ;  /* Particle Mass [10^10 Msol / h]           */
  TYPE_INT npart_tot ;  /* Particle total number                    */
  TYPE_INT npart     ;  /* Particle number                          */
  TYPE_INT ngrup     ;  /* Groups   number                          */
  TYPE_INT nseg      ;  /* Segments number                          */
  double   redshift  ;  /* Redshift                                 */
  double   aexp      ;  /*                                          */
  double   Hubble_a  ;  /*                                          */
  double   soft      ;  /* Softening [Kpc / h]                      */
} cp;

/* Input and output files */
extern struct SnapST
{
  int nfiles;
  char root[200], name[200];
  int num;
} snap;

extern struct gridst
{
  long ngrid;
  unsigned long nobj;
  TYPE_INT *icell;
} grid;

/* This structure holds all the information that is
 * stored for each particle of the simulation.
 */
#ifdef COLUMN

  extern struct particle_data 
  {
    TYPE_REAL      *x;
    TYPE_REAL      *y;
    TYPE_REAL      *z;
    #ifdef STORE_VELOCITIES
    TYPE_REAL      *vx;
    TYPE_REAL      *vy;
    TYPE_REAL      *vz;
    #endif
    #ifdef STORE_IDS
    TYPE_INT       *id;
    #endif
    TYPE_INT       *sub;
  } P;

#else

  extern struct particle_data 
  {
    TYPE_REAL pos[3];
    #ifdef STORE_VELOCITIES
    TYPE_REAL vel[3];
    #endif
    #ifdef STORE_IDS
    TYPE_INT       id;
    #endif
    TYPE_INT       sub;
  } *P;

#endif

extern struct grup_data
{
  TYPE_INT   save;
  TYPE_INT   id;
  TYPE_INT   npart;
  TYPE_REAL  *pos;
  TYPE_REAL  pcm[3];
  TYPE_REAL  aa, bb, cc;
  TYPE_REAL  evec[3][3];
#ifdef STORE_VELOCITIES  
  TYPE_REAL  *vel;
  TYPE_REAL  vcm[3], sig[3];
  TYPE_REAL  r200, m200, v200;
  TYPE_REAL  rvir, mvir, vvir;
  TYPE_REAL  vmax;
  TYPE_REAL  L[3]; 
  TYPE_REAL  aa_vel, bb_vel, cc_vel;
  TYPE_REAL  evec_vel[3][3];
#ifdef COMPUTE_EP
  TYPE_REAL  mostbound[3];
  TYPE_REAL  Ep, Ec, lambda;
#endif
#endif
} *Gr;

extern struct segmentstd
{
  TYPE_INT  size;
  TYPE_INT  *Ids_list;
  TYPE_REAL *Pos_list;
  TYPE_REAL *Pos_list_raw;
  TYPE_INT  flag;
  TYPE_REAL Mass[2];
  TYPE_REAL razon;
  TYPE_REAL len;
  TYPE_REAL cur;
  TYPE_REAL rms;
  TYPE_INT  size_raw;
  TYPE_REAL len_raw;
  TYPE_REAL cur_raw;
  TYPE_REAL rms_raw;
  TYPE_REAL vol;
  TYPE_REAL rho;
  TYPE_REAL mu;
  TYPE_REAL mass_part;
#ifdef STORE_VELOCITIES
  TYPE_REAL Vnodos[6];
  TYPE_REAL sigma_per;
  TYPE_REAL r_pos;
  int id_pos;
#endif
#ifdef WITH_EXTREMES
  TYPE_REAL Rvir[2];
#endif
} *Seg;

extern TYPE_INT  nfrac;
extern TYPE_REAL *fof;
extern char message[200];
#ifdef CHANGE_POSITION
  extern TYPE_REAL pmin[3], pmax[3];
#endif

#endif
