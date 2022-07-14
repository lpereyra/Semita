#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <math.h>

#include "colores.h"
#include "cosmoparam.h"
#include "variables.h"
#include "sph.h"
#include "grid.h"
#include "colores.h"
#include "leesnap.h"

static type_int  ***sum_fil;
static type_real ***sum_mean_par;
static type_real ***sum_mean_per;
static type_real ***sum_quad_per;
static type_int *id_fil;
static type_int *k_fil;

static struct info {
  type_int  id;
  type_real r_pos;
  type_real mod_vel;
};

static int compare_ids(const void *a, const void *b)
{
  int *c, *d;
  
  c = (int *) a;
  d = (int *) b;
  
  return (*c - *d);
}

static int compare_mod(const void *a, const void *b)
{

  if(((struct info *) a)->mod_vel < (((struct info *) b)->mod_vel))
    return -1;

  if(((struct info *) a)->mod_vel > (((struct info *) b)->mod_vel))
    return +1;

  return 0;
}

static int sgn(type_real a, type_real val) 
{
    return ((a > val) - (a < val));
}

static void search_cent(const type_int idpar, const type_int tid, const type_real h2) 
{	
  type_int  *unique_fil;
  long nfil,id_cent_seg;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox_grid[27], ibox;
  const type_real fac  = (type_real)grid.ngrid/cp.lbox;
  const type_real *Pos_cent = P[idpar].Pos;
  type_real r2, ss, mm, Posprima[3], vprima[3], vdir[3];
  type_real Posper[3]; 

  ixc  = (long)(Pos_cent[0]*fac);
  iyc  = (long)(Pos_cent[1]*fac);
  izc  = (long)(Pos_cent[2]*fac);

  ixci = ixc - 1;
  ixcf = ixc + 1;

  iyci = iyc - 1;
  iycf = iyc + 1;

  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= grid.ngrid ) ixcf = grid.ngrid - 1;
  if( iycf >= grid.ngrid ) iycf = grid.ngrid - 1;
  if( izcf >= grid.ngrid ) izcf = grid.ngrid - 1;
  #endif

  nfil = 0;
  ibox = 0;

  for(ixx = ixci; ixx <= ixcf; ixx++)
  {
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= (long)grid.ngrid) ix = ix - (long)grid.ngrid;
    if(ix < 0) ix = ix + grid.ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++)
    {
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= (long)grid.ngrid) iy = iy - (long)grid.ngrid;
      if(iy < 0) iy = iy + (long)grid.ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++)
      {
        iz = izz;
        #ifdef PERIODIC
        if(iz >= (long)grid.ngrid) iz = iz - (long)grid.ngrid;
        if(iz < 0) iz = iz + (long)grid.ngrid;
        #endif

        ibox_grid[ibox] = (ix * (long)grid.ngrid + iy) * (long)grid.ngrid + iz ;        
        nfil += grid.size[ibox_grid[ibox]];
        ibox++;
      }
    }  
  }

  unique_fil = (type_int *) malloc(nfil*sizeof(type_int));

  nfil = 0;
  for(ibox=0;ibox<27;ibox++)
  {
     id_cent_seg = grid.icell[ibox_grid[ibox]];

     while(id_cent_seg != grid.nobj)
     {
       unique_fil[nfil] = id_fil[id_cent_seg];
       id_cent_seg = grid.ll[id_cent_seg];
       nfil++;
     }
  }

  if(nfil>1)
  {
    qsort(unique_fil,nfil,sizeof(type_int),compare_ids);

    ix = 0;   
    for(iy=0; iy<nfil-1; iy++) 
      if(unique_fil[iy] != unique_fil[iy+1]) 
        unique_fil[ix++] = unique_fil[iy]; 
  
    unique_fil[ix++] = unique_fil[nfil-1];   
    nfil = ix;

    unique_fil = (type_int *) realloc(unique_fil,nfil*sizeof(type_int));
  }

  for(ibox=0;ibox<nfil;ibox++)
  {
    Posprima[0] = Seg[unique_fil[ibox]].Pos_list[0] - Pos_cent[0];
    Posprima[1] = Seg[unique_fil[ibox]].Pos_list[1] - Pos_cent[1];
    Posprima[2] = Seg[unique_fil[ibox]].Pos_list[2] - Pos_cent[2];

    #ifdef PERIODIC
    if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
    if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
    if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
    if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
    if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
    if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
    #endif

    r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

    if(r2 < Seg[unique_fil[ibox]].Rvir[0]) // uso el cuadrado   
      continue;

    Posprima[0] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size-1)]   - Pos_cent[0];
    Posprima[1] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size-1)+1] - Pos_cent[1];
    Posprima[2] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size-1)+2] - Pos_cent[2];

    #ifdef PERIODIC
    if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
    if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
    if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
    if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
    if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
    if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
    #endif

    r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

    if(r2 < Seg[unique_fil[ibox]].Rvir[1]) // uso el cuadrado   
      continue;

    /////////////////////////////////// REGION INTERNA ///////////////////////////

    mm = cp.lbox*cp.lbox;
    id_cent_seg = grid.nobj;

    for(ix=1;ix<Seg[unique_fil[ibox]].size;ix++)
    {
      Posprima[0] =  Pos_cent[0] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)]  ;
      Posprima[1] =  Pos_cent[1] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)+1];
      Posprima[2] =  Pos_cent[2] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)+2];

      #ifdef PERIODIC
      if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
      if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
      if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
      if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
      if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
      if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
      #endif

      vprima[1] = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2]; //r2

      if(vprima[1] <= h2) // entro en la esfera
      {

        vdir[0] = Seg[unique_fil[ibox]].Pos_list[3*ix]     - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)];
        vdir[1] = Seg[unique_fil[ibox]].Pos_list[3*ix + 1] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1) + 1];
        vdir[2] = Seg[unique_fil[ibox]].Pos_list[3*ix + 2] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1) + 2];
        #ifdef PERIODIC
        if(vdir[0]> 0.5*cp.lbox) vdir[0] -= cp.lbox;
        if(vdir[0]<-0.5*cp.lbox) vdir[0] += cp.lbox;
        if(vdir[1]> 0.5*cp.lbox) vdir[1] -= cp.lbox;
        if(vdir[1]<-0.5*cp.lbox) vdir[1] += cp.lbox;
        if(vdir[2]> 0.5*cp.lbox) vdir[2] -= cp.lbox;
        if(vdir[2]<-0.5*cp.lbox) vdir[2] += cp.lbox;
        #endif
        r2 = vdir[0]*vdir[0] + vdir[1]*vdir[1] + vdir[2]*vdir[2]; // para proyectar 
        r2 = 1.0/sqrt(r2);

        vdir[0] *= r2; // directo
        vdir[1] *= r2; // directo
        vdir[2] *= r2; // directo

        vprima[0] = Posprima[0]*vdir[0] + Posprima[1]*vdir[1] + Posprima[2] * vdir[2]; // dot

        vprima[2] = fabs(vprima[1] - vprima[0]*vprima[0]);

        if(vprima[2]<mm)
        {
          mm = vprima[2];
          id_cent_seg = ix-1; // me quedo con el anterior
        }
      }
    }

    if(id_cent_seg == grid.nobj || id_cent_seg == -1)
      continue;

    r2 = 0.0;
    for(ix=0;ix<3;ix++)
    {
      Posprima[ix] = Pos_cent[ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
      #ifdef PERIODIC
      if(Posprima[ix] >  0.5*cp.lbox) Posprima[ix] -= cp.lbox;
      if(Posprima[ix] < -0.5*cp.lbox) Posprima[ix] += cp.lbox;
      #endif

      vdir[ix] = Seg[unique_fil[ibox]].Pos_list[3*(id_cent_seg+1)+ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
      #ifdef PERIODIC
      if(vdir[ix] >  0.5*cp.lbox) vdir[ix] -= cp.lbox;
      if(vdir[ix] < -0.5*cp.lbox) vdir[ix] += cp.lbox;
      #endif

      r2 += (vdir[ix]*vdir[ix]);
    }

    r2 = 1.0/sqrt(r2);

    vdir[0] *= r2;
    vdir[1] *= r2;
    vdir[2] *= r2;
  
    mm = (Posprima[0]*vdir[0] + Posprima[1]*vdir[1] + Posprima[2]*vdir[2]);

    ss = r2 = 0.0f;
    for(ix=0;ix<3;ix++)
    {
      Posprima[ix] -= (mm*vdir[ix]);
    
      #ifdef PERIODIC
      if(Posprima[ix]>  0.5f*cp.lbox) Posprima[ix] -= cp.lbox;
      if(Posprima[ix]< -0.5f*cp.lbox) Posprima[ix] += cp.lbox;
      #endif

      r2 += (Posprima[ix]*Posprima[ix]);

      /////////////////////////////////////////////////////////

      vprima[ix] = P[idpar].Vel[ix]; // asigna

      //#ifdef CALCULA_VCM
        vprima[ix] -= Seg[unique_fil[ibox]].Vmedia[ix]; // asigna
      //#endif
      
      ss += (vprima[ix]*vdir[ix]);
    }
    r2 = 1.0/sqrt(r2);

    Posprima[0] *= r2;
    Posprima[1] *= r2;
    Posprima[2] *= r2;

    //mm = 0.0f;
    //for(ix=0;ix<3;ix++)
    //{
    //  vprima[ix] -= (ss*vdir[ix]);
    //  mm += (vprima[ix]*Posprima[ix]);
    //}

    mm = (vprima[0]*Posprima[0] + vprima[1]*Posprima[1] + vprima[2]*Posprima[2]);

    //Posper[0] = vdir[1]*Posprima[2]-vdir[2]*Posprima[1]; 
    //Posper[1] = vdir[2]*Posprima[0]-vdir[0]*Posprima[2];
    //Posper[2] = vdir[0]*Posprima[1]-vdir[1]*Posprima[0];

    //#ifdef PERIODIC
    //r2 = 0.0;
    //for(ix=0;ix<3;ix++)
    //{
    //  if(Posper[ix]>  0.5f*cp.lbox) Posper[ix] -= cp.lbox;
    //  if(Posper[ix]< -0.5f*cp.lbox) Posper[ix] += cp.lbox;
    //  r2 += Posper[ix]*Posper[ix];
    //}
    //r2 = 1.0/sqrt(r2);

    //Posper[0] *= r2;
    //Posper[1] *= r2;
    //Posper[2] *= r2;
    //#endif

    //r2 = (vprima[0]*Posper[0] + vprima[1]*Posper[1] + vprima[2]*Posper[2]);

    sum_mean_par[tid][unique_fil[ibox]][id_cent_seg] += ss;
    sum_mean_per[tid][unique_fil[ibox]][id_cent_seg] += mm;
    sum_quad_per[tid][unique_fil[ibox]][id_cent_seg] += (mm*mm); 
    //sum_mean_theta[tid][unique_fil[ibox]][id_cent_seg] +=  r2;
    //sum_quad_theta[tid][unique_fil[ibox]][id_cent_seg] += (r2*r2);  
    sum_fil[tid][unique_fil[ibox]][id_cent_seg]++;

  } //fin ibox

  free(unique_fil);
      
}

extern void calculo_rho(void)
{
  long i, j, k, size;
  type_real aux_real, r[3];
  int k_arvo;
  const type_int NTHREADS = omp_get_max_threads();

  aux_real  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.npart));
  aux_real *= cp.lbox;
  aux_real *= fof[1];
  size = 0;

  fprintf(stdout,"Init SPH\n");
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {
    if(i%10000==0)
    {
      fprintf(stdout,"Seg %u %.4f\n",i,(float)i/(float)cp.nseg);
      fflush(stdout);
    }

    Seg[i].Rvir[0] *= Seg[i].Rvir[0]; // Uso el cuadrado 
    Seg[i].Rvir[1] *= Seg[i].Rvir[1]; //
    size += Seg[i].size;
  }

  id_fil  = (type_int *) malloc(size*sizeof(type_int));
  k_fil   = (type_int *) malloc(size*sizeof(type_int));
  sum_fil = (type_int ***) calloc(NTHREADS,sizeof(type_int **));
  sum_mean_par = (type_real ***) calloc(NTHREADS,sizeof(type_real **));
  sum_mean_per = (type_real ***) calloc(NTHREADS,sizeof(type_real **));
  sum_quad_per = (type_real ***) calloc(NTHREADS,sizeof(type_real **));
  //sum_mean_theta = (type_real ***) calloc(NTHREADS,sizeof(type_real **));
  //sum_quad_theta = (type_real ***) calloc(NTHREADS,sizeof(type_real **));
  for(i=0;i<NTHREADS;i++)
  {
    sum_quad_per[i] = (type_real **) calloc(cp.nseg,sizeof(type_real *));
    sum_mean_per[i] = (type_real **) calloc(cp.nseg,sizeof(type_real *));
    //sum_quad_theta[i] = (type_real **) calloc(cp.nseg,sizeof(type_real *));
    //sum_mean_theta[i] = (type_real **) calloc(cp.nseg,sizeof(type_real *));
    sum_mean_par[i] = (type_real **) calloc(cp.nseg,sizeof(type_real *));
    sum_fil[i]      = (type_int **)  calloc(cp.nseg,sizeof(type_int  *));
    for(k=0;k<cp.nseg;k++)
    {
      sum_quad_per[i][k] = (type_real *) calloc(Seg[k].size,sizeof(type_real));
      sum_mean_per[i][k] = (type_real *) calloc(Seg[k].size,sizeof(type_real));
      //sum_quad_theta[i][k] = (type_real *) calloc(Seg[k].size,sizeof(type_real));
      //sum_mean_theta[i][k] = (type_real *) calloc(Seg[k].size,sizeof(type_real));
      sum_mean_par[i][k] = (type_real *) calloc(Seg[k].size,sizeof(type_real));
      sum_fil[i][k]      = (type_int *)  calloc(Seg[k].size,sizeof(type_int ));
    }
  }

  grid.nobj = size;
  grid.ngrid = (long)(cp.lbox/(sqrt(2.0)*R_SPH));

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build(id_fil, k_fil);

  read_gadget();

#ifdef CALCULA_MEDIA

  vmedia   = (type_real **) calloc(NTHREADS,sizeof(type_real *));
  sum_part = (type_int  **) calloc(NTHREADS,sizeof(type_real *));
  for(i=0;i<NTHREADS;i++)
  {
    vmedia[i]   = (type_real *) calloc(3*cp.nseg,sizeof(type_real *));
    sum_part[i] = (type_real *) calloc(cp.nseg,sizeof(type_real *));
  }

  #pragma omp parallel for schedule(dynamic) num_threads(NTHREADS) \
  default(none) private(i,k) shared(cp,sum_fil,sum_mean_par,\
  sum_mean_per,sum_quad_per,P,Seg,stdout) 
  //sum_mean_per,sum_quad_per,sum_mean_theta,sum_quad_theta,P,Seg,stdout) 
  for(i=0;i<cp.npart;i++)
  {
    if(i%10000000==0)
    {
      fprintf(stdout,"Mean Particles %u %.4f\r",i,(float)i/(float)cp.npart);
      fflush(stdout);
    }

    k = omp_get_thread_num();
    calc_fil_mean(i,k,R_SPH*R_SPH);
  }

  for(i=0;i<cp.nseg;i++)  
  {
    size = 0;
    for(k=0;k<NTHREADS;k++)
    {
      Seg[i].Vmedia[0] += vmedia[k][3*i]; 
      Seg[i].Vmedia[1] += vmedia[k][3*i+1];
      Seg[i].Vmedia[2] += vmedia[k][3*i+2];
      size += sum_part[k][i];
    }

    Seg[i].Vmedia[0] *= (size>0 ? 1.0f/(type_real)size : 0.0);
    Seg[i].Vmedia[1] *= (size>0 ? 1.0f/(type_real)size : 0.0);
    Seg[i].Vmedia[2] *= (size>0 ? 1.0f/(type_real)size : 0.0);
  }

  for(i=0;i<NTHREADS;i++)
  {
    free(vmedia[i]);
    free(sum_part[i]);
  }

  free(sum_part);
  free(vmedia);

#endif

  #pragma omp parallel for schedule(dynamic) num_threads(NTHREADS) \
  default(none) private(i,k) shared(cp,sum_fil,sum_mean_par,\
  sum_mean_per,sum_quad_per,P,Seg,stdout) 
  //sum_mean_per,sum_quad_per,sum_mean_theta,sum_quad_theta,P,Seg,stdout) 
  for(i=0;i<cp.npart;i++)
  {
    if(i%10000000==0)
    {
      fprintf(stdout,"Particles %u %.4f\r",i,(float)i/(float)cp.npart);
      fflush(stdout);
    }

    k = omp_get_thread_num();
    search_cent(i,k,R_SPH*R_SPH);
  }

  #pragma omp parallel for num_threads(NTHREADS) \
  default(none) private(i,j,k,k_arvo,size,aux_real,r) \
  shared(cp,Seg,sum_fil,sum_mean_par,sum_mean_per,sum_quad_per,\
  stdout)
  //sum_mean_theta,sum_quad_theta,stdout)
  for(i=0;i<cp.nseg;i++)
  {
    Seg[i].Rvir[0] = sqrt(Seg[i].Rvir[0]);
    Seg[i].Rvir[1] = sqrt(Seg[i].Rvir[1]);
    
    for(k=1;k<NTHREADS;k++)
    {
      for(j=0;j<Seg[i].size;j++)
      {
        //sum_quad_theta[0][i][j] += sum_quad_theta[k][i][j];
        //sum_mean_theta[0][i][j] += sum_mean_theta[k][i][j];
        sum_quad_per[0][i][j]   += sum_quad_per[k][i][j];
        sum_mean_per[0][i][j]   += sum_mean_per[k][i][j];
        sum_mean_par[0][i][j]   += sum_mean_par[k][i][j];
             sum_fil[0][i][j]   +=      sum_fil[k][i][j];
      }
    }

    size = 0;
    Seg[i].sigma_per    = 0.0f;
    //Seg[i].sigma_theta  = 0.0f;

#ifdef POR_CILINDRO

    for(j=0;j<Seg[i].size;j++)
    {

      size                  += sum_fil[0][i][j];

      sum_mean_per[0][i][j]   *= (sum_fil[0][i][j]>0 ? 1./(type_real)sum_fil[0][i][j] : 0.0f);
      sum_mean_par[0][i][j]   *= (sum_fil[0][i][j]>0 ? 1./(type_real)sum_fil[0][i][j] : 0.0f);
      //sum_mean_theta[0][i][j] *= (sum_fil[0][i][j]>0 ? 1./(type_real)sum_fil[0][i][j] : 0.0f);

      sum_quad_per[0][i][j]    = (sum_fil[0][i][j]>0 ? sum_quad_per[0][i][j] : 0.0f);
      //sum_quad_theta[0][i][j]  = (sum_fil[0][i][j]>0 ? sum_quad_theta[0][i][j] : 0.0f);

      Seg[i].sigma_per += sum_fil[0][i][j]>1 ? fabs(sum_quad_per[0][i][j] - \
      (type_real)sum_fil[0][i][j]*sum_mean_per[0][i][j]*sum_mean_per[0][i][j])* \  
      (1.f/(type_real)(sum_fil[0][i][j]-1)) : 0.0f;

      //Seg[i].sigma_theta += sum_fil[0][i][j]>1 ? fabs(sum_quad_theta[0][i][j] - \
      //(type_real)sum_fil[0][i][j]*sum_mean_theta[0][i][j]*sum_mean_theta[0][i][j])* \
      //(1.f/(type_real)(sum_fil[0][i][j]-1)) : 0.0f;

    }

    Seg[i].mass_part  = cp.Mpart*(type_real)size;
    Seg[i].sigma_per   = sqrt(Seg[i].sigma_per   / (type_real)(Seg[i].size));
    //Seg[i].sigma_theta = sqrt(Seg[i].sigma_theta / (type_real)(Seg[i].size));

    aux_real = Seg[i].len - Seg[i].Rvir[0] - Seg[i].Rvir[1];
    aux_real *= 1e-3;

    Seg[i].mu   = aux_real<0.0 ? 0.0f : Seg[i].mass_part/aux_real;
    Seg[i].rho  = aux_real<0.0 ? 0.0f : Seg[i].mass_part/Seg[i].vol;

#else

    type_real vel_per   = 0.0f;
    type_real vel_theta = 0.0f;

    for(j=0;j<Seg[i].size;j++)
    {
      size                  += sum_fil[0][i][j];
      vel_per               += sum_mean_per[0][i][j];
      //vel_theta             += sum_mean_theta[0][i][j];
      Seg[i].sigma_per      += sum_quad_per[0][i][j];
      //Seg[i].sigma_theta    += sum_quad_theta[0][i][j];
    }

    aux_real = Seg[i].len - Seg[i].Rvir[0] - Seg[i].Rvir[1];
    aux_real *= 1e-3;

    Seg[i].mass_part  = cp.Mpart*(type_real)size;
    Seg[i].rho = aux_real<0.0 ? 0.0f : Seg[i].mass_part/Seg[i].vol;
    Seg[i].mu   = aux_real<0.0 ? 0.0f : Seg[i].mass_part/aux_real;
    
    vel_per   *= ((size>1) ? (1./(type_real)size) : 0.0f);
    vel_theta *= ((size>1) ? (1./(type_real)size) : 0.0f);
    
    Seg[i].sigma_per = size>1 ? sqrt(fabs(Seg[i].sigma_per - \
    (type_real)size*vel_per*vel_per)*(1.f/(type_real)(size-1))) : 0.0f;

    //Seg[i].sigma_theta = size>1 ? sqrt(fabs(Seg[i].sigma_theta - \
    //(type_real)size*vel_theta*vel_theta)*(1.f/(type_real)(size-1))) : 0.0f;

#endif

    //Seg[i].Rvir[0] /= RVIR_FACTOR;
    //Seg[i].Rvir[1] /= RVIR_FACTOR;

    struct info *vel_silla;

    vel_silla = (struct info *) calloc(Seg[i].size,sizeof(struct info));

    size = 0;
    aux_real = 0.0f;
    k = -1; // para no repetir

    for(j=1;j<Seg[i].size;j++)
    {
      r[0] = Seg[i].Pos_list[3*j]   - Seg[i].Pos_list[3*(j-1)];
      r[1] = Seg[i].Pos_list[3*j+1] - Seg[i].Pos_list[3*(j-1)+1];
      r[2] = Seg[i].Pos_list[3*j+2] - Seg[i].Pos_list[3*(j-1)+2];

      #ifdef PERIODIC
      if(r[0]> 0.5*cp.lbox) r[0] -= cp.lbox;
      if(r[1]> 0.5*cp.lbox) r[1] -= cp.lbox;
      if(r[2]> 0.5*cp.lbox) r[2] -= cp.lbox;

      if(r[0]<-0.5*cp.lbox) r[0] += cp.lbox;
      if(r[1]<-0.5*cp.lbox) r[1] += cp.lbox;
      if(r[2]<-0.5*cp.lbox) r[2] += cp.lbox;
      #endif

      k_arvo = (sgn(sum_mean_par[0][i][j-1], 0.0f) == sgn(sum_mean_par[0][i][j], 0.0f)) ? \
      1 : 0; // cambio de signo

      if(sum_fil[0][i][j-1]>=10 && k_arvo==0 && (k != (j-1)))
      {
        vel_silla[size].id      = j-1;
        vel_silla[size].r_pos   = aux_real; 
        vel_silla[size].mod_vel = \
        sum_mean_par[0][i][j-1]*sum_mean_par[0][i][j-1]; // el cuadrado
        //sum_mean_per[0][i][j]*sum_mean_per[0][i][j]+\
        //sum_mean_par[0][i][j]*sum_mean_par[0][i][j]; // el cuadrado del mod vel
        size++;
      }

      aux_real += sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

      if(sum_fil[0][i][j]>=10 && k_arvo==0)
      {
        k = j;
        vel_silla[size].id      = j;
        vel_silla[size].r_pos   = aux_real; 
        vel_silla[size].mod_vel = \
        sum_mean_par[0][i][j]*sum_mean_par[0][i][j]; // el cuadrado
        //sum_mean_per[0][i][j]*sum_mean_per[0][i][j]+\
        //sum_mean_par[0][i][j]*sum_mean_par[0][i][j]; // el cuadrado del mod vel
        size++;
      }
    }

    if(size>0)
    {
      qsort(vel_silla,size,sizeof(struct info),compare_mod);

      Seg[i].r_pos = (vel_silla[0].id == (Seg[i].size-1)) ? 0.0 : vel_silla[0].r_pos/Seg[i].len;
      Seg[i].id_pos = \
      ((vel_silla[0].id == 0) || (vel_silla[0].id == (Seg[i].size-1))) ? -1 : vel_silla[0].id;
    }else{
      Seg[i].r_pos = 0.0;
      Seg[i].id_pos = -1;
    }

    #ifdef PERIODIC
    for(k=0;k<Seg[i].size;k++)
    {
      r[0] = Seg[i].Pos_list[3*k+0]+Seg[i].Pos_list[0];
      r[1] = Seg[i].Pos_list[3*k+1]+Seg[i].Pos_list[1];
      r[2] = Seg[i].Pos_list[3*k+2]+Seg[i].Pos_list[2];

      if(r[0]> 0.5*cp.lbox) Seg[i].Pos_list[3*k]   -= cp.lbox;
      if(r[1]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+1] -= cp.lbox;
      if(r[2]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+2] -= cp.lbox;

      if(r[0]<-0.5*cp.lbox) Seg[i].Pos_list[3*k]   += cp.lbox;
      if(r[1]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+1] += cp.lbox;
      if(r[2]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+2] += cp.lbox;

      Seg[i].Pos_list[3*k] = Seg[i].Pos_list[3*k] >= cp.lbox ? \
      Seg[i].Pos_list[3*k]-cp.lbox : Seg[i].Pos_list[3*k];
      Seg[i].Pos_list[3*k+1] = Seg[i].Pos_list[3*k+1] >= cp.lbox ? \
      Seg[i].Pos_list[3*k+1]-cp.lbox : Seg[i].Pos_list[3*k+1];
      Seg[i].Pos_list[3*k+2] = Seg[i].Pos_list[3*k+2] >= cp.lbox ? \
      Seg[i].Pos_list[3*k+2]-cp.lbox : Seg[i].Pos_list[3*k+2];

      Seg[i].Pos_list[3*k] = Seg[i].Pos_list[3*k] < 0.0f ? \
      Seg[i].Pos_list[3*k]+cp.lbox : Seg[i].Pos_list[3*k];
      Seg[i].Pos_list[3*k+1] = Seg[i].Pos_list[3*k+1] < 0.0f ? \
      Seg[i].Pos_list[3*k+1]+cp.lbox : Seg[i].Pos_list[3*k+1];
      Seg[i].Pos_list[3*k+2] = Seg[i].Pos_list[3*k+2] < 0.0f ? \
      Seg[i].Pos_list[3*k+2]+cp.lbox : Seg[i].Pos_list[3*k+2];
    }
    #endif

    free(vel_silla);
  }

  fprintf(stdout,"End SPH\n");
  fflush(stdout);

  for(i=0;i<NTHREADS;i++)
  {
    for(j=0;j<cp.nseg;j++)
    {
      free(sum_quad_per[i][j]);
      free(sum_mean_per[i][j]);
      //free(sum_quad_theta[i][j]);
      //free(sum_mean_theta[i][j]);
      free(sum_mean_par[i][j]);
      free(sum_fil[i][j]);
    }
    free(sum_fil[i]);
    free(sum_mean_par[i]);
    free(sum_mean_per[i]);
    free(sum_quad_per[i]);
    //free(sum_mean_theta[i]);
    //free(sum_quad_theta[i]);
  }
  free(sum_fil);
  free(sum_mean_par);
  free(sum_mean_per);
  free(sum_quad_per);
  //free(sum_mean_theta);
  //free(sum_quad_theta);

  free(id_fil);
  free(k_fil);
  grid_free();
  free(P);

  return;
}
