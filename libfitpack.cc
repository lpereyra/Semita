#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "variables.hh"
#include "colors.hh"
#include "libfitpack.hh"
#include "vol.hh"

#define IX(i,j) 3*i+j

static void smooth(int n, int idim, float *Pos, double s, int npoint)
{
  int k = 3, ipar = 0, iopt = 0, m = 0;  
  int i,j,ier, nest, nx, nc, lwrk;
  int *iwrk;
  double ub, ue, fp;
  double *c, *t, *wrk, *wa;
  double *x, *w, *u;

  assert(n > k);
  assert(k>=1 && k<=5);

  nest = n+2*k>2*k+3 ? n+2*k : 2*k+3;

  nx = idim*n;
  nc = idim*nest;
  lwrk = n*(k + 1) + nest*(6 + idim + 3*k);
  j = nc + 2*nest + lwrk;

  wa = (double *) malloc(j*sizeof(double));
  t = wa;
  c = t + nest;
  wrk = c + nc;
  iwrk = (int *)(wrk + lwrk);

  x = (double *) malloc(nx*sizeof(double));
  w = (double *) malloc(n*sizeof(double));
  u = (double *) malloc(n*sizeof(double));
  
  for(j=0;j<n;j++)
  {
    for(i=0;i<idim;i++)
    {
      x[idim*j+i] = (double)Pos[idim*j+i];
    }
    w[j] = 1.0f;
    u[j] = 0.0f;
  }

  parcur_(&iopt, &ipar, &idim, &n, u, &nx, x, w, &ub, &ue, &k, \
          &s, &nest, &m, t, &nc, c, &fp, wrk, &lwrk, iwrk, &ier);

  if(ier == 10) 
  {
    assert(0);
  }

  if(ier > 0 && m == 0) {
     n = 1;
  }

  free(u);
  free(x);
  free(w);
  
  nx = idim*npoint;
  x   = (double *) malloc(npoint*sizeof(double));
  w   = (double *) malloc(npoint*sizeof(double));

  fp = 1./(npoint-1);
  for(j=0;j<npoint;j++)
    x[j] = j*fp;

  for(i=0;i<idim;i++)
  {
    splev_(t,&m,c+i*m,&k,x,w,&npoint,&iopt,&ier);

    for(j=0;j<npoint;j++)
      Pos[j*idim+i] = (float)w[j];

    if(ier == 10) 
      assert(0);
  }

  free(x);
  free(w);
  free(wa);

  return;
}

static TYPE_REAL length(TYPE_REAL *Pos_list, TYPE_INT nsize)
{
	TYPE_REAL r, tmp[3], len;
	TYPE_INT k, idim;

	len = 0.0f;
  for(k=1;k<nsize;k++)
  {
    r = 0.0f;
    for(idim=0;idim<3;idim++)
    { 
      tmp[idim] = Pos_list[IX(k,idim)]-Pos_list[IX((k-1),idim)];
 
      #ifdef PERIODIC
      if(tmp[idim] >  0.5*cp.lbox) tmp[idim] -= cp.lbox;
      if(tmp[idim] < -0.5*cp.lbox) tmp[idim] += cp.lbox;
      #endif

      r += tmp[idim]*tmp[idim];
    }
    len += sqrt(r);
  }

	return len;
}

static TYPE_REAL curliness(TYPE_REAL *Pos_list, TYPE_INT nsize, TYPE_REAL len)
{
	TYPE_REAL r, tmp[3];
	TYPE_INT  idim;
	
	r = 0.0f;
	for(idim=0;idim<3;idim++)
	{ 
	  tmp[idim] = Pos_list[IX((nsize-1),idim)]-Pos_list[IX(0,idim)];
	
	  #ifdef PERIODIC
	  if(tmp[idim] >  0.5*cp.lbox) tmp[idim] -= cp.lbox;
	  if(tmp[idim] < -0.5*cp.lbox) tmp[idim] += cp.lbox;
	  #endif
	
	  r += tmp[idim]*tmp[idim];
	}
	
	r = sqrt(r);
	
	return r/len;
}

static TYPE_REAL rms(TYPE_REAL *Pos_list, TYPE_INT nsize)
{
	TYPE_REAL r, tmp[6], rms;
	TYPE_INT  k, idim;
	
	r = 0.0f;
	for(idim=0;idim<3;idim++)
	{ 
	  tmp[idim+3] = Pos_list[IX((nsize-1),idim)]-Pos_list[IX(0,idim)];
	
	  #ifdef PERIODIC
	  if(tmp[idim+3] >  0.5*cp.lbox) tmp[idim+3] -= cp.lbox;
	  if(tmp[idim+3] < -0.5*cp.lbox) tmp[idim+3] += cp.lbox;
	  #endif
	
	  r += tmp[idim+3]*tmp[idim+3];
	}
	r = sqrt(r);

	for(idim=0;idim<3;idim++)
	  tmp[idim+3] /= r;

  rms = 0.0f;
  for(k=1;k<nsize-1;k++)
  {
    r = 0.0f;
    for(idim=0;idim<3;idim++)
    { 
      tmp[idim] = Pos_list[IX(k,idim)]-Pos_list[IX(0,idim)];
 
      #ifdef PERIODIC
      if(tmp[idim] >  0.5*cp.lbox) tmp[idim] -= cp.lbox;
      if(tmp[idim] < -0.5*cp.lbox) tmp[idim] += cp.lbox;
      #endif
    }

    r = tmp[0]*tmp[3]+tmp[1]*tmp[4]+tmp[2]*tmp[5];              // dot
    r = tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2] - r*r;  // dist - dot
    r = fabs(r);  // por errores de redondeo
    rms += r;
  }

  rms /= (TYPE_REAL)nsize;
  rms = sqrt(rms);

	return rms;
}

extern void suaviza()
{
  TYPE_INT i,k,idim;
#ifdef PERIODIC
  TYPE_REAL diff;
#endif
#ifdef ITERA
  TYPE_INT j;
  TYPE_REAL *Pos_aux;
#endif
  TYPE_REAL fix_array[6];
  //const TYPE_INT NTHREADS = omp_get_max_threads();

  cp.lbox /= 1000.;

  for(i=0;i<cp.nseg;i++)
    for(k=0;k<Seg[i].size;k++)
      for(idim=0;idim<3;idim++)
        Seg[i].Pos_list[IX(k,idim)] /= 1000.;

  #ifdef ITERA

  GREEN("*********** ITERA ************\n");
  GREEN("******************************\n");
  fflush(stdout); 

  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].size==2) continue;
      
    Pos_aux = (TYPE_REAL *) malloc(3*Seg[i].size*sizeof(TYPE_REAL));

    for(j=0;j<NITERA;j++)
    {
      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          Pos_aux[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)];

          #ifdef PERIODIC
          diff = Seg[i].Pos_list[IX(k,idim)] - Seg[i].Pos_list[IX(0,idim)];
          if(diff> 0.5*cp.lbox) Pos_aux[IX(k,idim)] -= cp.lbox;
          if(diff<-0.5*cp.lbox) Pos_aux[IX(k,idim)] += cp.lbox;
          #endif
        }
      }

      for(k=1;k<Seg[i].size-1;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          //Seg[i].Pos_list[IX(k,idim)] = \
          //Pos_aux[IX(k,idim)];

          Seg[i].Pos_list[IX(k,idim)] = \
          0.25f*Pos_aux[IX((k-1),idim)]  + \
          0.50f*Pos_aux[IX(k,idim)]      + \
          0.25f*Pos_aux[IX((k+1),idim)];
  
          #ifdef PERIODIC
          while(Seg[i].Pos_list[IX(k,idim)] >= cp.lbox) 
            Seg[i].Pos_list[IX(k,idim)] -= cp.lbox;

          while(Seg[i].Pos_list[IX(k,idim)] < 0.0f) 
            Seg[i].Pos_list[IX(k,idim)] += cp.lbox;
          #endif
        }       
      }

    }

    free(Pos_aux);
  }

  GREEN("********* END ITERA **********\n");
  GREEN("******************************\n");
  fflush(stdout); 
  
  #endif

  GREEN("********** SUAVIZA ***********\n");
  GREEN("******************************\n");
  fflush(stdout); 

  for(i=0;i<cp.nseg;i++)
  {
    #ifdef FIX_NSMOOTH
      TYPE_INT Nsm = NSMOOTH;
    #else
      TYPE_INT Nsm = ceil(Seg[i].len/(TYPE_REAL)RSPACE);
    #endif

    if(Seg[i].size>=4)
    {

      if(Nsm<=Seg[i].size) Nsm = Seg[i].size;

      Seg[i].Pos_list = (TYPE_REAL *) realloc(Seg[i].Pos_list,3*Nsm*sizeof(TYPE_REAL));

      #ifdef PERIODIC
      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          diff = Seg[i].Pos_list[IX(k,idim)] - Seg[i].Pos_list[IX(0,idim)];
          if(diff> 0.5*cp.lbox) Seg[i].Pos_list[IX(k,idim)] -= cp.lbox;
          if(diff<-0.5*cp.lbox) Seg[i].Pos_list[IX(k,idim)] += cp.lbox;
        }
      }
      #endif
      
			for(idim=0;idim<3;idim++)
      {
        fix_array[idim]   = Seg[i].Pos_list[IX(0,idim)];
        fix_array[idim+3] = Seg[i].Pos_list[IX((Seg[i].size-1),idim)];
      }

      smooth((int)Seg[i].size,3,Seg[i].Pos_list,2.0,(int)Nsm);

      Seg[i].size = Nsm;

      for(idim=0;idim<3;idim++)
      {
        Seg[i].Pos_list[IX(0,idim)]               = fix_array[idim];
        Seg[i].Pos_list[IX((Seg[i].size-1),idim)] = fix_array[idim+3];
      }
     
      #ifdef PERIODIC
      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          Seg[i].Pos_list[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)] >= cp.lbox ? \
            Seg[i].Pos_list[IX(k,idim)]-cp.lbox : Seg[i].Pos_list[IX(k,idim)];

          Seg[i].Pos_list[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)] < 0.0f ? \
            Seg[i].Pos_list[IX(k,idim)]+cp.lbox : Seg[i].Pos_list[IX(k,idim)];
        }
      }
      #endif
    }

    for(k=0;k<Seg[i].size;k++)
      for(idim=0;idim<3;idim++)
        Seg[i].Pos_list[IX(k,idim)] *= 1000;
  }

  cp.lbox *= 1000.;

  GREEN("******** END SUAVIZADO ********\n");
  GREEN("*******************************\n");
  fflush(stdout); 

  GREEN("***** CALCULA PROPIEDADES *****\n");
  GREEN("*******************************\n");
  fflush(stdout); 

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  private(i,k,idim,fix_array) shared(cp,Seg,stdout)
  for(i=0;i<cp.nseg;i++)
  {
    /////////////////////////////////////////////////////////////////////  

    Seg[i].len = length(Seg[i].Pos_list,Seg[i].size);
		Seg[i].cur = curliness(Seg[i].Pos_list,Seg[i].size,Seg[i].len);
		Seg[i].rms = rms(Seg[i].Pos_list,Seg[i].size);
    
		///////////////////////////////////////////////////////////////////////  
    
		Seg[i].len_raw = length(Seg[i].Pos_list_raw,Seg[i].size_raw);
		Seg[i].cur_raw = curliness(Seg[i].Pos_list_raw,Seg[i].size_raw,Seg[i].len_raw);
		Seg[i].rms_raw = rms(Seg[i].Pos_list_raw,Seg[i].size_raw);

    ///////////////////////////////////////////////////////////////////// 

#ifdef WITH_EXTREMES
	 Seg[i].vol = volume(RADIUS_TUBE, Seg[i].size, Seg[i].Pos_list, Seg[i].Rvir);
#else
	 Seg[i].vol = volume(RADIUS_TUBE, Seg[i].size, Seg[i].Pos_list);
#endif

  }

  return;
}
