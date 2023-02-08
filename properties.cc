#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "variables.hh"
#include "properties.hh"
#include "octree.hh"
#include "deltas.hh"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

/* ============= Unidades =============== 
   Energy = 10^10 Msol (km/seg)^2 
	 Mom. Angular = 10^10 Msol kpc km/seg
	 M(200,vir) = 10^10 Msol/h 
	 R(200,vir) = kpc/h 
	 V(200,vir) = km/seg
	 ====================================== */

static void forma(struct grup_data *Prop)
{
  TYPE_INT i, j, k, kk;
  double data_pos[9];
#ifdef STORE_VELOCITIES  
  double data_vel[9];
#endif
  gsl_matrix *evec;
  gsl_vector *eval;
  gsl_eigen_symmv_workspace * w;
  gsl_matrix_view m;

  evec = gsl_matrix_alloc (3, 3);
  eval = gsl_vector_alloc (3);

  Prop->aa = 0.0;
  Prop->cc = 0.0;
  Prop->bb = 0.0;
  
#ifdef STORE_VELOCITIES  
  Prop->aa_vel = 0.0;
  Prop->cc_vel = 0.0;
  Prop->bb_vel = 0.0;
#endif

  /* Calc. 'algo asi como' el Tensor de Inercia */
  kk=0;
  for(i=0 ; i<3 ; i++)
    for(j=0 ; j<3 ; j++)
    {
#ifdef STORE_VELOCITIES
      data_pos[kk] = data_vel[kk] = 0. ;
#else
      data_pos[kk] = 0.f;
#endif
      for(k = 0; k < Prop->npart ; k++)
      {
        data_pos[kk] += (double)((Prop->pos[3*k+i]-Prop->pcm[i]) * (Prop->pos[3*k+j]-Prop->pcm[j]));
#ifdef STORE_VELOCITIES
        data_vel[kk] += (double)((Prop->vel[3*k+i]-Prop->vcm[i]) * (Prop->vel[3*k+j]-Prop->vcm[j])) ;
#endif
      }
      data_pos[kk] = data_pos[kk]/(double)Prop->npart;
#ifdef STORE_VELOCITIES  
      data_vel[kk] = data_vel[kk]/(double)Prop->npart;
#endif
      kk++;
    }

  m = gsl_matrix_view_array (data_pos, 3, 3);
  w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  Prop->cc = (TYPE_REAL) gsl_vector_get (eval, 0);
  Prop->bb = (TYPE_REAL) gsl_vector_get (eval, 1);
  Prop->aa = (TYPE_REAL) gsl_vector_get (eval, 2);

  Prop->cc = sqrt(Prop->cc);
  Prop->bb = sqrt(Prop->bb);
  Prop->aa = sqrt(Prop->aa);

  for(i=0;i<3;i++)
  {
    Prop->evec[2-i][0] = (TYPE_REAL)gsl_matrix_get(evec,0,i);
    Prop->evec[2-i][1] = (TYPE_REAL)gsl_matrix_get(evec,1,i);
    Prop->evec[2-i][2] = (TYPE_REAL)gsl_matrix_get(evec,2,i);
  }

#ifdef STORE_VELOCITIES  
  m = gsl_matrix_view_array (data_vel, 3, 3);
  w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  Prop->cc_vel = (TYPE_REAL) gsl_vector_get (eval, 0);
  Prop->bb_vel = (TYPE_REAL) gsl_vector_get (eval, 1);
  Prop->aa_vel = (TYPE_REAL) gsl_vector_get (eval, 2);

  Prop->cc_vel = sqrt(Prop->cc_vel);
  Prop->bb_vel = sqrt(Prop->bb_vel);
  Prop->aa_vel = sqrt(Prop->aa_vel);

  for(i=0;i<3;i++)
  {
    Prop->evec_vel[2-i][0] = (TYPE_REAL)gsl_matrix_get(evec,0,i);
    Prop->evec_vel[2-i][1] = (TYPE_REAL)gsl_matrix_get(evec,1,i);
    Prop->evec_vel[2-i][2] = (TYPE_REAL)gsl_matrix_get(evec,2,i);
  }
#endif

  gsl_matrix_free (evec);
  gsl_vector_free (eval);
  gsl_eigen_symmv_free (w);
}

extern void properties(struct grup_data *Prop)
{
#ifdef STORE_VELOCITIES  
	TYPE_REAL  dis,_vmax;
	TYPE_REAL  rho;
	TYPE_REAL  _t1, _t2;
	TYPE_INT   dim, i;
	int        indx, i200, ivir;
	TYPE_REAL  dx[3],dv[3];
	TYPE_REAL  a;
	gsl_vector_float *_dis;
	gsl_permutation  *_ind;
#ifdef COMPUTE_EP
	const float  Theta = 0.45;
  TYPE_REAL *Ec, *Ep;
	TYPE_REAL Epmin;
	TYPE_INT j;
  int iepmin;
#endif

	a = (1.0 / (1.0 + cp.redshift));

#ifdef COMPUTE_EP

  Ec = (TYPE_REAL *) malloc(Prop->npart*sizeof(TYPE_REAL));
  Ep = (TYPE_REAL *) malloc(Prop->npart*sizeof(TYPE_REAL));

	/* Si el halo tiene menos de 1000 particulas calcula la energia
	de forma directa (N^2). Si tiene mas de 1000 particulas usa 
	un octree. */
	if(Prop->npart < 1000)
	{
		for(i = 0; i < Prop->npart; i++)
	  {
			/* Calcula energia cinetica (velocidades en [km / seg]) */
			for(dim = 0; dim < 3; dim++)
			{
		  	dv[dim]  = cp.Hubble_a*a*(Prop->pos[3*i+dim]-Prop->pcm[dim])/1000.0f; /* vel de Hubble */
		   	dv[dim] += (TYPE_REAL)sqrt(a)*(Prop->vel[3*i+dim]-Prop->vcm[dim]);
			}

  		Ec[i]  = dv[0]*dv[0];
 		  Ec[i] += dv[1]*dv[1];
 			Ec[i] += dv[2]*dv[2];
 			Ec[i] *= 0.5;

      /* Calcula energia potencial */
			Ep[i]  = 0.;
		  for(j = 0; j < Prop->npart; j++)
		  {
			  if(j == i)continue;

        for(dim = 0; dim < 3; dim++)
          dx[dim] = Prop->pos[3*i+dim] - Prop->pos[3*j+dim];

			  dis  = dx[0]*dx[0];
			  dis += dx[1]*dx[1];
			  dis += dx[2]*dx[2];
				dis  = (TYPE_REAL)sqrt(dis);
				dis += (TYPE_REAL)cp.soft;

			  Ep[i] += 1./dis;
		  }
		  Ep[i] += (1./cp.soft);            /* Autoenergia */
		  Ep[i] *= (GCONS*cp.Mpart*Msol*1.E10/Kpc/a);
			Ep[i] *= (-1.);                   /* Cambio de signo para que Ep sea negativa */
		}

	}else{

    j = 2 * Prop->npart + 200;
		force_treeallocate(j);
		i = force_treebuild(Prop, Theta);

#ifdef DEBUG
    assert(i<j);
#endif

		for(i = 0; i < Prop->npart; i++)
		{

			/* Calcula la energia potencial */
      for(dim = 0; dim < 3; dim++)
        dx[dim] = Prop->pos[3*i+dim];

			Ep[i] = force_treeevaluate_potential(dx);
#ifdef DEBUG
      if(isnan(Ep[i]))exit(EXIT_FAILURE);
#endif
	 		Ep[i] -= cp.Mpart/cp.soft;    /* Autoenergia */
  		Ep[i] *= GCONS/a*Msol*1.E10/Kpc;

			/* Calcula la energia cinetica (velocidades en [km / seg]) */
			for(dim = 0; dim < 3; dim++)
			{
				dv[dim]  = cp.Hubble_a*a*(Prop->pos[3*i+dim]-Prop->pcm[dim])/1000.0f;  /* vel de Hubble */
				dv[dim] += (TYPE_REAL)sqrt(a)*(Prop->vel[3*i+dim]-Prop->vcm[dim]);
			}

	  	Ec[i]  = dv[0]*dv[0];
 			Ec[i] += dv[1]*dv[1];
 			Ec[i] += dv[2]*dv[2];
 			Ec[i] *= 0.5;
		}

		force_treefree();
	}

  Prop->Ep = Prop->Ec = 0.;
  iepmin = -1; Epmin = 0.;
	for(i = 0; i < Prop->npart; i++)
	{
		// Kinetic and potential halo energy
		Prop->Ep += (TYPE_REAL)Ep[i];
		Prop->Ec += (TYPE_REAL)Ec[i];

		if(Ep[i] < Epmin) 
		{
			Epmin  = (TYPE_REAL)Ep[i];
			iepmin = i;
      for(dim = 0; dim < 3; dim++)
        Prop->mostbound[dim] = Prop->pos[3*iepmin+dim];
		}
	}

	Prop->Ep *= 0.5f*(TYPE_REAL)cp.Mpart;
	Prop->Ec *= (TYPE_REAL)cp.Mpart;

#ifdef DEBUG
	assert((0 <= iepmin) && (iepmin<Prop->npart));
#endif

#endif

	/*Ordena por distancia al centro de potencial*/
 	_dis = gsl_vector_float_alloc(Prop->npart);
	_ind = gsl_permutation_alloc(Prop->npart);

	for(i = 0; i < Prop->npart; i++)
	{
		for(dim = 0; dim < 3; dim++)
		{
#ifdef COMPUTE_EP
			dx[dim] = Prop->pos[3*i+dim] - Prop->pos[3*iepmin+dim];
#else                                
      dx[dim] = Prop->pos[3*i+dim] - Prop->pcm[dim]         ;
#endif
			dv[dim] = Prop->vel[3*i + dim] - Prop->vcm[dim];
      Prop->sig[dim] += Prop->vel[3*i + dim] * Prop->vel[3*i + dim];
		}

		Prop->L[0] += dx[1]*dv[2] - dx[2]*dv[1];
    Prop->L[1] += dx[2]*dv[0] - dx[0]*dv[2];
    Prop->L[2] += dx[0]*dv[1] - dx[1]*dv[0];	

		dis  = dx[0]*dx[0];
		dis += dx[1]*dx[1];
		dis += dx[2]*dx[2];
		dis  = (TYPE_REAL)sqrt(fabs(dis));
		gsl_vector_float_set(_dis,i,dis);
	}

	/* Sigma y pasaje a coord fisicas del L, pibe! */
	for(dim = 0; dim < 3; dim++)
	{
    Prop->sig[dim]  = (Prop->sig[dim] - Prop->vcm[dim]*Prop->vcm[dim]*(TYPE_REAL)Prop->npart);
    Prop->sig[dim] /= (TYPE_REAL)(Prop->npart - 1);
    Prop->sig[dim]  = (TYPE_REAL)sqrt(Prop->sig[dim]) ;
		Prop->L[dim]   *= a*(TYPE_REAL)(sqrt(a)*cp.Mpart);
	}

	gsl_sort_vector_float_index(_ind,_dis);

#ifdef COMPUTE_EP
	/* Parametro adimensional de Spin */
	Prop->lambda  = (TYPE_REAL)sqrt(Prop->L[0]*Prop->L[0]+Prop->L[1]*Prop->L[1]+Prop->L[2]*Prop->L[2]);
	Prop->lambda *= (TYPE_REAL)sqrt(fabs(Prop->Ec + Prop->Ep));
	Prop->lambda /= (TYPE_REAL)(GCONS/Kpc*Msol*1.E10);
	Prop->lambda /= (TYPE_REAL)pow((double)Prop->npart*cp.Mpart,2.5);
#endif

	i200 = ivir = -1;
	Prop->vmax = 0.;
	_t1  = (TYPE_REAL)pow(cp.lbox/(TYPE_REAL)cp.npart,3);
	_t1 *= 0.75;
	_t1 /= M_PI;
  i = 4;
  do
	{
    if(i>=Prop->npart)printf("%d %d\n",i,Prop->npart);
    assert(i<Prop->npart);
		indx = gsl_permutation_get(_ind,i);
#ifdef DEBUG
    assert((0<=indx) && (indx<Prop->npart));
#endif
		dis = gsl_vector_float_get(_dis,indx);
		_t2  = (TYPE_REAL)(i+1);
    rho  = _t1;
		rho *= _t2;
		rho /= (TYPE_REAL)(dis*dis*dis);

		if (rho < 200./cp.omegam && i200 == -1) 
			i200 = i;

		if (rho < deltavir(cp.redshift) && ivir == -1) 
			ivir = i;

		_vmax = (TYPE_REAL)(GCONS/Kpc*Msol*(double)(i+1)*cp.Mpart)/dis;
		if(_vmax > Prop->vmax) Prop->vmax = _vmax;

		i++;
	}while(i < Prop->npart);

	/* R200 y M200 */
	if(i200 == -1) i200 = Prop->npart-1;
	indx = gsl_permutation_get(_ind,i200);
	dis  = gsl_vector_float_get(_dis,indx);
	Prop->r200 = dis;
	Prop->m200 = (TYPE_REAL)((double)(i200+1)*cp.Mpart);
	Prop->v200 = (TYPE_REAL)(GCONS/Kpc*Msol*Prop->m200/Prop->r200);
	Prop->v200 = (TYPE_REAL)(sqrt(Prop->v200)*1.E5);

	/* Rvir y Mvir */
	if(ivir == -1) ivir = Prop->npart-1;
	indx = gsl_permutation_get(_ind,ivir);
	dis  = gsl_vector_float_get(_dis,indx);
	Prop->rvir = dis;
	Prop->mvir = (TYPE_REAL)((double)(ivir+1)*cp.Mpart);
	Prop->vvir = (TYPE_REAL) (GCONS/Kpc*Msol*Prop->mvir/Prop->rvir);
	Prop->vvir = (TYPE_REAL) (sqrt(Prop->vvir)*1.E5);

	/* Velocidad circular maxima */
	Prop->vmax = (TYPE_REAL) (sqrt(Prop->vmax)*1.E5);

	gsl_vector_float_free(_dis);
	gsl_permutation_free(_ind);

#ifdef COMPUTE_EP
  free(Ec);
  free(Ep);
#endif

#endif

	/* Autovalores del tensor de forma */
	forma(Prop);

  return;
}
