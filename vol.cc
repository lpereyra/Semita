#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "variables.hh"
#include "vol.hh"

#define NPIX 10
#define NumberOfSides 1000

inline static double square( double f ) { return (f*f);}

inline static double Dot(double v[], double w[])
{
	int idim;
	double dot = 0.0;
	
	for(idim = 0; idim<3; idim++)
		dot += v[idim]*w[idim];

	return dot;
}

inline static double Normalize(double v[])
{
	int idim;
	double r;

	r = sqrt(Dot(v,v));
	for(idim = 0; idim<3; idim++)
		v[idim] *= 1.0/r;

	return r;
}	

inline static void Cross(const double vect_A[], const double vect_B[], double cross_P[])
{
  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

	return;
}

static void theoretical_volume(const double R_TUBE, const TYPE_INT numPts, double *points, double *loc_vol)
{
  TYPE_INT j;
	double l;
  double p[3];
  double SurfaceArea;
  double Volume;
	const double R_TUBE_2 = square(R_TUBE);
	
  // Initialize variables ...
  Volume = 0.0;
  SurfaceArea = 0.0;

	for(j=0; j<3; j++)
		p[j] = points[4+j] - points[j];
	
	l = sqrt(Dot(p,p));

	SurfaceArea += 2*M_PI*R_TUBE*l;
	Volume += M_PI*R_TUBE_2*l;

#ifdef WITH_EXTREMES

	for(j=0;j<2;j++)
	{
		l = sqrt(fabs(points[4*j+3]*points[4*j+3] - R_TUBE_2));
		SurfaceArea -= 2*M_PI*R_TUBE*l;
		Volume -= M_PI*R_TUBE_2*l;
		
		//Vol Spherical Cap
		l = points[4*j+3] - l;
		Volume -= (M_PI/6)*l*(3*R_TUBE_2 + l*l);
	}

#endif

	*loc_vol = *loc_vol + Volume;

	return;
}

static void local_volume(const double R_TUBE, const TYPE_INT numPts, double *points, double *loc_vol)
{
  TYPE_INT i, j, k, numTri;
  const double Theta = 2.0 * M_PI / NumberOfSides;
  double vol[3], kxyz[3];
  double munc[3], wxyz, wxy, wxz, wyz;
  double a, b, c, ss, area;
  double x[3], y[3], z[3];
  double ei[3], ej[3], ek[3], u[3], absu[3], length;
  double ii[3], jj[3], kk[3];
  double xavg, yavg, zavg;
  double pNext[3];
  double sNext[3];
  double sPrev[3];
  double p[3];
  double s[3];
  double w[3];
  double nP[3];
  double n[3];
  double default_normal[3] = { 0.0, 0.0, 1.0 };
  double Volume;
  double SurfaceArea;
	int i0, i1, i2;
	int tmp_id[3];
	double *newPts;
	double *normals;
	double normal[3];
	const double R_TUBE_2 = R_TUBE*R_TUBE;

  // Initialize variables ...
  Volume = 0.0;
  SurfaceArea = 0.0;

	if(numPts <= 1)
	{
		
		*loc_vol = *loc_vol + Volume;

		return;
	
	}else if(numPts == 2){

		for(i=0; i<3; i++)
			p[i] = points[3+i] - points[i];
	
		ss = sqrt(Dot(p,p));

		SurfaceArea += 2*M_PI*R_TUBE*ss;
		Volume += M_PI*R_TUBE_2*ss;

		*loc_vol = *loc_vol + Volume;

		return;

	}else{

		normals = (double *) malloc(3*numPts*sizeof(double));
		newPts = (double *) malloc(3*numPts*NumberOfSides*sizeof(double));
		
		for (j = 0; j < numPts; j++)
		{
			normals[3*j]   = default_normal[0];
			normals[3*j+1] = default_normal[1];
			normals[3*j+2] = default_normal[2];
		}

    for (i = 0; i < 3; i++)
		  sPrev[i] = points[3+i] - points[i];
		if (Normalize(sPrev) < 1e-8)
		{
		  fprintf(stdout,"Coincident points!\n");
			fflush(stdout);
		}

	  j = 0;
    while (++j < numPts)
    {
	    for (i = 0; i < 3; i++)
			  sNext[i] = points[3*(j+1)+i] - points[3*j+i];
			if (Normalize(sNext) < 1e-8)
			{
			  fprintf(stdout,"Coincident points!\n");
				fflush(stdout);
			}

      // now the starting normal should simply be the cross product
      // in the following if statement we check for the case where
      // the two segments are parallel, in which case, continue searching
      // for the next valid segment
		  Cross(sPrev, sNext, n);
			ss = sqrt(Dot(n,n));
      if (ss > 1.0E-3)
      {
	    	for (i = 0; i < 3; i++)
				{
        	normal[i] = n[i];
        	sPrev[i] = sNext[i];
				}
        break;
      }
    }

    if (j >= numPts) // only one valid segment
    {
      // a little trick to find othogonal normal
      for (i = 0; i < 3; ++i)
      {
        if (sPrev[i] != 0.0)
        {
          normal[(i + 2) % 3] = 0.0;
          normal[(i + 1) % 3] = 1.0;
          normal[i] = -sPrev[(i + 1) % 3] / sPrev[i];
          break;
        }
      }
    }
 
		if (Normalize(normal) < 1e-8)
		{
		  fprintf(stdout,"Bad normal n = (%g, %g, %g)\n", 
			normal[0], normal[1], normal[2]);
			fflush(stdout);
		}

    // compute remaining normals
    k = 0;
    while (++j < numPts)
    {
			if(j == numPts-1) 
			{
				break;
			}

	    for (i = 0; i < 3; i++)
			  sNext[i] = points[3*(j+1)+i] - points[3*j+i];
			if (Normalize(sNext) < 1e-8)
			{
				continue;
			}

      // compute rotation vector
		  Cross(sPrev, normal, w);
      if (Normalize(w) < 1e-8) // can't use this segment
      {
        continue;
      }

      // compute rotation of line segment
      Cross(sNext, sPrev, p);
      if (Normalize(p) < 1e-8) // can't use this segment
      {
        continue;
      }

      a = Dot(p, normal);
      b = 1.0 - square(a);
      if (b > 0.0)
      {
        b = sqrt(1.0 - square(a));
      }else{
        b = 0.0;
      }

		  for (i = 0; i < 3; i++)
			  s[i] = sNext[i] + sPrev[i];
      if (Normalize(s) < 1e-8) // can't use this segment
      {
        continue;
      }

			Cross(s, p, w);
      Cross(sPrev, p, s);
      if ((Dot(normal, s) * Dot(w, s)) < 0)
      {
        b = -1.0 * b;
      }

      // insert current normal before updating
      for (i = k; i < j; ++i)
      {
				normals[3*i]   = normal[0];
				normals[3*i+1] = normal[1];
				normals[3*i+2] = normal[2];
      }
      
			k = j;
		  for (i = 0; i < 3; i++)
			{
      	sPrev[i] = sNext[i];
      	// compute next normal
      	normal[i] = (a * p[i]) + (b * w[i]);
   		}
    }

    // insert last normal for the remaining points
    for (i = k; i < numPts; ++i)
    {
			normals[3*i]   = normal[0];
			normals[3*i+1] = normal[1];
			normals[3*i+2] = normal[2];
    }

		// Create points along each polyline that are
		// connected into NumberOfSides triangle strips. 
		// Use "averaged" segment to create beveled effect.
		// Watch out for first and last points.
		for (j = 0; j < numPts; j++)
		{
		  if (j == 0) // first point
		  {
		    p[0] = points[0];
		    p[1] = points[1];
		    p[2] = points[2];
		
				pNext[0] = points[3];
		    pNext[1] = points[4];
		    pNext[2] = points[5];
		
		    for (i = 0; i < 3; i++)
		    {
		      sNext[i] = pNext[i] - p[i];
		      sPrev[i] = sNext[i];
		    }
		  }
		  else if (j == (numPts - 1)) // last point
		  {
		    for (i = 0; i < 3; i++)
		    {
		      sPrev[i] = sNext[i];
		      p[i] = pNext[i];
		    }
		  }
		  else
		  {
		    for (i = 0; i < 3; i++)
		    {
		      p[i] = pNext[i];
		    }
	
				pNext[0] = points[3*(j+1)];
		    pNext[1] = points[3*(j+1)+1];
		    pNext[2] = points[3*(j+1)+2];
		
		    for (i = 0; i < 3; i++)
		    {
		      sPrev[i] = sNext[i];
		      sNext[i] = pNext[i] - p[i];
		    }
		  }
		
		  if (Normalize(sNext) < 1e-8)
		  {
		    fprintf(stdout,"Coincident points!\n");
				fflush(stdout);
		  }

		  for (i = 0; i < 3; i++)
			{
				normal[i] = normals[3*j+i];
		    s[i] = (sPrev[i] + sNext[i]) * 0.5; // average vector
			}

		  // if s is zero then just use sPrev cross n
		  if (Normalize(s) < 1e-8)
		  {
		    fprintf(stdout,"Using alternate bevel vector\n");
				fflush(stdout);
		    
				Cross(sPrev, normal, s);
		    if (Normalize(s) < 1e-8)
		    {
		      fprintf(stdout,"Using alternate bevel vector\n");
					fflush(stdout);
		    }
		  }
		
		  Cross(s, normal, w);
		  if (Normalize(w) < 1e-8)
		  {
		    fprintf(stdout,"Bad normal s = (%g, %g, %g), n = (%g, %g, %g)\n",
				s[0], s[1], s[2], normal[0], normal[1], normal[2]);
				fflush(stdout);
		  }
		
		  Cross(w, s, nP); // create orthogonal coordinate system
		  Normalize(nP);
		
		  for (k = 0; k < NumberOfSides; k++)
		  {
		    for (i = 0; i < 3; i++)
		    {
		      n[i] = w[i] * cos((double)k * Theta) + nP[i] * sin((double)k * Theta);
		      u[i] = p[i] + R_TUBE * n[i];
					newPts[3*(j*NumberOfSides+k)+i] = u[i];
		    }
		  } // for each side
		} // for all points in polyline

	  wxyz = 0;
	  wxy = 0.0;
	  wxz = 0.0;
	  wyz = 0.0;
	  for(i = 0; i < 3; i++)
	  {
	    munc[i] = 0.0;
	    vol[i]  = 0.0;
	    kxyz[i] = 0.0;
	  }
	
		numTri = 0;
		for(k = 0; k < NumberOfSides; k++)
		{
		  tmp_id[0] = k % NumberOfSides;
		  tmp_id[1] = (k + 1) % NumberOfSides;
	
		  for (i = 0; i < (numPts-1); i++)
		  {
		  	tmp_id[2] = i*NumberOfSides;
				for(j=0;j<2;j++)
				{
					if(j == 0) // Oriented Triangle
					{
						i0 = tmp_id[0] + tmp_id[2];
						x[0] = newPts[3*i0];
						y[0] = newPts[3*i0+1];
						z[0] = newPts[3*i0+2];
		
						i1 = tmp_id[1] + tmp_id[2];
						x[1] = newPts[3*i1];
						y[1] = newPts[3*i1+1];
						z[1] = newPts[3*i1+2];
	
						i2 = tmp_id[0] + tmp_id[2] + NumberOfSides;
						x[2] = newPts[3*i2];
						y[2] = newPts[3*i2+1];
						z[2] = newPts[3*i2+2];
					}else{ //Not Oriented Triangle
	
						i0 = tmp_id[0] + tmp_id[2] + NumberOfSides;
						x[0] = newPts[3*i0];
						y[0] = newPts[3*i0+1];
						z[0] = newPts[3*i0+2];
		
						i1 = tmp_id[1] + tmp_id[2];
						x[1] = newPts[3*i1];
						y[1] = newPts[3*i1+1];
						z[1] = newPts[3*i1+2];
	
						i2 = tmp_id[1] + tmp_id[2] + NumberOfSides;
						x[2] = newPts[3*i2];
						y[2] = newPts[3*i2+1];
						z[2] = newPts[3*i2+2];
					}
	
					// get i j k vectors ...
					ei[0] = (x[1] - x[0]);
					ej[0] = (y[1] - y[0]);
					ek[0] = (z[1] - z[0]);
					ei[1] = (x[2] - x[0]);
					ej[1] = (y[2] - y[0]);
					ek[1] = (z[2] - z[0]);
					ei[2] = (x[2] - x[1]);
					ej[2] = (y[2] - y[1]);
					ek[2] = (z[2] - z[1]);
		
					// cross product between two vectors, to determine normal vector
					u[0] = (ej[0] * ek[1] - ek[0] * ej[1]);
					u[1] = (ek[0] * ei[1] - ei[0] * ek[1]);
					u[2] = (ei[0] * ej[1] - ej[0] * ei[1]);
		
					// normalize normal vector to 1
					length = sqrt(Dot(u, u));
					if (length != 0.0)
					{
					  u[0] /= length;
					  u[1] /= length;
					  u[2] /= length;
					}
					else
					{
					  u[0] = u[1] = u[2] = 0.0;
					}
		
					// determine max unit normal component...
					absu[0] = fabs(u[0]);
					absu[1] = fabs(u[1]);
					absu[2] = fabs(u[2]);
		
					if ((absu[0] > absu[1]) && (absu[0] > absu[2]))
					{
					  munc[0]++;
					}
					else if ((absu[1] > absu[0]) && (absu[1] > absu[2]))
					{
					  munc[1]++;
					}
					else if ((absu[2] > absu[0]) && (absu[2] > absu[1]))
					{
					  munc[2]++;
					}
					else if ((absu[0] == absu[1]) && (absu[0] == absu[2]))
					{
					  wxyz++;
					}
					else if ((absu[0] == absu[1]) && (absu[0] > absu[2]))
					{
					  wxy++;
					}
					else if ((absu[0] == absu[2]) && (absu[0] > absu[1]))
					{
					  wxz++;
					}
					else if ((absu[1] == absu[2]) && (absu[0] < absu[2]))
					{
					  wyz++;
					}
					else
					{
					  fprintf(stdout,"Unpredicted situation...!\n");
						fflush(stdout);
					}
		
					// This is reduced to ...
					ii[0] = ei[0] * ei[0];
					ii[1] = ei[1] * ei[1];
					ii[2] = ei[2] * ei[2];
					jj[0] = ej[0] * ej[0];
					jj[1] = ej[1] * ej[1];
					jj[2] = ej[2] * ej[2];
					kk[0] = ek[0] * ek[0];
					kk[1] = ek[1] * ek[1];
					kk[2] = ek[2] * ek[2];
		
					// area of a triangle...
					a = sqrt(ii[1] + jj[1] + kk[1]);
					b = sqrt(ii[0] + jj[0] + kk[0]);
					c = sqrt(ii[2] + jj[2] + kk[2]);
					ss = 0.5 * (a + b + c);
					area = sqrt(fabs(ss * (ss - a) * (ss - b) * (ss - c)));
				  // Surface Area ...
					SurfaceArea += area;
					
					// volume elements ...
					zavg = (z[0] + z[1] + z[2]) / 3.0;
					yavg = (y[0] + y[1] + y[2]) / 3.0;
					xavg = (x[0] + x[1] + x[2]) / 3.0;
		
					vol[2] += (area * u[2] * zavg);
					vol[1] += (area * u[1] * yavg);
					vol[0] += (area * u[0] * xavg);
					numTri++;
				}
		  }
		}

		// Weighting factors in Discrete Divergence theorem for volume calculation.
		//
		//numTri = 2*NumberOfSides*(numPts-1);
		kxyz[0]  = (munc[0] + (wxyz / 3.0) + ((wxy + wxz) / 2.0)) / (double)numTri;
		kxyz[1]  = (munc[1] + (wxyz / 3.0) + ((wxy + wyz) / 2.0)) / (double)numTri;
		kxyz[2]  = (munc[2] + (wxyz / 3.0) + ((wxz + wyz) / 2.0)) / (double)numTri;
		Volume  += fabs(kxyz[0] * vol[0] + kxyz[1] * vol[1] + kxyz[2] * vol[2]);

		free(normals);
		free(newPts);
	
		*loc_vol = *loc_vol + Volume;

		return;
	}
}

#ifdef WITH_EXTREMES

static double sphere_line_intersection(const double *points, const TYPE_INT numPts, const TYPE_INT index, const TYPE_INT flag, const double r, const double R_TUBE)
{
  double a, b, c, i, mu;
  double u[3], p[3];
	double p1[3], p2[3];
	double *tmp = NULL;
	double vol = 0.0;
	TYPE_INT k, l;

	if(flag == 0)
	{
		for(l = 0; l < 3; l++)
		{
  		p1[l] = points[3*index+l]     - points[l];
#ifdef PERIODIC
 	    if(p1[l]> 0.5*cp.lbox) p1[l] -= cp.lbox;
  	  if(p1[l]<-0.5*cp.lbox) p1[l] += cp.lbox;
#endif
  		p2[l] = points[3*(index+1)+l] - points[l];
#ifdef PERIODIC
 	    if(p2[l]> 0.5*cp.lbox) p2[l] -= cp.lbox;
  	  if(p2[l]<-0.5*cp.lbox) p2[l] += cp.lbox;
#endif
		}

	}else{

		for(l = 0; l < 3; l++)
		{
  		p1[l] = points[3*index+l]     - points[3*(numPts-1)+l];
#ifdef PERIODIC
 	    if(p1[l]> 0.5*cp.lbox) p1[l] -= cp.lbox;
  	  if(p1[l]<-0.5*cp.lbox) p1[l] += cp.lbox;
#endif
  		p2[l] = points[3*(index-1)+l] - points[3*(numPts-1)+l];
#ifdef PERIODIC
 	    if(p2[l]> 0.5*cp.lbox) p2[l] -= cp.lbox;
  	  if(p2[l]<-0.5*cp.lbox) p2[l] += cp.lbox;
#endif
		}
	}

  u[0] = (p2[0] - p1[0]);
  u[1] = (p2[1] - p1[1]); 
  u[2] = (p2[2] - p1[2]);

  a =  square(u[0]) + square(u[1]) + square(u[2]);
  b =  2 * Dot(u, p1);
  c =  square(p1[0]) + square(p1[1]) + square(p1[2]) - square(r) ;
  i =  square(b) - 4 * a * c ;
	
  if ( i > 0.0 )
  {
    // two intersections
    mu = -b;
		if(mu > 0)
		{	
			mu -= sqrt(i);
		}else{
			mu += sqrt(i);
		}
		mu /= (2*a);

		for(l=0; l<3; l++)
	    p[l] = p1[l] + mu*u[l];

  }else if ( i < 0.0 ){
    // no intersection

		fprintf(stdout,"Error: No intersection!\n");
		fflush(stdout);
		return 0;
 
  }else{
    // one intersection
    mu = -b/(2*a);
		for(l=0; l<3; l++)
	    p[l] = p1[l] + mu*u[l];
  }

	if(flag == 0)
	{
		tmp = (double *) malloc(3*(index+1)*sizeof(double));
    
		for(k=0;k<index;k++)
		{
			for(l = 0; l < 3; l++)
			{
  			tmp[3*k+l] = points[3*k+l] - points[l];
#ifdef PERIODIC
 	    	if(tmp[3*k+l]> 0.5*cp.lbox) tmp[3*k+l] -= cp.lbox;
  	  	if(tmp[3*k+l]<-0.5*cp.lbox) tmp[3*k+l] += cp.lbox;
#endif
			}
		}

		for(l=0; l<3; l++)
		{
			tmp[3*index+l] = p[l] - tmp[l];
#ifdef PERIODIC
 	    if(tmp[3*index+l]> 0.5*cp.lbox) tmp[3*index+l] -= cp.lbox;
  	  if(tmp[3*index+l]<-0.5*cp.lbox) tmp[3*index+l] += cp.lbox;
#endif
		}

		vol = 0.0;
		local_volume(R_TUBE, index+1, tmp, &vol);

	}else{

		tmp = (double *) malloc(3*(numPts-index+1)*sizeof(double));
    
		for(k=numPts; k>index; k--)
		{
			for(l = 0; l < 3; l++)
			{
  			tmp[3*(numPts-k)+l] = points[3*(k-1)+l] - points[3*(numPts-1)+l];
#ifdef PERIODIC
 	    	if(tmp[3*(numPts-k)+l]> 0.5*cp.lbox) tmp[3*(numPts-k)+l] -= cp.lbox;
  	  	if(tmp[3*(numPts-k)+l]<-0.5*cp.lbox) tmp[3*(numPts-k)+l] += cp.lbox;
#endif
			}
		}

		for(l=0; l<3; l++)
		{
			tmp[3*(numPts-index)+l] = p[l] - tmp[l];
#ifdef PERIODIC
 	    if(tmp[3*(numPts-index)+l]> 0.5*cp.lbox) tmp[3*(numPts-index)+l] -= cp.lbox;
  	  if(tmp[3*(numPts-index)+l]<-0.5*cp.lbox) tmp[3*(numPts-index)+l] += cp.lbox;
#endif
		}

		vol = 0.0;
		local_volume(R_TUBE, numPts-index+1, tmp, &vol);
	}

	free(tmp);

  return vol;
}

#endif

#ifdef WITH_EXTREMES
	extern TYPE_REAL volume(const double R_TUBE, const TYPE_INT nsize, TYPE_REAL *Pos, TYPE_REAL Rvir[2])
#else
	extern TYPE_REAL volume(const double R_TUBE, const TYPE_INT nsize, TYPE_REAL *Pos)
#endif
{
  TYPE_INT j, k, l;
	double r[3];
  double *points, vv;
	TYPE_REAL vol;

	if(nsize == 2)
	{
		points = (double *) malloc(4*nsize*sizeof(double));

  	for(k=0;k<nsize;k++)
		{
			points[4*k+3] = 0.0;
		 	for(l = 0; l < 3; l++)
			{
  	 		points[4*k+l] = Pos[3*k+l];
#ifdef PERIODIC
    		r[l] = Pos[3*k+l] - Pos[l];
 	    	if(r[l]> 0.5*cp.lbox) points[4*k+l] -= cp.lbox;
  	  	if(r[l]<-0.5*cp.lbox) points[4*k+l] += cp.lbox;
#endif
			}
#ifdef WITH_EXTREMES
  	 	points[4*k+3] += Rvir[k];
#endif
		}

		vv = 0.0;
		theoretical_volume(R_TUBE, nsize, points, &vv);
		vol = (TYPE_REAL)vv;

	}else{

		points = (double *) malloc(3*(NPIX*(nsize-1) + 1)*sizeof(double));
    
		for(k=0;k<nsize;k++)
		{
	  	for(l = 0; l < 3; l++)
			{
  	  	points[3*NPIX*k+l] = Pos[3*k+l] - Pos[l];
#ifdef PERIODIC
    		r[l] = points[3*NPIX*k+l];
 	    	if(r[l]> 0.5*cp.lbox) points[3*NPIX*k+l] -= cp.lbox;
  	  	if(r[l]<-0.5*cp.lbox) points[3*NPIX*k+l] += cp.lbox;
#endif
			}
		}

    for(k=1;k<nsize;k++)
    {
		  for(l = 0; l < 3; l++)
  	  	r[l] = points[3*NPIX*k+l] - points[3*NPIX*(k-1)+l];

      vv = Normalize(r);

			for(j=1;j<NPIX;j++)
		  	for(l = 0; l < 3; l++)
  	  		points[3*(NPIX*(k-1)+j)+l] = points[3*NPIX*(k-1)+l] + ((double)j*vv/(double)NPIX)*r[l];
    }

		vv = 0.0;
		local_volume(R_TUBE, NPIX*(nsize-1) + 1, points, &vv);
		vol = (TYPE_REAL)vv;
		
		/////////////////////////////////////////

#ifdef WITH_EXTREMES

		for(k=0; k<(NPIX*(nsize-1)); k++)
		{
			j = l = 0;
			r[0] = points[3*k]   - points[0];
			r[1] = points[3*k+1] - points[1];
			r[2] = points[3*k+2] - points[2];

			vv = sqrt(Dot(r,r));

			if(vv < Rvir[0])
			{
				j = 1;
			}

			r[0] = points[3*(k+1)]   - points[0];
			r[1] = points[3*(k+1)+1] - points[1];
			r[2] = points[3*(k+1)+2] - points[2];

			vv = sqrt(Dot(r,r));

			if(vv < Rvir[0])
			{
				l = 1;
			}
			
			if(j == 1 && l == 0)
				break;
		}

		if(j == 1 && l == 1)
		{
			// All points inside Rvir[0]
			vol  = 0.0f;			
		}else{
			vv = sphere_line_intersection(points, NPIX*(nsize-1) + 1, k, 0, Rvir[0], R_TUBE);
			vol -= (TYPE_REAL)vv;
		}
	
		/////////////////////////////////////////////
	
		for(k=NPIX*(nsize-1)+1; k>1; k--)
		{
			j = l = 0;
			r[0] = points[3*(k-1)]   - points[3*NPIX*(nsize-1)];
			r[1] = points[3*(k-1)+1] - points[3*NPIX*(nsize-1)+1];
			r[2] = points[3*(k-1)+2] - points[3*NPIX*(nsize-1)+2];

			vv = sqrt(Dot(r,r)); 

			if(vv < Rvir[1])
			{
				j = 1;
			}

			r[0] = points[3*(k-2)]   - points[3*NPIX*(nsize-1)];
			r[1] = points[3*(k-2)+1] - points[3*NPIX*(nsize-1)+1];
			r[2] = points[3*(k-2)+2] - points[3*NPIX*(nsize-1)+2];

			vv = sqrt(Dot(r,r)); 

			if(vv < Rvir[1])
			{
				l = 1;
			}

			if(j == 1 && l == 0)
				break;
		}

		if(j == 1 && l == 1)
		{
			// All points inside Rvir[1]
			vol  = 0.0f;			
		}else{
			vv = sphere_line_intersection(points, NPIX*(nsize-1) + 1, k-1, 1, Rvir[1], R_TUBE);
			vol -= (TYPE_REAL)vv;
		}

#endif 

	}

	vol = vol < 0.0 ? 0.0 : vol;

  free(points);

  return vol;
}
