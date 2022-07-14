#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "variables.hh"
#include "deltas.hh"

extern double deltavir( const double z )
{
  double q ;
  double as, x ;
	double tmp;
  as=1./(1.+z);
  if((cp.omegam-1. < 1.E-9) & (cp.omegal < 1.E-9))
  {
    q = 18.*M_PI*M_PI;      
    /*Einstein - de Sitter Model Peebles, 1980. Eq 19.50*/
  }
  if((cp.omegam < 1.) & (cp.omegal < 1.E-9))
  {
    x = as * 2. *(1./cp.omegam-1.);
		tmp = sqrt(x*x+2.*x);
		tmp-= log((1.+x)+sqrt(x*x+2.*x));
    q = 4.*M_PI*M_PI;
		q*= (x*x*x);
		q/= (tmp*tmp);
    /*Maoz, 1990 Eq 12. Lacey & Cole, 1993 Eq A15. Oukbir & Blanchard, 1997*/
  }
  if((cp.omegam < 1.) & (cp.omegak < 1.E-9))
  {
		x = cp.omegam*pow(1.+z,3);
		x+= cp.omegak*pow(1+z,2);
		x+= cp.omegal;
		x = cp.omegam*pow(1+z,3)/x;
		q = (18.*M_PI*M_PI+82.*(x-1.)-39.*pow(x-1.,2))/x;
		/* Bryan & Norman 1998ApJ,495,80B, Bolshoi appendix*/
  }
  return q;
}
