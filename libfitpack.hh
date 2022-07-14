#ifndef SMOOTH_H
#define SMOOTH_H

extern "C" {
  extern void splev_(double*,int*,double*,int*,double*,double*,int*,int*,int*);
  extern void parcur_(int*,int*,int*,int*,double*,int*,double*, \
              double*,double*,double*,int*,double*,int*,int*, \
              double*,int*,double*,double*,double*,int*,int*,int*);
}

extern void suaviza(void);

#endif
