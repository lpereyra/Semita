#ifndef VOL_H 
#define VOL_H

#ifdef WITH_EXTREMES
	extern TYPE_REAL volume(const double R_TUBE, const TYPE_INT nsize, TYPE_REAL *Pos, TYPE_REAL Rvir[2]);
#else
	extern TYPE_REAL volume(const double R_TUBE, const TYPE_INT nsize, TYPE_REAL *Pos);
#endif

#endif
