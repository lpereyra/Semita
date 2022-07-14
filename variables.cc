#include "variables.hh"

struct cosmoparam cp;
struct SnapST snap;
struct gridst grid;
#ifdef COLUMN
  struct particle_data P;
#else
  struct particle_data *P;
#endif
struct grup_data *Gr;
struct segmentstd *Seg;
TYPE_INT  nfrac;
TYPE_REAL *fof;
#ifdef CHANGE_POSITION
  TYPE_REAL pmin[3], pmax[3];
#endif
char message[200];
