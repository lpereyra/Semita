#ifndef ALLOCATE_H
#define ALLOCATE_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

#ifdef COLUMN
  extern int allocate_particles(struct particle_data  *, const TYPE_INT);
  extern int reallocate_particles(struct particle_data  *, const TYPE_INT);
  extern void free_particles(struct particle_data  *);
#else
  extern int allocate_particles(struct particle_data **, const TYPE_INT);
  extern int reallocate_particles(struct particle_data **, const TYPE_INT);
  extern void free_particles(struct particle_data **);
#endif

#endif
