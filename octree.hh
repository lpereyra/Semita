#ifndef OCTREE_H
#define OCTREE_H

struct NODE{ 
  TYPE_REAL s[3];                     /* center of mass */
  long  partind;
  TYPE_REAL center[3],len;            /* center and sidelength of treecubes */
  float mass;                         /* mass*/
  float oc;                           /* variable for opening criter*/
  struct NODE *next,*sibling,*father,*suns[8];
};

extern void force_treeallocate(int);
extern void force_treefree(void);

extern int force_treebuild(struct grup_data *, float);
extern TYPE_REAL force_treeevaluate_potential(const TYPE_REAL []);

#endif
