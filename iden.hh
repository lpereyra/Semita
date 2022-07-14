#ifndef IDENTIFICADOR_H
#define IDENTIFICADOR_H

struct iden_st
{
  double *r0;        /*Linking Length para el FoF*/
  TYPE_INT nobj;
  TYPE_INT ngrupos;
};

struct temporary 
{
  TYPE_INT *head;
  TYPE_INT *ll;
  TYPE_INT *npgrup;
};

extern TYPE_INT Find(TYPE_INT, TYPE_INT * __restrict__);
extern void identification(void);

#endif
