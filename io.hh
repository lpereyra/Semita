#ifndef IO_H
#define IO_H

extern void init_variables(int , char **);
extern void write_properties(FILE *, struct grup_data);
#ifdef FILE_ASCII
  extern void write_properties_ascii(FILE *, struct grup_data);
#endif
extern void write_seg(const TYPE_INT NNN, const TYPE_REAL * fof);
 
#endif
