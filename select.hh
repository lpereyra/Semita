#ifndef PRUNE 
#define PRUNE 

#ifdef PRUNED
#ifdef BRANCH_SURVIVE
	extern void Podado(const TYPE_INT level, std::vector<std::vector<TYPE_INT> > &vec, const TYPE_INT NumPartCut);
#else
  extern void Podado(const TYPE_INT level, std::vector<std::vector<TYPE_INT> > &vec);
#endif
#endif
	extern void Selection(const TYPE_INT NumPartCut, std::vector<std::vector<TYPE_INT> > &adjacency_list);

#endif
