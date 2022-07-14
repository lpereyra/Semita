#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.hh"
#include "timer.hh"
#include "colors.hh"
#include "io.hh"
#include "allocate.hh"
#include "leesnap.hh"
#include "iden.hh"
#include "voronoi.hh"
#include "kruskal.hh"
#include "select.hh"
#include "libfitpack.hh"

int main(int argc, char **argv)
{
  double start, end;
  double MassCut;

  omp_set_nested(1);

  TIMER(start);
  
  init_variables(argc,argv);

  MassCut = atof(argv[2]);

  /* Read Simulation */
  read_gadget();

#ifdef CHANGE_POSITION
  RED("\nBegins Change Positions\n");
  fflush(stdout);
  change_positions(cp.npart);
  GREEN("\nEnd Change Positions\n");
  fflush(stdout);
#endif

  fprintf(stdout, "\nBegins Identification\n");
  fflush(stdout);
  identification();
  fprintf(stdout, "\nEnds Identification\n");
  fflush(stdout);
  
	std::vector<std::pair<TYPE_REAL,std::pair<TYPE_INT,TYPE_INT> > > edges;
  std::vector<std::vector<TYPE_INT> > adjacency_list;
	TYPE_INT NumPartCut;

	Voronoi_Grupos(fof[0],edges);

  Kruskal(edges,adjacency_list);

  NumPartCut = (TYPE_INT)(MassCut/cp.Mpart) ;  // Masa de la partícula [10^10 Msol / h]
  
  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Mass cut %g -> Nodos Npart %d Mass %g\n", \
  MassCut*1.e10,NumPartCut,NumPartCut*cp.Mpart*1e10);
	RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

	#ifdef PRUNED
  	#ifdef BRANCH_SURVIVE
	  	Podado(LEVEL_PRUNED,adjacency_list,NumPartCut);  
		  fprintf(stdout,"Sobreviven en ramas los nodos con %d\n",NumPartCut);
  		fflush(stdout);
		#else
	  	Podado(LEVEL_PRUNED,adjacency_list);  
  	#endif
  #endif

  Selection(NumPartCut,adjacency_list);
	suaviza();
	write_seg(NumPartCut, fof);

  free(Gr);
  free_particles(&P);
  edges.clear();
  adjacency_list.clear();
  TIMER(end);
  
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}
