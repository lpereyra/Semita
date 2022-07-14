#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <functional>
#include <algorithm>

#include "variables.hh"
#include "leesnap.hh"
#include "kruskal.hh"
#include "colors.hh"
#include "iden.hh"

extern void Kruskal(std::vector<std::pair<TYPE_REAL,std::pair<TYPE_INT,TYPE_INT> > > &edges,\
										std::vector<std::vector<TYPE_INT> > &adjacency_list)
{
  TYPE_INT i,j,k,id,idv,count;
  TYPE_INT *Root, *Rank;
  char filename[200];
  FILE *pfout;

  RED("Allocating Adjacency List..\n");
	Root = (TYPE_INT *) malloc(cp.ngrup*sizeof(TYPE_INT));
  Rank =  (TYPE_INT *) malloc(cp.ngrup*sizeof(TYPE_INT));

  for(i=0;i<cp.ngrup;i++)
  {
    Root[i] = i;
    Rank[i] = 0;
    adjacency_list.push_back(std::vector<TYPE_INT>());
  }
  
	RED("Build MST\n");

  sprintf(filename,"%.2d_mst_%.2f_%.2f.bin",snap.num,fof[0],fof[1]);

  count = 0;
  pfout=fopen(filename,"w");
  fwrite(&count,sizeof(TYPE_INT),1,pfout);

  sort(edges.begin(),edges.end(), std::greater<std::pair<TYPE_REAL,std::pair<TYPE_INT,TYPE_INT> > >());

  while(!edges.empty()) 
  {
    j = edges.back().second.first;
    k = edges.back().second.second;
    edges.pop_back();

    id = Find(j,Root);
    idv = Find(k,Root);

    if(id!=idv)
    {
      if(Rank[id] < Rank[idv])
        Root[id] = idv;
      else if(Rank[idv] < Rank[id])
        Root[idv] = id;
      else
      {
        Root[idv] = id;
        Rank[id]++;
      }

	    adjacency_list[j].push_back(k);
  	  adjacency_list[k].push_back(j);

      fwrite(&j,sizeof(TYPE_INT),1,pfout);
      fwrite(&k,sizeof(TYPE_INT),1,pfout);
      count++;
    } 

  }

  fprintf(stdout,"%u NumEdges in MST\n",count);
  fflush(stdout);

  rewind(pfout);
  fwrite(&count,sizeof(TYPE_INT),1,pfout);
  fclose(pfout);
  
	free(Root);
  free(Rank);

  return;
}
