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
#include "colors.hh"
#include "leesnap.hh"
#include "select.hh"

#ifdef PRUNED

	#ifdef BRANCH_SURVIVE
	  static TYPE_INT DFS(const TYPE_INT i, std::vector<std::vector<TYPE_INT> > &vec, \
	  TYPE_INT cut, const TYPE_INT NumPartCut)
	#else
    static TYPE_INT DFS(const TYPE_INT i, std::vector<std::vector<TYPE_INT> > &vec, TYPE_INT cut)
  #endif
	{
	  TYPE_INT j,k,id,idv;
	
	  j = 1;
	  id = i;
	  idv = vec[id][0];
	
	  while(vec[idv].size()==2)
	  {
	
	    #ifdef BRANCH_SURVIVE
	      if(Gr[idv].npart>NumPartCut)
	      {
	        j += cut+1;
	      }
	    #endif
	   
	    k   = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];
	    id  = idv;
	    idv = k;
	    j++;        
	  }
		
		#ifdef BRANCH_SURVIVE
	    if(Gr[idv].npart>NumPartCut)
	    {
	      j += cut+1;
	    }
	  #endif
	
	  if(j<=cut)
	    return 1;
	  else 
	    return 0;
	
	}
	
	static void Delete_Branch(const TYPE_INT i, std::vector<std::vector<TYPE_INT> > &vec)
	{
	  TYPE_INT k, id, idv;
	
	  id = i;
	  idv = vec[id][0];
	  vec[id].clear();
	 
	  while(vec[idv].size()==2)
	  {
	    k = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];
	
	    vec[idv].clear();
	
	    id = idv;
	    idv = k;
	  }

	  if(vec[idv].size()==1)
	    vec[idv].clear();

	  return;
	}
	
#ifdef BRANCH_SURVIVE
	extern void Podado(const TYPE_INT level, std::vector<std::vector<TYPE_INT> > &vec, const TYPE_INT NumPartCut)
#else
	extern void Podado(const TYPE_INT level, std::vector<std::vector<TYPE_INT> > &vec)
#endif
	{
	  TYPE_INT i, k, cut, N_threads, itera;
	  std::vector<TYPE_INT> aux;
	
	  #ifdef NTHREADS
	  omp_set_dynamic(0);
	  omp_set_num_threads(NTHREADS);
	  #endif
	
	  N_threads = NTHREADS;
	
	  fprintf(stdout,"pruned level %d\n",level);
	  fflush(stdout);
	
	  cut = 1;
	
	  do{
	
	    itera = 0;
	
#ifdef BRANCH_SURVIVE
	    #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	    shared(cp,Gr,vec,cut,itera,NumPartCut)
#else
	    #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	    shared(cp,Gr,vec,cut,itera)
#endif
	    for(i=0;i<cp.ngrup;i++)
	    {
	
	      if(vec[i].size()==1)
	      {
	        
	      #ifdef BRANCH_SURVIVE
	        if(Gr[i].npart>NumPartCut) continue;
	      #endif
	
        #ifdef BRANCH_SURVIVE
	        if(DFS(i,vec,cut,NumPartCut))
        #else
	        if(DFS(i,vec,cut))
        #endif
	        {
	          Delete_Branch(i,vec);
	          itera = 1;
	        }
	  
	      }
	
	    }
	
	    if(itera==0) break;
	
	    #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	    private(aux,k) shared(cp,vec) 
	    for(i=0;i<cp.ngrup;i++)
	    {
	      if(vec[i].size()>2)
	      {
	        for(k=0;k<vec[i].size();k++)
	        {
	          if(!vec[vec[i][k]].empty())
	            aux.push_back(vec[i][k]);
	        }
	
	        swap(aux,vec[i]);
	        aux.clear();
	
	      }
	    }
	
	    if(cut==level) break;
	
	    cut++;
	
	  }while(1);
	
	  fprintf(stdout,"Termina el podado\n");
	  fflush(stdout);
	
	  return;
	  
	}

#endif

static void DLU(const TYPE_INT id, const TYPE_INT pre, std::vector<std::vector<TYPE_INT> > &adj, TYPE_INT *Padre, TYPE_INT *Rank)
{
  TYPE_INT k, idv;

  Padre[id] = pre;
  Rank[id]  = adj[id].size();

  for(k=0;k<Rank[id];k++)
  {
    idv = adj[id][k];
    if(idv != pre)
	  {
			DLU(idv,id,adj,Padre,Rank);
	  }
  }

  adj[id].clear();

  return;

}

static void DL(std::vector<std::pair<TYPE_INT,TYPE_INT> > &vec_sort, std::vector<std::vector<TYPE_INT> > &vec, \
TYPE_INT * __restrict__ Padre, TYPE_INT * __restrict__ Rank)
{
  TYPE_INT i;

  for(i=0;i<cp.ngrup;i++)
  {
    Padre[i] = cp.ngrup;
    Rank[i]  = 0;

    if(vec[i].size()>0)
      vec_sort.push_back(std::make_pair(Gr[i].npart,i));
  }

  sort(vec_sort.begin(),vec_sort.end());

  for(i=vec_sort.size();i>0;i--)
  {
    TYPE_INT k = vec_sort[i-1].second;

    if(Padre[k]==cp.ngrup)
		{
      DLU(k,Padre[k],vec,Padre,Rank);
  	}
  }

  return;
}

extern void Selection(const TYPE_INT NumPartCut, std::vector<std::vector<TYPE_INT> > &adjacency_list)
{
  TYPE_INT i,j,k,id;
  TYPE_INT *Padre, *Rank;
  std::vector<TYPE_INT> aux;
  std::vector<std::vector<TYPE_INT> > segmentos;
  std::vector<std::pair<TYPE_INT,TYPE_INT> > vec_orden;
	
	Padre = (TYPE_INT *) malloc(cp.ngrup*sizeof(TYPE_INT));
  Rank =  (TYPE_INT *) malloc(cp.ngrup*sizeof(TYPE_INT));

	DL(vec_orden,adjacency_list,Padre,Rank);

  j = 0;
  while(!vec_orden.empty())
  {
    i = vec_orden.back().second;

    if(Padre[i]>=cp.ngrup || Rank[i] >= cp.ngrup)
		{
    	vec_orden.pop_back();
			continue;
		}

    aux.push_back(i);
    Rank[i] += cp.ngrup;
    id = Padre[i];

    while(id<cp.ngrup)
    {
      if(Rank[id]>=cp.ngrup)
      {
        aux.push_back(id);
    		Rank[id] += cp.ngrup;
        break;
      }

      if(Gr[id].npart>NumPartCut)
      {
        aux.push_back(id);
        Rank[id] += cp.ngrup;
        segmentos.push_back(aux);
        aux.clear();
        j++;

				if(Padre[id]>=cp.ngrup)
					break;
      }

      aux.push_back(id);
      Rank[id] += cp.ngrup;
      id = Padre[id];               
    }

		if(aux.size() > 1) 
		{
			segmentos.push_back(aux);
			j++;
		}

    aux.clear();
    vec_orden.pop_back();
  }

  assert(segmentos.size()==j);

  //////for(i=0;i<j;i++)
  //////{
	//////  if(Gr[segmentos[i][0]].npart>Gr[segmentos[i].back()].npart)
	//////	{
	//////  	id = segmentos[i].size()-1;
	//////		aux.push_back(segmentos[i][id]);
	//////  	for(k=segmentos[i].size()-1;k>1;k--)
	//////  	{
	//////  	  if(Gr[segmentos[i][id]].npart<Gr[segmentos[i][k-1]].npart)
	//////			{
	//////				aux.push_back(segmentos[i][k-1]);
	//////				segmentos.push_back(aux);
	//////				aux.clear();
	//////				id = k-1;
	//////				aux.push_back(segmentos[i][id]);
	//////			}else{
	//////				aux.push_back(segmentos[i][k-1]);
	//////			}
	////// 		}
	//////		aux.push_back(segmentos[i][0]);

	//////		if(aux.size() != segmentos[i].size())
	//////			segmentos[i].swap(aux);
	//////		aux.clear();
	//////		
	//////	}else{

	//////  	id = 0;
	//////		aux.push_back(segmentos[i][id]);
	//////  	for(k=1;k<segmentos[i].size()-1;k++)
	//////  	{
	//////  	  if(Gr[segmentos[i][id]].npart<Gr[segmentos[i][k]].npart)
	//////			{
	//////				aux.push_back(segmentos[i][k]);
	//////				segmentos.push_back(aux);
	//////				aux.clear();
	//////				id = k;
	//////				aux.push_back(segmentos[i][id]);
	//////			}else{
	//////				aux.push_back(segmentos[i][k]);
	//////			}
	////// 		}
	//////		aux.push_back(segmentos[i][segmentos[i].size()-1]);

	//////		if(aux.size() != segmentos[i].size())
	//////			segmentos[i].swap(aux);
	//////		aux.clear();

	//////	}
	//////}

	fprintf(stdout,"Number Seg %lu\n", segmentos.size());
	cp.nseg = segmentos.size();

  #ifdef SORT

	  BLUE("********** Importante ***********\n");
  	sprintf(message,"Ordena el nodo de la derecha es más grande\n");RED(message);
	  BLUE("*********************************\n");

  #endif

  Seg  = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  for(i=0;i<cp.nseg;i++)
  {
    #ifdef SORT
    if(Gr[segmentos[i][0]].npart>Gr[segmentos[i].back()].npart)
      reverse(segmentos[i].begin(),segmentos[i].end());
    #endif

    aux = segmentos[i];
    
    k = aux.back();

		Seg[i].size = Seg[i].size_raw = aux.size();
    Seg[i].flag = 0;
    if(Gr[aux[0]].npart>NumPartCut) Seg[i].flag++;
    if(Gr[k].npart>NumPartCut)      Seg[i].flag++;

    Seg[i].Mass[0] = cp.Mpart*(TYPE_REAL)Gr[aux[0]].npart;
    Seg[i].Mass[1] = cp.Mpart*(TYPE_REAL)Gr[k].npart; 
    Seg[i].razon = Seg[i].Mass[1]<Seg[i].Mass[0] ? Seg[i].Mass[1]/Seg[i].Mass[0] : Seg[i].Mass[0]/Seg[i].Mass[1];

#ifdef WITH_EXTREMES
    Seg[i].Rvir[0] = Gr[aux[0]].rvir;
    Seg[i].Rvir[1] = Gr[k].rvir; 
#endif
 
#ifdef STORE_VELOCITIES
	  for(j=0;j<3;j++)
    {
      Seg[i].Vnodos[j]   = Gr[aux[0]].vcm[j];
      Seg[i].Vnodos[j+3] = Gr[k].vcm[j];
    }
#endif

		Seg[i].Ids_list     = (TYPE_INT  *) malloc(Seg[i].size_raw*sizeof(TYPE_INT));
		Seg[i].Pos_list     = (TYPE_REAL *) malloc(3*Seg[i].size_raw*sizeof(TYPE_REAL));
		Seg[i].Pos_list_raw = (TYPE_REAL *) malloc(3*Seg[i].size_raw*sizeof(TYPE_REAL));

    for(k=0;k<Seg[i].size_raw;k++)
    {
      Seg[i].Ids_list[k] = aux[k];
			for(j=0;j<3;j++)
			{
	      Seg[i].Pos_list[3*k+j] = \
				Seg[i].Pos_list_raw[3*k+j] = \
				Gr[aux[k]].pcm[j];
    	}
    }
  }

  while(!segmentos.empty())
    segmentos.pop_back();

	free(Padre);
  free(Rank);

  return;
}
