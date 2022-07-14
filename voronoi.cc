#include <assert.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <functional>
#include <algorithm>
#include <vector>

#include "variables.hh"
#include "grid.hh"
#include "voronoi.hh"
#include "voro++/voro++.hh"

static TYPE_INT control(const TYPE_REAL rsph, const TYPE_REAL rsph2, const TYPE_INT idg, const TYPE_REAL Pos_cent[])
{
	
  TYPE_INT i, ilim;
  long ixc, iyc, izc, ibox;
  TYPE_REAL Posprima[3];
  TYPE_REAL dd;

  ilim = (long)(rsph*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox))+1;
  ixc  = (long)(Pos_cent[0]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(Pos_cent[1]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(Pos_cent[2]*(TYPE_REAL)grid.ngrid*(1.f/cp.lbox));

  #ifndef PERIODIC
    for(long ixx = ((ixc-ilim<0) ? 0 : ixc-ilim); ixx <= ((ixc+ilim >= grid.ngrid) ? grid.ngrid-1 : ixc+ilim); ixx++)
  #else
    for(long ixx = ixc-ilim; ixx <= ixc+ilim; ixx++)
  #endif
  {
    #ifndef PERIODIC
      for(long iyy = ((iyc-ilim<0) ? 0 : iyc-ilim); iyy <= ((iyc+ilim >= grid.ngrid) ? grid.ngrid-1 : iyc+ilim); iyy++)
    #else
      for(long iyy = iyc-ilim ; iyy <= iyc+ilim ; iyy++)
    #endif
    {
      #ifndef PERIODIC
        for(long izz = ((izc-ilim<0) ? 0 : izc-ilim); izz <= ((izc+ilim >= grid.ngrid) ? grid.ngrid-1 : izc+ilim); izz++)
      #else
        for(long izz = izc-ilim ; izz <= izc+ilim ; izz++)
      #endif
      {

      	#ifdef PERIODIC
          ibox = igrid(( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) ),\
                       ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ),\
                       ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) ),\
                       (long)grid.ngrid);
        #else
          ibox = igrid(ixx,iyy,izz,(long)grid.ngrid);
        #endif

        for(i = grid.icell[ibox]; i < grid.icell[ibox+1]; i++)
        {

          if(P[i].sub==idg)
          {
 
#ifdef COLUMN
            Posprima[0] = P.x[i] - Pos_cent[0];
            Posprima[1] = P.y[i] - Pos_cent[1];
            Posprima[2] = P.z[i] - Pos_cent[2];
#else
            Posprima[0] = P[i].pos[0] - Pos_cent[0];
            Posprima[1] = P[i].pos[1] - Pos_cent[1];
            Posprima[2] = P[i].pos[2] - Pos_cent[2];
#endif
            #ifdef PERIODIC
            if(Posprima[0] >  cp.lbox*0.5f) Posprima[0] = Posprima[0] - cp.lbox;
            if(Posprima[1] >  cp.lbox*0.5f) Posprima[1] = Posprima[1] - cp.lbox;
            if(Posprima[2] >  cp.lbox*0.5f) Posprima[2] = Posprima[2] - cp.lbox;
            if(Posprima[0] < -cp.lbox*0.5f) Posprima[0] = Posprima[0] + cp.lbox;
            if(Posprima[1] < -cp.lbox*0.5f) Posprima[1] = Posprima[1] + cp.lbox;
            if(Posprima[2] < -cp.lbox*0.5f) Posprima[2] = Posprima[2] + cp.lbox;
            #endif

            dd = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

            if(dd<rsph2)
            {
              return 1;
            }
          }
        } /*end particles loop*/
      } /*end izz*/
    } /*end iyy*/
  } /*end ixx*/  

  return 0;

}

static TYPE_INT inside_fof(const TYPE_INT id, const TYPE_REAL r0, const TYPE_REAL r0_2, const TYPE_REAL r, const TYPE_REAL vdir[])
{
   TYPE_INT  c;
   TYPE_INT  itera = (TYPE_INT)(r/r0);
   TYPE_REAL Pos_cent[3];
   TYPE_REAL rsep = r0;

   //fprintf(stdout,"%d\n",itera);

   Pos_cent[0] = Gr[id].pcm[0];
   Pos_cent[1] = Gr[id].pcm[1];
   Pos_cent[2] = Gr[id].pcm[2];

   if(itera==0)
   {     
      Pos_cent[0] += (0.5*r)*vdir[0];
      Pos_cent[1] += (0.5*r)*vdir[1];
      Pos_cent[2] += (0.5*r)*vdir[2];

      #ifdef PERIODIC
      Pos_cent[0] = Pos_cent[0]<0.0f    ? Pos_cent[0]+cp.lbox : Pos_cent[0];
      Pos_cent[1] = Pos_cent[1]<0.0f    ? Pos_cent[1]+cp.lbox : Pos_cent[1];
      Pos_cent[2] = Pos_cent[2]<0.0f    ? Pos_cent[2]+cp.lbox : Pos_cent[2];

      Pos_cent[0] = Pos_cent[0]>cp.lbox ? Pos_cent[0]-cp.lbox : Pos_cent[0];
      Pos_cent[1] = Pos_cent[1]>cp.lbox ? Pos_cent[1]-cp.lbox : Pos_cent[1];
      Pos_cent[2] = Pos_cent[2]>cp.lbox ? Pos_cent[2]-cp.lbox : Pos_cent[2];
      #endif

      //if(!control(0.5*r,r0_2,rsph,Gr[id].save,vdir,Pos_cent))
      if(control(r0, r0_2, Gr[id].save, Pos_cent) == 0)
        return 0; // primer punto evaluado
   }

   Pos_cent[0] += (0.5f*rsep)*vdir[0];
   Pos_cent[1] += (0.5f*rsep)*vdir[1];
   Pos_cent[2] += (0.5f*rsep)*vdir[2];

   #ifdef PERIODIC
   Pos_cent[0] = Pos_cent[0]<0.0f    ? Pos_cent[0]+cp.lbox : Pos_cent[0];
   Pos_cent[1] = Pos_cent[1]<0.0f    ? Pos_cent[1]+cp.lbox : Pos_cent[1];
   Pos_cent[2] = Pos_cent[2]<0.0f    ? Pos_cent[2]+cp.lbox : Pos_cent[2];

   Pos_cent[0] = Pos_cent[0]>cp.lbox ? Pos_cent[0]-cp.lbox : Pos_cent[0];
   Pos_cent[1] = Pos_cent[1]>cp.lbox ? Pos_cent[1]-cp.lbox : Pos_cent[1];
   Pos_cent[2] = Pos_cent[2]>cp.lbox ? Pos_cent[2]-cp.lbox : Pos_cent[2];
   #endif

   if(control(r0, r0_2, Gr[id].save, Pos_cent) == 0)
     return 0; // primer punto evaluado

   for(c=1;c<itera;c++)
   {

     Pos_cent[0] += rsep*vdir[0];
     Pos_cent[1] += rsep*vdir[1];
     Pos_cent[2] += rsep*vdir[2];

     #ifdef PERIODIC
     Pos_cent[0] = Pos_cent[0]<0.0f    ? Pos_cent[0]+cp.lbox : Pos_cent[0];
     Pos_cent[1] = Pos_cent[1]<0.0f    ? Pos_cent[1]+cp.lbox : Pos_cent[1];
     Pos_cent[2] = Pos_cent[2]<0.0f    ? Pos_cent[2]+cp.lbox : Pos_cent[2];
    
     Pos_cent[0] = Pos_cent[0]>cp.lbox ? Pos_cent[0]-cp.lbox : Pos_cent[0];
     Pos_cent[1] = Pos_cent[1]>cp.lbox ? Pos_cent[1]-cp.lbox : Pos_cent[1];
     Pos_cent[2] = Pos_cent[2]>cp.lbox ? Pos_cent[2]-cp.lbox : Pos_cent[2];
     #endif

     if(control(r0, r0_2, Gr[id].save, Pos_cent) == 0)
       return 0; // primer punto evaluado

   }

   Pos_cent[0] += (0.5*rsep)*vdir[0];
   Pos_cent[1] += (0.5*rsep)*vdir[1];
   Pos_cent[2] += (0.5*rsep)*vdir[2];

   rsep = r - (TYPE_REAL)itera*rsep; // cilindro pequeño centrado en la diferencia

   Pos_cent[0] += (0.5*rsep)*vdir[0];
   Pos_cent[1] += (0.5*rsep)*vdir[1];
   Pos_cent[2] += (0.5*rsep)*vdir[2];

   #ifdef PERIODIC
   Pos_cent[0] = Pos_cent[0]<0.0f    ? Pos_cent[0]+cp.lbox : Pos_cent[0];
   Pos_cent[1] = Pos_cent[1]<0.0f    ? Pos_cent[1]+cp.lbox : Pos_cent[1];
   Pos_cent[2] = Pos_cent[2]<0.0f    ? Pos_cent[2]+cp.lbox : Pos_cent[2];
   
   Pos_cent[0] = Pos_cent[0]>cp.lbox ? Pos_cent[0]-cp.lbox : Pos_cent[0];
   Pos_cent[1] = Pos_cent[1]>cp.lbox ? Pos_cent[1]-cp.lbox : Pos_cent[1];
   Pos_cent[2] = Pos_cent[2]>cp.lbox ? Pos_cent[2]-cp.lbox : Pos_cent[2];
   #endif
  
   if(control(r0, r0_2, Gr[id].save, Pos_cent) == 0)
      return 0; // ultimo punto evaluado

   return 1; 
}

static void init_edges(std::vector<std::pair<TYPE_INT,TYPE_INT> > &pares)
{

  TYPE_INT i, Ngrid;
  const double gap = 1.0e-10; // uso un pequeño gap
  double x_min, y_min, z_min;	  
  double x_max, y_max, z_max;	  
  bool xbool, ybool, zbool;
  voro::voronoicell_neighbor cell;

  Ngrid = (TYPE_INT)pow((float)cp.ngrup/5.0,1.0/3.0);
  #ifdef PERIODIC
    xbool = ybool = zbool = true;
  #else
    xbool = ybool = zbool = false;
  #endif

  x_min = y_min = z_min = 0.0-gap;	  
  x_max = y_max = z_max = cp.lbox+gap;	  

  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,Ngrid,Ngrid,Ngrid,xbool,ybool,zbool,8);

  for(i=0;i<cp.ngrup; i++)
  {
    con.put(i,Gr[i].pcm[0],Gr[i].pcm[1],Gr[i].pcm[2]);	  
  }
  assert(cp.ngrup==(TYPE_INT)con.total_particles());

  voro::c_loop_all clo(con);

  if(clo.start()) do if(con.compute_cell(cell,clo))
  {
    std::vector<int>  vec;
    int id = clo.pid();
    cell.neighbors(vec);

    for(i=0; i<vec.size(); i++)
    {
      int idv = vec[i];
    
      #ifndef PERIODIC
      if(idv<0) continue;
      #endif

      if(Gr[id].save != Gr[idv].save) continue;

      if(Gr[id].id > Gr[idv].id)
        pares.push_back(std::make_pair((TYPE_INT)idv,(TYPE_INT)id));

    }

    vec.clear();

  }while(clo.inc());

  con.clear();

  return;
}

extern void Voronoi_Grupos(const TYPE_REAL fof, std::vector<std::pair<TYPE_REAL,std::pair<TYPE_INT,TYPE_INT> > > &edges)
{

  TYPE_INT i, idv, id, Tid;
  TYPE_REAL r,r0,r0_2;
  TYPE_REAL vdir[3];
  std::vector<std::pair<TYPE_INT,TYPE_INT> > pares;
  std::vector<std::vector<std::pair<TYPE_REAL,std::pair<TYPE_INT,TYPE_INT> > > > lados(NTHREADS);
  #ifdef WRITE_WEIGHT
    char filename[200];
    TYPE_INT  *c;
    FILE **pfweight;
  #endif

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  ///////////////////////////////////////////

  init_edges(pares);
  fprintf(stdout,"Npares %lu\n",pares.size());

  ///////////////////////////////////////////
 
  //// Init Grid ////
  r0  = cbrt(1./(1.+fof));
  r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.f;
  grid.nobj  = cp.npart;
  grid.ngrid = (TYPE_INT)(1.1f*cp.lbox/r0);
  if(grid.ngrid > NGRIDMAX)
  {
    grid.ngrid = NGRIDMAX;
    fprintf(stdout,"Using NGRIDMAX = %lu, r0 = %f, lbox = %f\n",grid.ngrid,r0,cp.lbox);
  }else{
    fprintf(stdout,"Ngrid = %lu, r0 = %f, lbox = %f\n",grid.ngrid,r0,cp.lbox);
  }

  grid_init();
  grid_build();
  r0_2 = r0*r0; // Para usar el cuadrado
  r0 /= 100.; // Para usar el cuadrado
  ///////////////////////////////////////////

  #ifdef WRITE_WEIGHT
    c       = (TYPE_INT *)  malloc(NTHREADS*sizeof(TYPE_INT));
    pfweight = (FILE **) malloc(NTHREADS*sizeof(FILE));

    for(i=0;i<NTHREADS;i++)
    {
      c[i] = 0;
      sprintf(filename,"%.2d_weight.%.2d.bin",snap.num,i);
      pfweight[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(TYPE_INT),1,pfweight[i]);    
    }
  #endif

  #ifdef WRITE_WEIGHT

    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,r,vdir) \
    shared(r0,r0_2,cp,Gr,grid,pares,lados,pfweight,c,stdout)

  #else

    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,r,vdir) \
    shared(r0,r0_2,cp,Gr,grid,pares,lados,stdout)

  #endif
  for(i=0;i<pares.size();i++)
  {

    Tid = omp_get_thread_num();
    id  = pares[i].first;
    idv = pares[i].second;

    if(i%1000000==0) fprintf(stdout,"%u %u %.4f\n",Tid,i,(float)i/(float)pares.size());

    vdir[0] = Gr[idv].pcm[0] - Gr[id].pcm[0];
    vdir[1] = Gr[idv].pcm[1] - Gr[id].pcm[1];
    vdir[2] = Gr[idv].pcm[2] - Gr[id].pcm[2];

    #ifdef PERIODIC
    vdir[0] = vdir[0] >= cp.lbox*0.5 ? vdir[0]-cp.lbox : vdir[0];
    vdir[0] = vdir[0] < -cp.lbox*0.5 ? vdir[0]+cp.lbox : vdir[0];

    vdir[1] = vdir[1] >= cp.lbox*0.5 ? vdir[1]-cp.lbox : vdir[1];
    vdir[1] = vdir[1] < -cp.lbox*0.5 ? vdir[1]+cp.lbox : vdir[1];

    vdir[2] = vdir[2] >= cp.lbox*0.5 ? vdir[2]-cp.lbox : vdir[2];
    vdir[2] = vdir[2] < -cp.lbox*0.5 ? vdir[2]+cp.lbox : vdir[2];
    #endif

    r = sqrt(vdir[0]*vdir[0]+vdir[1]*vdir[1]+vdir[2]*vdir[2]);

    vdir[0] *= (1./r);
    vdir[1] *= (1./r);
    vdir[2] *= (1./r);

    if(inside_fof(id,r0,r0_2,r,vdir)==0) continue;
                  //    (id, r0, r0_2, const TYPE_REAL r, const TYPE_REAL vdir[])

    #ifdef WRITE_WEIGHT
    fwrite(&id,sizeof(TYPE_INT),1,pfweight[Tid]);
    fwrite(Gr[id].pcm,sizeof(TYPE_REAL),3,pfweight[Tid]);
    fwrite(&Gr[id].npart,sizeof(TYPE_INT),1,pfweight[Tid]);
    fwrite(&idv,sizeof(TYPE_INT),1,pfweight[Tid]);
    fwrite(Gr[idv].pcm,sizeof(TYPE_REAL),3,pfweight[Tid]);
    fwrite(&Gr[idv].npart,sizeof(TYPE_INT),1,pfweight[Tid]);
    fwrite(&r,sizeof(TYPE_REAL),1,pfweight[Tid]); // distancia
    #endif

    r = -((TYPE_REAL)Gr[id].npart*(TYPE_REAL)Gr[idv].npart)/(r*r);
    lados[Tid].push_back(std::make_pair((TYPE_REAL)r,pares[i]));

    #ifdef WRITE_WEIGHT
    fwrite(&r,sizeof(TYPE_REAL),1,pfweight[Tid]); // peso
    c[Tid]++;
    #endif

  }

  pares.clear();

  #ifdef WRITE_WEIGHT
  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfweight[i]);
    fwrite(&c[i],sizeof(int),1,pfweight[i]);
    fclose(pfweight[i]);
  }

  free(c);
  #endif

  for(i=0;i<NTHREADS;i++)
  {
    edges.insert(edges.end(),lados[i].begin(),lados[i].end());
    lados[i].clear();
  }

	fprintf(stdout,"%lu NumEdges\n",edges.size());
  fflush(stdout);

  grid_free();
  
  return;
}
