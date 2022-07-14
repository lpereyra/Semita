#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "variables.hh"
#include "properties.hh"
#include "colors.hh"

extern void init_variables(int argc, char **argv)
{
  FILE *pfin;
  char filename[200];
  TYPE_INT i;

  RED("Initializing variables...\n");

  sprintf(filename,"%s",argv[1]);
  if(!(pfin=fopen(filename,"r")))
  {
    sprintf(message,"can't open file `%s` \n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&snap.nfiles))
  {
    sprintf(message,"can't read file `%s`\nneed # of snapshots\n",filename);RED(message);
    exit(0);
  }
    
  if(!fscanf(pfin,"%s  \n",snap.root))
  {
    sprintf(message,"can't read file `%s`\nneed snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.name))
  {
    sprintf(message,"can't read file `%s`\nneed snapname\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf \n",&cp.soft))
  {
    sprintf(message,"can't read file `%s`\nneed softening of simulation\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%u  \n",&nfrac))
  {
    sprintf(message,"can't read file `%s`\nneed identification steps\n",filename);RED(message);
    exit(0);
  }

  fof = (TYPE_REAL *) malloc(nfrac*sizeof(TYPE_REAL));

  for(i=0;i<nfrac;i++)
  {
    #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf  \n",&fof[i]))
    #else
    if(!fscanf(pfin,"%f   \n",&fof[i]))
    #endif
    {
      sprintf(message,"can't read file `%s`\nneed %d identification step\n",filename,i);RED(message);
      exit(0);
    }
  }

  fclose(pfin);

  /////////////////////////////////////////////////////////////////////////
  char *p = snap.name;
  while (*p) 
  { 
    if(isdigit(*p)) 
    { // Upon finding a digit, ...
        snap.num = strtol(p, &p, 10); // Read a number, ...
    } else { // Otherwise, move on to the next character.
        p++;
    }
  }
  ///////////////////////////////////////////////////////////////////////

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"Snapname Num:            %d\n",snap.num);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
#ifdef READ_LEVELS
  BLUE("************* Read Levels ************\n");
  for(i=0;i<nfrac;i++)
  {
    sprintf(message,"%d overdensity %.2f\n",i,fof[i]);BLUE(message);
    fof[i] = cbrt(1./(1.+fof[i]));
  }
#endif

  BLUE("********** Makefile Options ***********\n");
  #ifdef DEBUG
  BLUE("  DEBUG\n");
  #endif
  #ifdef PERIODIC
  BLUE("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  BLUE("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  BLUE("  LONGIDS\n");
  #endif
  #ifdef POSFACTOR
  sprintf(message,"  POSFACTOR = %f\n",POSFACTOR);BLUE(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef LOCK
  BLUE("  USING LOCKS\n");
  #endif

  GREEN("END\n");
}

extern void write_properties(FILE *pfout, struct grup_data Prop)
{
	TYPE_INT   jj;

	fwrite(&Prop.save,sizeof(TYPE_INT),1,pfout);
	fwrite(&Prop.id,sizeof(TYPE_INT),1,pfout);
	fwrite(&Prop.npart,sizeof(TYPE_INT),1,pfout);
  fwrite(&Prop.pcm,sizeof(TYPE_REAL),3,pfout);
#ifdef STORE_VELOCITIES  
  fwrite(&Prop.vcm,sizeof(TYPE_REAL),3,pfout);
#ifdef COMPUTE_EP
  fwrite(&Prop.mostbound,sizeof(TYPE_REAL),3,pfout);
#endif
  fwrite(&Prop.sig,sizeof(TYPE_REAL),3,pfout);
  fwrite(&Prop.L,sizeof(TYPE_REAL),3,pfout);
#ifdef COMPUTE_EP
  fwrite(&Prop.lambda,sizeof(TYPE_REAL),1,pfout);
#endif
  fwrite(&Prop.m200,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.r200,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.v200,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.mvir,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.rvir,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.vvir,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.vmax,sizeof(TYPE_REAL),1,pfout);
#ifdef COMPUTE_EP
  fwrite(&Prop.Ep,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.Ec,sizeof(TYPE_REAL),1,pfout);
#endif  
#endif  
  fwrite(&Prop.aa,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.bb,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.cc,sizeof(TYPE_REAL),1,pfout);
  for(jj=0;jj<3;jj++)
    fwrite(&Prop.evec[jj],sizeof(TYPE_REAL),3,pfout);
#ifdef STORE_VELOCITIES  
  fwrite(&Prop.aa_vel,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.bb_vel,sizeof(TYPE_REAL),1,pfout);
  fwrite(&Prop.cc_vel,sizeof(TYPE_REAL),1,pfout);
  for(jj=0;jj<3;jj++)
    fwrite(&Prop.evec_vel[jj],sizeof(TYPE_REAL),3,pfout);
#endif  

  return;
}

#ifdef FILE_ASCII

  extern void write_properties_ascii(FILE *pfout, struct grup_data Prop)
  {
	  TYPE_INT   jj;

    fprintf(pfout,"%d %f %f %f ",
         Prop.npart,Prop.pcm[0],Prop.pcm[1],Prop.pcm[2]);
  #ifdef STORE_VELOCITIES  
    fprintf(pfout,"%f %f %f ",
         Prop.vcm[0],Prop.vcm[1],Prop.vcm[2]); 
    fprintf(pfout,"%f %f %f %f %f %f ",
         Prop.sig[0],Prop.sig[1],Prop.sig[2],Prop.L[0],Prop.L[1],Prop.L[2]);
  #ifdef COMPUTE_EP
    fprintf(pfout,"%f ",
         Prop.lambda);
  #endif
    fprintf(pfout,"%f %f %f %f %f %f %f ",
         Prop.m200,Prop.r200,Prop.v200,Prop.mvir,Prop.rvir,Prop.vvir,Prop.vmax);
  #ifdef COMPUTE_EP
    fprintf(pfout,"%f %f ",
         Prop.Ep,Prop.Ec);
  #endif
  #endif
    fprintf(pfout,"%f %f %f",
         Prop.aa,Prop.bb,Prop.cc);
    for(jj=0;jj<3;jj++)
      fprintf(pfout," %f %f %f",Prop.evec[jj][0],Prop.evec[jj][1],Prop.evec[jj][2]);
  #ifdef STORE_VELOCITIES  
    fprintf(pfout," %f %f %f",
         Prop.aa_vel,Prop.bb_vel,Prop.cc_vel);
    for(jj=0;jj<3;jj++)
      fprintf(pfout," %f %f %f",Prop.evec_vel[jj][0],Prop.evec_vel[jj][1],Prop.evec_vel[jj][2]);
  #endif
    fprintf(pfout,"\n");
    fflush(pfout);
  }

#endif

static void set_name(const char * prefix, char * name, const TYPE_INT NNN, const TYPE_REAL * fof, const TYPE_INT flag)
{
  sprintf(name,"%.2d_%.4d_%s",snap.num,NNN,prefix);

#ifdef PARTICLES
  if(flag == 1)
    sprintf(name,"%s_%.2f_%.2f_%.2d_RVIR.bin",name,fof[0],fof[1],(TYPE_INT)RVIR_FACTOR);
  else
#endif
    sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

	return;
}

extern void write_seg(const TYPE_INT NNN, const TYPE_REAL * fof)
{
  TYPE_INT i,j;
  char filename[200];
  FILE *pfout, *pfpropiedades;
  FILE *pfout_raw, *pfpropiedades_raw;

  fprintf(stdout,"WRITE...\n");
  fflush(stdout);

	set_name("smooth_segmentos",filename,NNN,fof,0);
  pfout=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(TYPE_INT),1,pfout);

	set_name("smooth_propiedades",filename,NNN,fof,1);
  pfpropiedades=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(TYPE_INT),1,pfpropiedades);

	set_name("raw_segmentos",filename,NNN,fof,0);
  pfout_raw=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(TYPE_INT),1,pfout_raw);

	set_name("raw_propiedades",filename,NNN,fof,0);
  pfpropiedades_raw=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(TYPE_INT),1,pfpropiedades_raw);

  for(i=0;i<cp.nseg;i++)
  {
    fwrite(&Seg[i].flag,sizeof(TYPE_INT),1,pfpropiedades);
    fwrite(&Seg[i].size,sizeof(TYPE_INT),1,pfpropiedades);
    fwrite(&Seg[i].Mass[0],sizeof(TYPE_REAL),2,pfpropiedades); 
#ifdef WITH_EXTREMES
    fwrite(&Seg[i].Rvir[0],sizeof(TYPE_REAL),2,pfpropiedades); 
#endif
#ifdef STORE_VELOCITIES
    fwrite(&Seg[i].Vnodos[0],sizeof(TYPE_REAL),6,pfpropiedades); 
#endif
    fwrite(&Seg[i].razon,sizeof(TYPE_REAL),1,pfpropiedades); 
    fwrite(&Seg[i].len,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].cur,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].rms,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].vol,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].mass_part,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].rho,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].mu,sizeof(TYPE_REAL),1,pfpropiedades);
#ifdef STORE_VELOCITIES
    fwrite(&Seg[i].sigma_per,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].r_pos,sizeof(TYPE_REAL),1,pfpropiedades);
    fwrite(&Seg[i].id_pos,sizeof(int),1,pfpropiedades);
#endif

    fwrite(&Seg[i].size,sizeof(TYPE_INT),1,pfout);
    for(j=0;j<Seg[i].size;j++)
    {
      fwrite(&Seg[i].Pos_list[3*j],  sizeof(TYPE_REAL),1,pfout); 
      fwrite(&Seg[i].Pos_list[3*j+1],sizeof(TYPE_REAL),1,pfout); 
      fwrite(&Seg[i].Pos_list[3*j+2],sizeof(TYPE_REAL),1,pfout); 
    }
		
		//////////////////////////////////////////////////////////////////
	
    fwrite(&Seg[i].flag,sizeof(TYPE_INT),1,pfpropiedades_raw);
    fwrite(&Seg[i].size_raw,sizeof(TYPE_INT),1,pfpropiedades_raw);
    fwrite(&Seg[i].Mass[0],sizeof(TYPE_REAL),2,pfpropiedades_raw); 
#ifdef STORE_VELOCITIES
    fwrite(&Seg[i].Vnodos[0],sizeof(TYPE_REAL),6,pfpropiedades_raw); 
#endif
    fwrite(&Seg[i].razon,sizeof(TYPE_REAL),1,pfpropiedades_raw); 
    fwrite(&Seg[i].len_raw,sizeof(TYPE_REAL),1,pfpropiedades_raw);
    fwrite(&Seg[i].cur_raw,sizeof(TYPE_REAL),1,pfpropiedades_raw);
    fwrite(&Seg[i].rms_raw,sizeof(TYPE_REAL),1,pfpropiedades_raw);

    fwrite(&Seg[i].size_raw,sizeof(TYPE_INT),1,pfout_raw);
    for(j=0;j<Seg[i].size_raw;j++)
      fwrite(&Seg[i].Ids_list[j],sizeof(TYPE_INT),1,pfout_raw); 
  }

  fclose(pfout);
  fclose(pfpropiedades);
  fclose(pfout_raw);
  fclose(pfpropiedades_raw);

  return;

}
