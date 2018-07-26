///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//          Read-write routines for the periodic-box mode            //
/*********************************************************************/
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"
#include "hdf5.h"

typedef struct {
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
  // fills to 256 Bytes
} gad_header;

typedef struct {
  char label[4];
  int size;
} gad_title;

static float wrap_double(float x)
{
  //////
  // Returns x mod(l_box)
  if(x<0) 
    return wrap_double(x+l_box);
  else if(x>=l_box) 
    return wrap_double(x-l_box);
  else 
    return x;
}

static size_t my_fread(void *p,size_t size,
		       size_t nmemb,FILE *stream)
{
  //////
  // Self-checked binary reading routine
  size_t nread;

  if((nread=fread(p,size,nmemb,stream))!=nmemb) {
    fprintf(stderr,"CUTE: error reading binary file \n");
    exit(1);
  }
  return nread;
}

void make_CF(unsigned long long DD[],int nD,
	     double corr[],double ercorr[])
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double *edd;
  double rho_av=nD/(l_box*l_box*l_box);
  printf("rho_av=%f\n", rho_av); 
  int ii;

  edd=(double *)malloc(sizeof(double)*NB_R);
  if(edd==NULL)
    error_mem_out();
  
#ifndef _LOGBIN
  DD[0]-=nD; //Substract diagonal
#endif //_LOGBIN
  for(ii=0;ii<NB_R;ii++)
    edd[ii]=1./sqrt((double)DD[ii]);

  for(ii=0;ii<NB_R;ii++) {
    if(DD[ii]==0) {
      corr[ii]=0;
      ercorr[ii]=0;
    }
    else {
      double r0,r1,vr,rho_r;
#ifdef _LOGBIN
      r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
      r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
      r0=ii/(I_R_MAX*NB_R);
      r1=(ii+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
      vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      rho_r=DD[ii]/(nD*vr);
      corr[ii]=rho_r/rho_av-1;
      ercorr[ii]=(1+corr[ii])*edd[ii];
    }
  }

  free(edd);
}

static void check_params(void)
{
  //////
  // Check all parameters are there and are sensible

  //vital files
  if(!strcmp(fnameData,"default")) {
    fprintf(stderr,"CUTE: Data catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameOut,"default")) {
    fprintf(stderr,"CUTE: Output filename was not provided \n");
    exit(1);
  }
  if(l_box<0) {
    fprintf(stderr,"CUTE: Box size was not provided \n");
    exit(1);
  }
  //input format
  if((input_format!=0)&&(input_format!=1)&&(input_format!=2)&&(input_format!=5)) {
    fprintf(stderr,"CUTE: wrong input format. Using standard ASCII file \n");
    input_format=0;
  }
  if((input_format>0)&&(n_objects!=-1)) {
    fprintf(stderr,"CUTE: can't select #objects for GADGET input format.");
    fprintf(stderr," Reading all objects \n");
    n_objects=-1;
  }
  //numbers of objects from data and random
  if((n_objects!=-1)&&(n_objects<0)) {
    fprintf(stderr,"CUTE: Wrong wumber of lines from the data file. ");
    fprintf(stderr,"CUTE: All objects in the catalog will be read.");
    n_objects=-1;
  }
  if((use_pm)&&(use_tree)) {
    fprintf(stderr,"CUTE: Tree and PM algorithms CAN'T be combined.");
    fprintf(stderr," Using tree \n");
    use_pm=0;
    use_tree=1;
  }

  if(use_tree<0) {
    fprintf(stderr,"CUTE: No Tree option was provided \n");
    exit(1);
  }
  else if(use_tree) {
    if(max_tree_order<0) {
      fprintf(stderr,"CUTE: Maximum tree order was not provided. Using 7 \n");
      max_tree_order=7;
    }
    if(max_tree_nparts<0) {
      fprintf(stderr,"CUTE: Maximum #particles per leaf was not provided.");
      fprintf(stderr," Using 10\n");
      max_tree_nparts=10;
    }
#ifndef _LOGBIN //Check resolution
    double binsize=1/(I_R_MAX*NB_R);
    double leafsize=l_box/pow(2,max_tree_order);

    if(binsize<=leafsize) {
      fprintf(stderr,"CUTE: Warning! binsize is smaller than cell size (%.3lf < %.3lf). ",
	      binsize,leafsize);
      fprintf(stderr," Tree may be suboptimal \n");
    }
#endif //_LOGBIN
  }

  if((use_pm!=0)&&(use_pm!=1)) {
    fprintf(stderr,"CUTE: No PM option was provided \n");
    exit(1);
  }
  if(use_pm) {
    if(n_grid<=0) {
      fprintf(stderr,"CUTE: No PM grid size was provided \n");
      exit(1);
    }

    //Check resolution
    double cellsize=l_box/n_grid;
#ifdef _LOGBIN
    fprintf(stderr,"CUTE: Warning! Logarithmic binning with PM algorithm. ");
    fprintf(stderr,"Scales below cellsize (%.3lf) won't be correctly calculated \n",cellsize);
#else //_LOGBIN
    double binsize=1/(I_R_MAX*NB_R);
    if(binsize<=cellsize) {
      fprintf(stderr,"CUTE: Warning! binsize is smaller than cell size (%.3lf < %.3lf). ",
	      binsize,cellsize);
      fprintf(stderr," Using PM is not recommended \n");
    }
#endif //_LOGBIN
  }
}

void read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  
  printf("*** Reading run parameters \n");
  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"data_filename="))
      sprintf(fnameData,"%s",s2);
    else if(!strcmp(s1,"num_lines=")) {
      if(!strcmp(s2,"all"))
	n_objects=-1;
      else
	n_objects=atoi(s2);
    }
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"box_size=")) {
      l_box=atof(s2);
      l_box_half=0.5*l_box;
    }
    else if(!strcmp(s1,"use_tree="))
      use_tree=atoi(s2);
    else if(!strcmp(s1,"max_tree_order="))
      max_tree_order=atoi(s2);
    else if(!strcmp(s1,"max_tree_nparts="))
      max_tree_nparts=atoi(s2);
    else if(!strcmp(s1,"use_pm="))
      use_pm=atoi(s2);
    else if(!strcmp(s1,"n_grid_side="))
      n_grid=atoi(s2);
    else
      fprintf(stderr,"CUTE: Unknown parameter %s\n",s1);
  }
  fclose(fi);
  
  check_params();
  
  printf("\n");
}

static void gad_check_block(int b1,int b2)
{
  //////
  // Checks that a block from a snapshot is
  // consistent from its begin/end values
  if(b1!=b2) {
    fprintf(stderr,"CUTE: Corrupted block!\n");
    exit(1);
  }
}

static int gad_seek_block(FILE *snap,char name[])
{
  //////
  // Seeks block from title
  gad_title tit;
  int block1,block2;

  rewind(snap);

  while(1>0) {
    if(!(fread(&block1,sizeof(int),1,snap))||
       feof(snap)||ferror(snap)) {
      fprintf(stderr,"CUTE: Block %s not found!!\n",name);
      exit(1);
    }
    my_fread(&tit,sizeof(gad_title),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);
    if(strncmp(tit.label,name,3)!=0)
      fseek(snap,tit.size,SEEK_CUR);
    else
      break;
  }

  return 0;
}

static Catalog read_ascii(char *fname,lint *np)
{ 
  //////
  // Reads catalog from ascii file with
  // default format.
  Catalog cat;
  lint ii,n_lin;
  FILE *fd;

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    n_lin=linecount(fd);
  else
    n_lin=n_objects;
  rewind(fd);
      
#ifdef _VERBOSE
  printf("  %ld objects will be read \n",(long)n_lin);
#endif
      
  //Allocate catalog memory
  *np=n_lin;
  cat.np=n_lin;
  cat.pos=(float *)malloc(3*cat.np*sizeof(float));
  if(cat.pos==NULL)
    error_mem_out();
  
  rewind(fd);
  //Read galaxies in mask
  for(ii=0;ii<n_lin;ii++) {
    float xx,yy,zz;
    char s0[1024];
    int sr;
    if(fgets(s0,sizeof(s0),fd)==NULL)
      error_read_line(fname,ii+1);
    sr=sscanf(s0,"%f %f %f",&xx,&yy,&zz);
    if(sr!=3)
      error_read_line(fname,ii+1);
    cat.pos[3*ii]=xx;
    cat.pos[3*ii+1]=yy;
    cat.pos[3*ii+2]=zz;
  }
  fclose(fd);

  return cat;
}

static int check_num_files(char *prefix)
{
  FILE *fil;

  fil=fopen(prefix,"r");
  if(fil!=NULL) {
    fclose(fil);
    return 1;
  }
  else {
    int nfils=0;
    while(nfils>=0) {
      char fname[256];
      sprintf(fname,"%s.%d",prefix,nfils);
      fil=fopen(fname,"r");
      if(fil!=NULL) {
	fclose(fil);
	nfils++;
      }
      else {
	if(nfils==0) {
	  fprintf(stderr,"CUTE: can't find file %s or %s.x\n",prefix,prefix);
	  return -1;
	}
	else if(nfils==1) {
	  fprintf(stderr,"CUTE: only file %s found. Weird.\n",fname);
	  return -1;
	}
	else {
	  return nfils;
	}
      }
    }
  }
  
  fprintf(stderr,"CUTE: this shouldn't have happened \n");
  return -1;
}

static Catalog read_snapshot_single(char *fname,lint *np,int input)
{
  //////
  // Creates catalog from a single snapshot file
  Catalog cat;
  lint ii;
  gad_header head;
  int block1,block2;

  FILE *snap=fopen(fname,"r");
  if(snap==NULL) error_open_file(fname);
  
  //Read header
  if(input==2)
    gad_seek_block(snap,"HEAD");
  my_fread(&block1,sizeof(int),1,snap);
  my_fread(&head,sizeof(gad_header),1,snap);
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);

  if(head.num_files!=1) {
    fprintf(stderr,"CUTE: Multi-file input not expected \n");
    exit(1);
  }
  
#ifdef _VERBOSE
  printf("  The cosmological model is:\n");
  printf("   - Omega_M = %.3lf\n",head.Omega0);
  printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
  printf("   - h = %.3lf\n",head.HubbleParam);
  printf("  This file contains: \n");
  for(ii=0;ii<6;ii++) {
    printf("   - %d particles of type %d with mass",
	   head.npart[ii],(int)ii);
    printf(" %.3lE (%d in total)\n",
	   head.mass[ii],head.npartTotal[ii]);
  }
  printf("  The box size is %.3lf\n",head.BoxSize);
  printf("  Redshift z = %.3lf \n",head.redshift);
#endif //_VERBOSE

  l_box=head.BoxSize;
  l_box_half=l_box*0.5;
  cat.np=0;
  for(ii=0;ii<6;ii++) {
    cat.np+=head.npart[ii];
    if(head.npart[ii]!=head.npartTotal[ii]) {
      fprintf(stderr,"CUTE: error reading snapshot \n");
      exit(1);
    }
  }
  *np=cat.np;

  cat.pos=(float *)malloc(3*cat.np*sizeof(float));
  if(cat.pos==NULL)
    error_mem_out();

  if(input==2)
    gad_seek_block(snap,"POS");
  my_fread(&block1,sizeof(int),1,snap);
  for(ii=0;ii<cat.np;ii++) {
    float pos[3];
    my_fread(pos,sizeof(float),3,snap);
    cat.pos[3*ii]=(float)(pos[0]);
    cat.pos[3*ii+1]=(float)(pos[1]);
    cat.pos[3*ii+2]=(float)(pos[2]);
  }
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);
  fclose(snap);

  return cat;
}



static Catalog read_hdf5(char *prefix, lint *np) {
	// change this to 1 if we want to compute redshift space correlation function from millgas
	int redshift_space = 0;
	int nfiles=512;
	int i;
	Catalog cat;
	char fname[256];
	printf("  Reading %d hdf5 files \n",nfiles);

	int rank;
	hid_t file_id;
	hid_t dataset_id;
	hid_t data_type;
	hid_t attr_id;
	herr_t status;
	hid_t filespace;
	hid_t memspace;
	const char *attrname = "Time";
	hsize_t dims[2];
	// This is the number of particles in JH's down-sampled hdf5 file for p-millennium
	//long long int numparts = 10000214;
	//double sampling_rate = 1.0;
	// these were used for measuring the millennium correlation function 
	//long long int numparts = 10077696000;
	long long int numparts = 100776960;
	float sampling_rate = 0.001;
	float *coords;
	float *vel;
	double scalefact;
	long long int partind=0;
	long long int j=0;
	float num;

	cat.np=(lint) (numparts);
	printf("reading %ld particles\n", cat.np);
	cat.pos=(float *)malloc(3*cat.np*sizeof(float));
	if(cat.pos==NULL)
	  error_mem_out();

	for (i=0; i<nfiles; i++) {
		sprintf(fname,"%s.%d.hdf5",prefix, i);
		//sprintf(fname,"%s.hdf5",prefix);
		file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (file_id < 0) {
		  printf("Cannot open %s \n",fname);
		  break;
		}
		else {
		  printf("Opened %s \n",fname);
		}

		dataset_id=H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
		filespace=H5Dget_space(dataset_id);

		rank=H5Sget_simple_extent_ndims(filespace);
		status=H5Sget_simple_extent_dims(filespace, dims, NULL);

		memspace=H5Screate_simple(rank, dims, NULL);
		coords=(float *) malloc(dims[0]*dims[1]*sizeof(float));
		if(coords==NULL)
		  error_mem_out();

		status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, coords);

		H5Sclose(memspace);
		data_type=H5Dget_type(dataset_id);
		status=H5Tclose(data_type);
		status=H5Dclose(dataset_id);

		if (redshift_space) {
			dataset_id=H5Dopen2(file_id, "/PartType1/Velocities", H5P_DEFAULT);
			data_type=H5Dget_type(dataset_id);
		
			filespace=H5Dget_space(dataset_id);
		
			rank=H5Sget_simple_extent_ndims(filespace);
			status=H5Sget_simple_extent_dims(filespace, dims, NULL);

			memspace=H5Screate_simple(rank, dims, NULL);
			vel=(float *) malloc(dims[0]*dims[1]*sizeof(float));
			status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, vel);

		
			H5Sclose(memspace);
			status=H5Tclose(data_type);
			status=H5Dclose(dataset_id);


			// get scale factor from snapshot file
			attr_id=H5Aopen_by_name(file_id, "/Header", attrname, H5P_DEFAULT, H5P_DEFAULT);

			status=H5Aread(attr_id, H5T_NATIVE_DOUBLE, &scalefact);
		}

		srand(time(0)); // Use current time as a seed
		for (j=0; j<dims[0]; j++) {
		  if (partind == cat.np) {
		    break;
		  }
		  num = (float)rand()/(float)RAND_MAX ;
		  // randomly down-sample particles
		  if (num < sampling_rate) {
		    cat.pos[3*partind]=coords[3*j] ;
		    cat.pos[3*partind+1]=coords[3*j+1];
		    cat.pos[3*partind+2]=coords[3*j+2];
		    
		    if (redshift_space) {
		      cat.pos[3*partind]=coords[3*j]+vel[3*j]/100.*sqrt(scalefact);
		      cat.pos[3*partind]=wrap_double(cat.pos[3*partind]);
		    }
		    
		    partind++;
		    if (partind == cat.np) {
		      break;
		    }
		  }

		}

		free(coords);
		if (redshift_space) {
			free(vel);
		}
		status=H5Fclose(file_id);


	}

	printf("finished reading files. Number of particles read: %lld, %ld\n", partind, cat.np);
	printf("%f %f %f \n", cat.pos[0], cat.pos[1], cat.pos[2]);

	cat.np=partind;
	*np=partind;
	return cat;

}



static Catalog read_gadget(char *prefix,lint *np,int input)
{
  Catalog cat;
  int nfils=check_num_files(prefix);
  if(nfils<=0) exit(1);

#ifdef _VERBOSE
  printf("  Reading from GADGET snapshot format \n");
#endif //_VERBOSE
  
  if(nfils==1) {
    printf("  Reading single snapshot file\n");
    cat=read_snapshot_single(prefix,np,input);
    return cat;
  }
  else {
    lint ii;
    char fname[256];
    gad_header head;
    int block1,block2;
    FILE *snap;

    printf("  Reading %d snapshot files \n",nfils);
    //    fprintf(stderr,"CUTE: multi-file input not supported \n");
    //    exit(1);

    sprintf(fname,"%s.0",prefix);
    snap=fopen(fname,"r");
    if(snap==NULL) error_open_file(fname);
  
    //Read header
    if(input==2)
      gad_seek_block(snap,"HEAD");
    my_fread(&block1,sizeof(int),1,snap);
    my_fread(&head,sizeof(gad_header),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);

    if(head.num_files!=nfils) {
      fprintf(stderr,
	      "CUTE: Header and existing files do not match %d != %d.\n",
	      nfils,head.num_files);
      fprintf(stderr,"      There may be some files missing\n");
      exit(1);
    }

#ifdef _VERBOSE
    printf("  The cosmological model is:\n");
    printf("   - Omega_M = %.3lf\n",head.Omega0);
    printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
    printf("   - h = %.3lf\n",head.HubbleParam);
    printf("  This file contains: \n");
    for(ii=0;ii<6;ii++) {
      printf("   - %d particles of type %d with mass %.3lE\n",
	     head.npartTotal[ii],(int)ii,head.mass[ii]);
    }
    printf("  The box size is %.3lf\n",head.BoxSize);
    printf("  Redshift z = %.3lf \n",head.redshift);
#endif //_VERBOSE

    l_box=head.BoxSize;
    l_box_half=l_box*0.5;
    cat.np=0;
    for(ii=0;ii<6;ii++)
      cat.np+=head.npartTotal[ii];
    *np=cat.np;

    cat.pos=(float *)malloc(3*cat.np*sizeof(float));
    if(cat.pos==NULL)
      error_mem_out();
    fclose(snap);

    lint np_read=0;
    for(ii=0;ii<nfils;ii++) {
      lint np_new;
      lint jj;

      sprintf(fname,"%s.%d",prefix,(int)ii);
      snap=fopen(fname,"r");
      if(snap==NULL) error_open_file(fname);
#ifdef _VERBOSE
      printf("  Reading file  %s \n",fname);
#endif //_VERBOSE

      //Read header
      if(input==2)
	gad_seek_block(snap,"HEAD");
      my_fread(&block1,sizeof(int),1,snap);
      my_fread(&head,sizeof(gad_header),1,snap);
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      np_new=0;
      for(jj=0;jj<6;jj++)
	np_new+=head.npart[jj];
      printf("  %ld parts in file %ld \n",(long)np_new,(long)ii);

      if(np_read+np_new>cat.np) {
	fprintf(stderr,
		"CUTE: files seem to contain too many particles\n");
	fprintf(stderr,"      file %s, %ld > %ld \n",
		fname,(long)(np_read+np_new),(long)(cat.np));
	exit(1);
      }

      if(input==2)
	gad_seek_block(snap,"POS");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=np_read;jj<np_read+np_new;jj++) {
	float pos[3];
	my_fread(pos,sizeof(float),3,snap);
	cat.pos[3*jj]=(float)(pos[0]);
	cat.pos[3*jj+1]=(float)(pos[1]);
	cat.pos[3*jj+2]=(float)(pos[2]);
      }
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);
      fclose(snap);

      np_read+=np_new;
    }
    if(np_read!=cat.np) {
      fprintf(stderr,
	      "CUTE: #particles read disagrees with header: %ld != %ld\n",
	      (long)np_read,(long)(cat.np));
      exit(1);
    }
    
    return cat;
  }
}
// Need to edit this
Catalog read_catalog(char *fname,lint *np)
{
  //////
  // Creates catalog from file fname
  lint ii;
  float x_mean=0,y_mean=0,z_mean=0;
  Catalog cat;

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");

  if((input_format==1) || (input_format==2))
    cat=read_gadget(fname,np,input_format);
  else if (input_format == 5)
    cat=read_hdf5(fname, np);
  else
    cat=read_ascii(fname,np);
  
  //Correct particles out of bounds and calculate CoM
  for(ii=0;ii<cat.np;ii++) {
    float xx,yy,zz;
    xx=cat.pos[3*ii];
    yy=cat.pos[3*ii+1];
    zz=cat.pos[3*ii+2];
    if((xx<0)||(xx>=l_box)) xx=wrap_double(xx);
    if((yy<0)||(yy>=l_box)) yy=wrap_double(yy);
    if((zz<0)||(zz>=l_box)) zz=wrap_double(yy);
    cat.pos[3*ii]=xx;
    cat.pos[3*ii+1]=yy;
    cat.pos[3*ii+2]=zz;
    x_mean+=xx/cat.np;
    y_mean+=yy/cat.np;
    z_mean+=zz/cat.np;
  }

#ifdef _VERBOSE
  printf("  The center of mass is (%.3lf,%.3lf,%.3lf) \n",
	 x_mean,y_mean,z_mean);
#endif //_VERBOSE
  
  printf("\n");
  return cat;
}
