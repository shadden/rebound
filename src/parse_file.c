#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tools.h"
#include "particle.h"

int parse_param_data(FILE *fi,int* Npl,int* inputType,double* dt,double* tmax,int* N_dt_out,double* Mcen,double* boxsize)
{
	if(fscanf(fi,"%d",Npl)!=1){
		printf("Improper value, line 1");
		return 0;
	}

	if(fscanf(fi,"%d",inputType)!=1){
		printf("Improper value, line 2");
		return 0;
	}

	if(fscanf(fi,"%lf",dt)!=1){
		printf("Improper value, line 3");
		return 0;
	}

	if(fscanf(fi,"%lf",tmax)!=1){
		printf("Improper value, line 4");
		return 0;
	}
	
	if(fscanf(fi,"%d",N_dt_out)!=1){
		printf("Improper value, line 5");
		return 0;
	}
		
	if(fscanf(fi,"%lf",Mcen)!=1){
		printf("Improper value, line 6");	
		return 0;
	}
	
	if(fscanf(fi,"%lf",boxsize)!=1){
		printf("Improper value, line 7");	
		return 0;
	}
	
	return 1;
}

void parse_cartesian_data(FILE *fi,int Npl, double masses[Npl], double ** pos, double **vel){
	double mass;
	double x,y,z;
	double vx,vy,vz;
	int i =0;
	for(i=0; i<Npl;i++){
		if ( fscanf(fi,"%lf %lf %lf %lf %lf %lf %lf",&mass,&x,&y,&z,&vx,&vy,&vz)!=7){
			printf("Error in planet file, line %d\n",i+1);
			exit(1);
		}
		masses[i] = mass;
		pos[i][0] = x;	pos[i][1] = y;	pos[i][2] = z;
		vel[i][0] = vx;	vel[i][1] = vy;	vel[i][2] = vz;
	}
}

void parse_orbel2D_data(FILE *fi,int Npl,double Mcen){
	double mass,a,e,omega,M;
	
	int i =0;
	for(i=0; i<Npl;i++){
		if ( fscanf(fi,"%lf %lf %lf %lf %lf",&mass,&a,&e,&omega,&M)!=5){
			printf("Error in planet file, line %d\n",i+1);
			exit(1);
		}
		
		double f = tools_MeanAnom2TrueAnom(M,e);
		particles_add( tools_init_orbit2d(Mcen,mass,a,e,omega,f) );
	}
}

