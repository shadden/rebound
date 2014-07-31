#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void parse_param_data(FILE *fi,int* npl,double* dt,int* N_dt_out,int* N_dt_timing){
	if(fscanf(fi,"%d",npl)!=1)
		printf("Improper value, line 1");
	if(fscanf(fi,"%lf",dt)!=1)
		printf("Improper value, line 2");
	if(fscanf(fi,"%d",N_dt_out)!=1)
		printf("Improper value, line 3");
	if(fscanf(fi,"%d",N_dt_timing)!=1)
		printf("Improper value, line 4");
	
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

