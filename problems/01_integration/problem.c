/**
 * @file 	problem.c
 * @brief 	Template file for new problems.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include <assert.h>

int i,Npl,N_dt_out,N_dt;
double dtdump;
char *parfile,*posfile;

// Parsing Functions
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
//

void problem_init(int argc, char* argv[]){

	assert(argc==3);	
	parfile = argv[1];
	posfile = argv[2];
	boxsize=100;
	init_box();

	parse_param_data(fopen(parfile,"r"),&Npl,&dt,&N_dt,&N_dt_out);
	// Setup constants
	tmax = N_dt * dt;
	dtdump = N_dt_out *dt;

	double **pos,**vel,*masses;
	pos = (double **)malloc(sizeof(double *)*Npl);
	vel = (double **)malloc(sizeof(double *)*Npl);
	for(i=0; i < Npl; i++) {
	  pos[i] = (double *)malloc(sizeof(double)*3);
	  vel[i] = (double *)malloc(sizeof(double)*3);
	}
	masses = (double *)malloc(sizeof(double)*Npl);

	parse_cartesian_data(fopen(posfile,"r"),Npl,masses,pos,vel);

	for (int i=0; i < Npl; i++){
		struct particle p;
		p.x  = pos[i][0]; 		p.y  = pos[i][1];	 	p.z  = pos[i][2];
		p.vx = vel[i][0]; 		p.vy = vel[i][1];	 	p.vz = vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = masses[i];
		particles_add(p); 
	}
	
	system("rm -f orbits.txt");
	for(i=0; i < Npl; i++) {
	  free(pos[i]);
	  free(vel[i]);
	}
	free(pos);
	free(vel);
	free(masses);
	
}

void problem_inloop(){

}

void problem_output(){

	if (output_check(dtdump)){
		output_timing();
		output_append_orbits("orbits.txt");
	}
}

void problem_finish(){
}
