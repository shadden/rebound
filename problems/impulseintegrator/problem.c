/**
 * @file 	problem.c
 * @brief 	Example problem: self-gravity disc.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A self-gravitating disc is integrated using
 * the leap frog integrator. This example is also compatible with 
 * the Wisdom Holman integrator. Collisions are not resolved.
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
#include "input.h"
#include "tools.h"

extern int Nmax;
int writeBest; 

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	softening 	= 0.01;		
	dt		= input_get_double(argc,argv,"dt",0.0001);
	writeBest		= input_get_int(argc,argv,"write",0);
	boxsize 	= 100;
	tmax		= 1;
	root_nx = 1; root_ny = 1; root_nz = 1;
	nghostx = 0; nghosty = 0; nghostz = 0; 		
	init_box();
	
	/**
	 * This function sets up a Plummer sphere.
	 * @details This function is based on a routine from the NEMO package, P. Teuben (1995).
	 * For details on the implementation see the Appendix of Aarseth, Henon and Wielen (1974). 
	 * @param _N Number of particles in the plummer sphere.
	 * @param mlow Lower mass fraction cutoff (can be 0).
	 * @param rfrac Upper radius cutoff (the Plummer sphere is formally an inifitely large object).
	 * @param quiet Noisyness of the model, 0=noise, 1=medium, 2=quiet.
	 * @param scale Scales the final model before adding it to the simulation.
	 * @param shift Shift the entire sphere in position and velocity space (6 values). 
	 */
	//double shift[6] = {0,0,0,0,0,0};
	srand(0);
	//tools_init_plummer(100, 0., 1, 0, 1, shift);
	int totalN = 100;
	do{
		struct particle star;
		star.m = 1;
		double r = tools_uniform(0,1) ;
		double _r;
		do{
			star.x = tools_uniform(-1,1);
			star.y = tools_uniform(-1,1);
			star.z = tools_uniform(-1,1);
			_r = sqrt(star.x*star.x + star.y*star.y  + star.z*star.z);
		}while(_r>1.);
		star.x *= r/_r;
		star.y *= r/_r;
		star.z *= r/_r;
	
		double _v;	
		do{
			star.vx = tools_uniform(-1,1);
			star.vy = tools_uniform(-1,1);
			star.vz = tools_uniform(-1,1);
			_v = sqrt(star.vx*star.vx + star.vy*star.vy  + star.vz*star.vz);
		}while(_v>1.);
		double v = 0.1*sqrt(2.*(double)totalN/r);
		star.vx *= v/_v;
		star.vy *= v/_v;
		star.vz *= v/_v;
		
		particles_add(star);
	}while(N<totalN);

	double v = 0;
	for (int i=0;i<N;i++){
		v += sqrt(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz);
	}
	printf("System size: %f\n", 1.);
	printf("Smoothing length: %f\n", softening);
	printf("Interparticle separation: %f\n", 1./powf((double)N,1./3.));
	printf("Characteristic velocity: %f\n", sqrt(G*(double)N*particles[0].m/1.));
	printf("Mean velocity: %f\n", v/(double)N);
	printf("Crossing time: %f\n", 1./(v/(double)N));
	printf("Encounter time: %f\n", 1./(v/(double)N)/powf((double)N,1./3.));
	printf("Timestep: %f\n", dt);
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.0*dt)) output_timing();
}


void calculate_error(){
	FILE* inf = fopen("best.bin","rb"); 
	double* pos = malloc(3*N*sizeof(double));
	double* vel = malloc(3*N*sizeof(double));
	for (int i=0;i<N;i++){
		fread(&(pos[3*i]),sizeof(double)*3,1,inf);
		fread(&(vel[3*i]),sizeof(double)*3,1,inf);
	}
	fclose(inf);
	double dif_pos = 0;
	double dif_posabs = 0;
	double dif_vel = 0;
	for (int i=0;i<N;i++){
		double dx = pos[3*i+0] - particles[i].x;
		double dy = pos[3*i+1] - particles[i].y;
		double dz = pos[3*i+2] - particles[i].z;
		dif_pos += sqrt(dx*dx + dy*dy + dz*dz )/(double)N;
		dif_posabs += fabs(dx)/(double)N;
		dif_posabs += fabs(dy)/(double)N;
		dif_posabs += fabs(dz)/(double)N;
		double dvx = vel[3*i+0] - particles[i].vx;
		double dvy = vel[3*i+1] - particles[i].vy;
		double dvz = vel[3*i+2] - particles[i].vz;
		dif_vel += sqrt(dvx*dvx + dvy*dvy + dvz*dvz )/(double)N;
	}
	double energy_best = 0;
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			if (i<j){
				double dx = pos[3*i+0] - pos[3*j+0];
				double dy = pos[3*i+1] - pos[3*j+1];
				double dz = pos[3*i+2] - pos[3*j+2];
				double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
				energy_best -= 1./r;
			}
		}
		double dvx = vel[3*i+0];
		double dvy = vel[3*i+1];
		double dvz = vel[3*i+2];
		energy_best += 1./2. * (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	double energy = 0;
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			if (i<j){
				double dx = particles[i].x - particles[j].x;
				double dy = particles[i].y - particles[j].y;
				double dz = particles[i].z - particles[j].z;
				double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
				energy -= 1./r;
			}
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		energy += 1./2. * (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	double dif_energy = fabs(energy-energy_best)/(double)N;
	FILE* of = fopen("error.txt","a+"); 
	fprintf(of,"%e\t",dt);
	fprintf(of,"%e\t",dif_pos);
	fprintf(of,"%e\t",dif_vel);
	fprintf(of,"%e\t",dif_energy);
	fprintf(of,"%e\t",dif_posabs);
	fprintf(of,"\n");
	fclose(of);
}

struct vec3o {
	double x;
	double y;
	double z;
};

void output_binaryxv(char* filename){
	FILE* of = fopen(filename,"wb"); 
	for (int i=0;i<N;i++){
		struct vec3o v;
		v.x = particles[i].x;
		v.y = particles[i].y;
		v.z = particles[i].z;
		fwrite(&(v),sizeof(struct vec3o),1,of);
		v.x = particles[i].vx;
		v.y = particles[i].vy;
		v.z = particles[i].vz;
		fwrite(&(v),sizeof(struct vec3o),1,of);
	}
	fclose(of);
}

void problem_finish(){
	printf("\nFinished. Time: %f\n",t);
	if (writeBest){
		output_binaryxv("best.bin");
	}
	calculate_error();
}
