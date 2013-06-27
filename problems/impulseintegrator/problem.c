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
	softening 	= 0.0;		
	dt		= input_get_double(argc,argv,"dt",0.0001);
	writeBest	= input_get_int(argc,argv,"write",0);

	boxsize 	= 10;
	tmax		= 1;
	root_nx = 1; root_ny = 1; root_nz = 1;
	nghostx = 0; nghosty = 0; nghostz = 0; 		
	init_box();
	
	int _N = 100;
	for (int i=0;i<_N;i++){
		struct particle star;
		double r = powf(powf(tools_uniform(0,1),-2./3.)-1.,-1./2.);
		double x2 = tools_uniform(0,1);
		double x3 = tools_uniform(0,2.*M_PI);
		star.z = (1.-2.*x2)*r;
		star.x = sqrt(r*r-star.z*star.z)*cos(x3);
		star.y = sqrt(r*r-star.z*star.z)*sin(x3);

		particles_add(star);

	}

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

double calculate_energy(struct particle* _particles, int _N){
	double energy = 0;
	for (int i=0;i<_N;i++){
		for (int j=0;j<i;j++){
			double dx = _particles[i].x - particles[j].x;
			double dy = _particles[i].y - particles[j].y;
			double dz = _particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy += -G*_particles[i].m*_particles[j].m/r;
		}
		double dvx = _particles[i].vx;
		double dvy = _particles[i].vy;
		double dvz = _particles[i].vz;
		energy += 1./2. * _particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	return energy;
}


void calculate_error(){
	FILE* inf = fopen("best.bin","rb"); 
	int _N;
	double _t;
	fread(&_N,sizeof(int),1,inf);
	fread(&_t,sizeof(double),1,inf);
	struct particle* _particles = malloc(_N*sizeof(struct particle));
	for (int i=0;i<_N;i++){
		fread(&(_particles[i]),sizeof(struct particle),1,inf);
	}
	fclose(inf);
	double dif_pos = 0;
	double dif_posabs = 0;
	double dif_vel = 0;
	for (int i=0;i<N;i++){
		double dx = _particles[i].x - particles[i].x;
		double dy = _particles[i].y - particles[i].y;
		double dz = _particles[i].z - particles[i].z;
		dif_pos += sqrt(dx*dx + dy*dy + dz*dz )/(double)N;
		dif_posabs += fabs(dx)/(double)N;
		dif_posabs += fabs(dy)/(double)N;
		dif_posabs += fabs(dz)/(double)N;
		double dvx = _particles[i].vx - particles[i].vx;
		double dvy = _particles[i].vy - particles[i].vy;
		double dvz = _particles[i].vz - particles[i].vz;
		dif_vel += sqrt(dvx*dvx + dvy*dvy + dvz*dvz )/(double)N;
	}
	double energy_best = calculate_energy(_particles,_N);
	double energy      = calculate_energy(particles,N);
	double dif_energy = fabs((energy-energy_best)/energy_best);
#ifdef INTEGRATOR_LEAPFROG
	FILE* of = fopen("error_leapfrog.txt","a+"); 
#endif
#ifdef INTEGRATOR_EULER
	FILE* of = fopen("error_euler.txt","a+"); 
#endif
#ifdef INTEGRATOR_MODIFIEDEULER
	FILE* of = fopen("error_modifiedeuler.txt","a+"); 
#endif
#ifdef INTEGRATOR_RADAU15
	FILE* of = fopen("error_radau15.txt","a+"); 
#endif
#ifdef INTEGRATOR_IMPULSE
	FILE* of = fopen("error_impulse.txt","a+"); 
#endif
#ifdef INTEGRATOR_IMPULSE2NDORDER
	FILE* of = fopen("error_impulse2ndorder.txt","a+"); 
#endif
	fprintf(of,"%e\t",dt);
	fprintf(of,"%e\t",dif_pos);
	fprintf(of,"%e\t",dif_vel);
	fprintf(of,"%e\t",dif_energy);
	fprintf(of,"%e\t",dif_posabs);
	fprintf(of,"\n");
	fclose(of);
	FILE* ofp = fopen("position.txt","a+"); 
	fprintf(ofp,"%e\t%e\t%e\t%e\n",dt,particles[0].x,particles[0].y,particles[0].z);
	fclose(ofp);
	
}

struct vec3o {
	double x;
	double y;
	double z;
};


void problem_finish(){
	printf("\nFinished. Time: %e. Difference to exact time: %e. Timestep: %e\n",t,t-tmax,dt);
	if (writeBest){
		output_binary("best.bin");
	}
	calculate_error();
}
