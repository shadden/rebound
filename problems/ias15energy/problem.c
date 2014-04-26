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
#include <gmp.h>
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

double energy();
double energy_init = 0; 
double energy_max = 0; 
double ecc =0;
extern double integrator_epsilon;
extern double integrator_error;
double timing_start;
double timing_stop;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		

	// Setup homog. sphere
	tmax		= 1e4*sqrt(2.)*M_PI;
	dt		= tmax/input_get_double(argc,argv,"timesteps",1000);
#ifdef INTEGRATOR_IAS15
	//integrator_epsilon = input_get_double(argc,argv,"epsilon",0);
#endif // INTEGRATOR_IAS15
	boxsize = 500;	
	init_box();

	ecc=input_get_double(argc,argv,"e",0);

	// star
	struct particle star;
	star.x = 0.0; 
	star.y = 0.0; 
	star.z = 0; 
	star.vx = 0; 
	star.vy = 0.; 
	star.vz = 0; 
	star.m = 1;
	particles_add(star);
	
	{
		struct particle p2;
		p2.x = 1.+ecc; 
		p2.y = 0; 
		p2.z = 0; 
		p2.vx = 0; 
		p2.vy = sqrt((1.-ecc)/(1.+ecc)*1./p2.x); 
		p2.vz = 0; 
		p2.m = 1e-3;
		particles_add(p2);
 	}
	{
		struct particle p2;
		p2.x = 2.+ecc; 
		p2.y = 0; 
		p2.z = 0; 
		p2.vx = 0; 
		p2.vy = sqrt((1.-ecc)/(1.+ecc)*1./p2.x); 
		p2.vz = 0; 
		p2.m = 1e-3;
		particles_add(p2);
 	}

#ifndef INTEGRATOR_WH
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH
	mpf_set_default_prec(512);
	energy_init = energy();
	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_start = tim.tv_sec+(tim.tv_usec/1000000.0);
}

double energy(){
	mpf_t energy_kinetic;
	mpf_init(energy_kinetic);
	mpf_t energy_potential;
	mpf_init(energy_potential);
	mpf_t mass;
	mpf_init(mass);
	mpf_t vx;
	mpf_init(vx);
	mpf_t vy;
	mpf_init(vy);
	mpf_t vz;
	mpf_init(vz);
	for (int i=0;i<N;i++){
		mpf_t temp;
		mpf_init(temp);
		mpf_set_d(temp,particles[i].vx*particles[i].m);
		mpf_add(vx, vx, temp);
		mpf_set_d(temp,particles[i].vy*particles[i].m);
		mpf_add(vy, vy, temp);
		mpf_set_d(temp,particles[i].vz*particles[i].m);
		mpf_add(vz, vz, temp);

		mpf_set_d(temp,particles[i].m);
		mpf_add(mass, mass, temp);
	}
	mpf_div(vx,vx,mass);
	mpf_div(vy,vy,mass);
	mpf_div(vz,vz,mass);

	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			mpf_t temp;
			mpf_init(temp);
			mpf_set_d(temp,-G*particles[i].m*particles[j].m/r);
			mpf_add(energy_potential, energy_potential,temp);
			
		}
	
		double dvx = particles[i].vx-mpf_get_d(vx);
		double dvy = particles[i].vy-mpf_get_d(vy);
		double dvz = particles[i].vz-mpf_get_d(vz);
		mpf_t temp;
		mpf_init(temp);
		
		mpf_set_d(temp,1./2.*particles[i].m * (dvx*dvx + dvy*dvy + dvz*dvz));
		mpf_add(energy_kinetic, energy_kinetic,temp);
		
	}
	mpf_add(energy_kinetic,energy_kinetic,energy_potential);

	return mpf_get_d(energy_kinetic);
}
void problem_inloop(){
}

//FILE* ofe = NULL; 
void problem_finish(){
	printf("\nFinished. Time: %e. Difference to exact time: %e. Timestep: %e\n",t,t-tmax,dt);
#ifdef INTEGRATOR_LEAPFROG
	FILE* of = fopen("energy_leapfrog.txt","a+"); 
#endif
#ifdef INTEGRATOR_IAS15
	FILE* of = fopen("energy_ias15.txt","a+"); 
#endif
#ifdef INTEGRATOR_WH
	FILE* of = fopen("energy_wh.txt","a+"); 
#endif
	fprintf(of,"%e\t",dt);
	fprintf(of,"%e\t",integrator_epsilon);
	fprintf(of,"%e\t",ecc);
	fprintf(of,"%e\t",energy_max);
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_stop = tim.tv_sec+(tim.tv_usec/1000000.0);
	fprintf(of,"%e\t",timing_stop-timing_start);
	fprintf(of,"\n");
	fclose(of);
}

void problem_output(){
	if(output_check(tmax/10000)){
		double en = fabs((energy()-energy_init)/energy_init);
		if (en>energy_max) energy_max = en;
	}
}

