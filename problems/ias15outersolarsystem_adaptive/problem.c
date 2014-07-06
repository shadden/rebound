/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate the outer planets of the solar system. The initial 
 * conditions are taken from Applegate et al 1986. Pluto is a test
 * particle.
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
#include "main.h"
#include "output.h"
#include "input.h"
#include "integrator.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

double energy();
double energy_init;
void input_binary_special(char* filename);
int init_N;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",1.);			// days
	const double k 	= 0.01720209895;	// Gaussian constant 
	G		= k*k;

	integrator_epsilon = input_get_double(argc,argv,"integrator_epsilon",0.01);
	integrator_force_is_velocitydependent = 0;
	
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(150); 			// Init box with width 200 astronomical units

	input_binary_special("particles.bin");
	init_N = N;

#ifndef INTEGRATOR_WH
	// Move to barycentric frame
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH
	mpf_set_default_prec(512);
	energy_init = energy();
}

void problem_inloop(){
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
void problem_output(){
	/*
	if(output_check(3.)){
		struct orbit o = tools_p2orbit(particles[1], particles[0]);
		FILE* of = fopen("orbit.txt","a+"); 
		fprintf(of,"%e\t%e\t%e\n",t,o.a,o.e);
		fclose(of);
	}
	*/
}

void problem_finish(){
	FILE* of = fopen("energy.txt","w"); 
	double rel_energy = fabs((energy()-energy_init)/energy_init);
	if (init_N!=N){
		rel_energy = -1;
	}
	fprintf(of,"%e",rel_energy);
	fclose(of);
}

void input_binary_special(char* filename){
	FILE* inf = fopen(filename,"rb"); 
	long objects = 0;
	int _N;
	objects += fread(&_N,sizeof(int),1,inf);
	objects += fread(&tmax,sizeof(double),1,inf);
	for (int i=0;i<_N;i++){
		struct particle p;
		objects += fread(&p,sizeof(struct particle),1,inf);
		particles_add(p);
	}
	fclose(inf);
}
