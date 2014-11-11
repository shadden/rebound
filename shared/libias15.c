/**
 * @file 	libias15.c
 * @brief 	Declarations and functions for shared library access to the IAS15 integrator.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Dave Spiegel <dave@ias.edu>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel
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
#include <string.h>
#include "particle.h"
#include "integrator.h"
#include "gravity.h"
#include "particle.h"
#include "main.h"
#include "boundaries.h"

// Default values of parameters and constants
double dt 	= 0.01;	
double t 	= 0;
double tmax	= 0;
double G 	= 1;
double softening = 0;

// Function pointer to additional forces
void (*problem_additional_forces) () = NULL;

// Particle structure
struct particle* particles = NULL;

// Particle getter/setter methods.
void setp(struct particle* _p){
	free(particles);
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle particle_get(int i){
	return particles[i];
}
struct particle* particles_get(){
	return particles;
}
void set_additional_forces(void (* _cb)()){
	problem_additional_forces = _cb;
}

// Integrate for 1 step
void step(){ 
	if (N<=0){
		fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
		return;
	}
	integrator_part1();
	gravity_calculate_acceleration();
	if (problem_additional_forces) problem_additional_forces();
	integrator_part2();
}

// Integrate until t=_tmax
void integrate(double _tmax){
	tmax = _tmax;
	double dt_last_done = dt;
	while(t<tmax){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return;
		}
		step();
		if (t+dt>=tmax){
			dt = tmax-t;
		}else{
			dt_last_done = dt;
		}
	}
	dt = dt_last_done;
}
	 
// The following code is needed to leave the original REBOUND files unchanged. 
// Currently the libias15 library only supports an infinitely large box.
// Infinite box size
int root_nx = 1;
int root_ny = 1;
int root_nz = 1;
double boxsize = 0;
double boxsize_x = 0;
double boxsize_y = 0;
double boxsize_z = 0;

int boundaries_particle_is_in_box(struct particle p){
	return 1;
}

// No ghost boxes for now.
int nghostx = 0;	
int nghosty = 0;
int nghostz = 0;

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = 0;
	gb.shifty = 0;
	gb.shiftz = 0;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}
