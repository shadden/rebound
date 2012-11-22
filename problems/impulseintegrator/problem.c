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
#include "tools.h"

extern int Nmax;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	softening 	= 0.1;		
	dt 		= 3e-1;
	boxsize 	= 10;
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
	double shift[6] = {0,0,0,0,0,0};
	tools_init_plummer(100, 0., 1, 0, 1, shift);
	/*
	struct particle star = {
		.m=1,
		.x=0,.y=0,.z=0.,
		.vx=0,.vy=0,.vz=0.
		};
	particles_add(star);
	struct particle planet1 = {
		.m=0.001,
		.x=1,.y=0,.z=0.,
		.vx=0,.vy=1,.vz=0.
		};
	particles_add(planet1);
	struct particle planet2 = {
		.m=0.001,
		.x=0,.y=-1,.z=0.,
		.vx=1,.vy=0,.vz=0.
		};
	particles_add(planet2);
	struct particle planet3 = {
		.m=0.001,
		.x=0,.y=1,.z=0.,
		.vx=-1,.vy=0,.vz=0.
		};
	particles_add(planet3);
	struct particle planet4 = {
		.m=0.001,
		.x=-1,.y=0,.z=0.,
		.vx=0,.vy=-1,.vz=0.
		};
	particles_add(planet4);
	*/
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.0*dt)) output_timing();
}

void problem_finish(){
}
