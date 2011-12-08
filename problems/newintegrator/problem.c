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

double energy();
double energy_init;

void problem_init(int argc, char* argv[]){
	boxsize = 100;
	softening = 0;//boxsize/100.;
	dt = 0.01;
	init_box();
	double velocity = 3;
	tmax = boxsize_x/4./velocity;


	while(N<200){
		struct particle pt;
		pt.x 	= tools_uniform(-boxsize_x/16.,boxsize_x/16.);
		pt.y 	= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
		if (fabs(pt.y)<boxsize_y/16.) continue;
		pt.z 	= 0;
		pt.vx 	= 0; 	pt.vy 	= 0; 	pt.vz 	= 0;
		pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
		pt.m 	= 1;
		particles_add(pt);
	}
	N_active = N;
	
	while(N<300){
		struct particle pt;
		pt.x 	= -boxsize_x/8.;
		pt.y 	= tools_uniform(-boxsize_y/32.,boxsize_y/32.);
		pt.z 	= 0;
		pt.vx 	= velocity; 	pt.vy 	= 0; 	pt.vz 	= 0;
		pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
		pt.m 	= 0;
		particles_add(pt);
	}

	energy_init = energy();
}

double energy(){
	double e = 0;
	for (int i=N_active; i<N; i++){
		for (int j=0; j<N_active; j++){
			if (i==j) continue;
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			e += -G/r*particles[j].m;
		}
	}
	return e;
}




void problem_inloop(){
}

void problem_output(){
	output_timing();
}

void problem_finish(){
	double energy_finish = energy();
	printf("\nrelative energy error = %e\n", (energy_finish-energy_init)/(energy_finish+energy_init)*2.);
}
