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
double energy_error_sum = 0;
int    energy_error_N = 0;

void problem_init(int argc, char* argv[]){
	boxsize   = 1;
	softening = 0.;
	init_box();
	double velocity = 30;
	tmax = 0.03;

	if (argc>1){ 	
		dt	= atof(argv[1]);
	}else{
		dt      = 1e-6;
	}

	while(N<100){
		struct particle pt;
		pt.x 	= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
		pt.y 	= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
		pt.z 	= tools_uniform(-boxsize_z/2.,boxsize_z/2.);
		if (fabs(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z)>boxsize_y*boxsize_y/128.) continue;
		pt.vx 	= 0; 	pt.vy 	= 0; 	pt.vz 	= 0;
		pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
		pt.m 	= 1;
		particles_add(pt);
	}
	
	energy_init = energy();
}

double energy(){
	double e = 0;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			if (i==j) continue;
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			e += -1./(r)*particles[j].m;
		}
		double vx = particles[i].vx;
		double vy = particles[i].vy;
		double vz = particles[i].vz;
		e += 0.5*(vx*vx + vy*vy + vz*vz);
	}
	return e/(double)(N-N_active);
}




void problem_inloop(){
}

void problem_output(){
	output_timing();
	double e = energy();
	energy_error_sum += fabs((e-energy_init)/energy_init);
	energy_error_N++;
}

void problem_finish(){
#ifdef INTEGRATOR_FIVESTEP
	FILE* f = fopen("error_fivestep.txt","a+");
#else  // FIVESTEPINTEGRATOR
	FILE* f = fopen("error_leapfrog.txt","a+");
#endif // FIVESTEPINTEGRATOR
	fprintf(f,"%e\t%e\n", dt,energy_error_sum/(double)energy_error_N);
	fclose(f);
}
