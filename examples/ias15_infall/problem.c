/**
 * @file 	problem.c
 * @brief 	Example problem: Infall.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "main.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	G 		= 1;
	dt 		= 1e-6;				// Initial timestep 	
	integrator_epsilon = 1e-3;			// Accuracy
	softening 	= 0.1;

	init_boxwidth(10); 

	int _N = 17;					// Number of particles
	
	double v = 0;
	if (0){ 					// Circular orbit instead of infall
		for (int i=1;i<_N;i++){
			v += 1./sin(M_PI*(double)i/(double)_N);
		}
		v = 1./2.*sqrt(v); 			// times sqrt(G*m/r) (=1)
	}

	for (int i=0;i<_N;i++){
		struct particle p;
		p.m  = 1; 
		p.x  =  cos((double)i/(double)_N*2.*M_PI);
		p.y  =  sin((double)i/(double)_N*2.*M_PI);
		p.z  = 0;
		p.vx = -v*sin((double)i/(double)_N*2.*M_PI);
		p.vy =  v*cos((double)i/(double)_N*2.*M_PI);
		p.vz = 0;
		particles_add(p); 
	}
	struct particle p;
	p.m = 1;
	p.x = 0; p.y = 0; p.z = 0;
	p.vx = 0; p.vy = 0; p.vz = 0;
	// This will break it:
	particles_add(p); 
}

void problem_inloop(){
}

void problem_output(){
	printf("%e\n",dt);
	output_timing();
}

void problem_finish(){
}
