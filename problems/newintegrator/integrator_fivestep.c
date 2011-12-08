/**
 * @file 	integrator.c
 * @brief 	McLachlan ABABA integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements a leap-frog style method of order
 * two. The first error term is O(epsilon^2 t^3), where epsilon is the
 * small parameter that measures how important gravity is.
 * Note that in a standard leapfrog algorithm, the error term is
 * O(epsilon t^3). 
 * For more details see McLauchlan (2003). 
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
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "integrator.h"


const int integrator_substep_N   = 5;
const enum integrator_substep_type integrator_substeps[5] = {IST_DRIFT, IST_KICK, IST_DRIFT, IST_KICK, IST_DRIFT};

void integrator_part0();
void integrator_part1();
void integrator_part2();

void integrator_part(int part){
	switch (part){
		case 0:
		case 4:
			integrator_part0();
			break;
		case 1:
		case 3:
			integrator_part1();
			break;
		case 2:
			integrator_part2();
			break;
	}
}

void integrator_part0(){
	const double c = (3.-sqrt(3.))/6.;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += c* dt * particles[i].vx;
		particles[i].y  += c* dt * particles[i].vy;
		particles[i].z  += c* dt * particles[i].vz;
	}
}

void integrator_part1(){
	const double c = 0.5;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += c* dt * particles[i].ax;
		particles[i].vy += c* dt * particles[i].ay;
		particles[i].vz += c* dt * particles[i].az;
	}
	t+=dt/2.;
}
	
void integrator_part2(){
	const double c = 1./sqrt(3.);
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += c* dt * particles[i].vx;
		particles[i].y  += c* dt * particles[i].vy;
		particles[i].z  += c* dt * particles[i].vz;
	}
}

