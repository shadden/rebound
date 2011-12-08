/**
 * @file 	integrator.c
 * @brief 	Euler integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the Euler integration scheme.  
 * There is practically no use for this first order integrator. For 
 * non-rotating coordinates systems integrator_leapfrog should be used. 
 * For shearing sheet boundary conditions integrator_sei is well suited.
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

const int integrator_substep_N   = 2;
const enum integrator_substep_type integrator_substeps[5] = {IST_DRIFT, IST_KICK};

void integrator_part0();
void integrator_part1();

void integrator_part(int part){
	switch (part){
		case 0:
			integrator_part0();
			break;
		case 1:
			integrator_part1();
			break;
	}
}

// Drift
void integrator_part1(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += dt * particles[i].vx;
		particles[i].y  += dt * particles[i].vy;
		particles[i].z  += dt * particles[i].vz;
	}
	t+=dt;
}

// Kick
void integrator_part2(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
	}
}
	

