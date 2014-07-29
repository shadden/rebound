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

int Npl,N_dt_out,N_dt_timing;
char* parfile,posfile;

void problem_init(int argc, char* argv[]){
	if (argc>1){ 						// Try to read boxsize from command line
		parfile = argv[1];
		posfile = argv[2];
	}else{
		boxsize = 100;
	}
	init_box();

	parse_param_data(fopen(parfile,"r"),Npl,N_dt_out,N_dt_timing);
	
	double mass[Npl],pos[Npl][3],vel[Npl][3];
	parse_cartesian_data(fopen(posefile,"r"),Npl,mass,pos,vel);

	for (int i=0; i < Npl; i++){
		struct particle p;
		p.x  = pos[i][0]; 		p.y  = pos[i][1];	 	p.z  = pos[i][2];
		p.vx = vel[i][0]; 		p.vy = vel[i][1];	 	p.vz = vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = mass[i];
		particles_add(p); 
	}
	
	system("rm -f orbits.txt");
	
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(N_dt_timing*dt){
		output_timing();
	}
	if (output_check(N_dt_out*dt){
		output_append_orbits("orbits.txt");
	}
}

void problem_finish(){
}
