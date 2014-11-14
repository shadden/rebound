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
#include "parse_file.h"
#include <assert.h>

#define CARTESIAN 1
#define ORBEL2D 2

int i,Npl,N_dt_out;
int inputType;
double dtdump;
double Mcen;
char *parfile,*posfile;

void problem_init(int argc, char* argv[]){

	assert(argc==3);	
	parfile = argv[1];
	posfile = argv[2];
	
	if (parse_param_data(fopen(parfile,"r"),&Npl,&inputType,&dt,&tmax,&N_dt_out,&Mcen,&boxsize) != 1)
		exit(-1);
	
	// Initialize box
	init_box();
	// Setup constants
	dtdump = N_dt_out * dt;

	if(inputType==CARTESIAN)
		parse_cartesian_data(fopen(posfile,"r"),Npl,Mcen);
	else if(inputType==ORBEL2D)
		parse_orbel2D_data(fopen(posfile,"r"),Npl,Mcen);
	else{
		printf("Unsupported input type!");
		exit(-1);
	}
	tools_move_to_center_of_momentum();
	system("rm -f orbits.txt");	
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(dtdump)){
		output_timing();
		output_append_orbits("orbits.txt");
	}
}

void problem_finish(){
}
