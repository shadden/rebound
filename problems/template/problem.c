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

int count = 0;

void problem_init(int argc, char* argv[]){
	if (argc>1){ 						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}else{
		boxsize = 10;
	}
	init_box();
	
	struct particle pt;
	double vv,mm,aa;

	dt = 0.1;
	mm = 2.;
	aa = 1.;
	vv = sqrt(G*(mm+mm)/aa)/4;
	printf("  G: %lf; m: %lf; a: %lf; v: %lf\n",G,mm,aa,vv);

	// particle 1
	pt.x 	=-aa/2;	pt.y 	= 0; 	pt.z 	= 0;
	pt.vx 	= 0; 	pt.vy 	= vv; 	pt.vz 	= 0;
	pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
	pt.m 	= mm;
	particles_add(pt);

	// particle 2
	pt.x 	= aa/2;	pt.y 	= 0; 	pt.z 	= 0;
	pt.vx 	= 0; 	pt.vy 	= -vv; 	pt.vz 	= 0;
	pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
	pt.m 	= mm;
	particles_add(pt);
}

void problem_inloop(){
}

void problem_output(){
  //        if(output_check(4000.*dt)){
  //              output_timing();
  //      }
  printf("%d; %lf\n",count++,dt);
  output_append_orbits("orbits.txt");
  //      if(output_check(12.)){
  //      }
}

void problem_finish(){
}
