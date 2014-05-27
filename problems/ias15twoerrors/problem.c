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
#include "input.h"
#include "tools.h"

extern double integrator_epsilon;
extern double integrator_error;
double timing_start;
double ecc;
double timing_stop;

double velocity(){
	struct particle p = particles[1];
	double vx = p.vx;
	double vy = p.vy;
	double vz = p.vz;
	double v2 = vx*vx + vy*vy + vz*vz;
	return sqrt(v2);
}

double energy(){
	double energy = 0;
	struct particle p = particles[1];
	double vx = p.vx;
	double vy = p.vy;
	double vz = p.vz;
	double v2 = vx*vx + vy*vy + vz*vz;
	energy += 0.5 *v2;
	struct particle p2 = particles[0];
	double rx = p.x-p2.x;
	double ry = p.y-p2.y;
	double rz = p.z-p2.z;
	double r2 = rx*rx + ry*ry + rz*rz;
	energy -= G*p.m/sqrt(r2);
	return energy;
}
double energy_init;
double velocity_init;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		

	// Setup homog. sphere
	tmax		= 1e2*2.*M_PI;
	dt			= tmax/input_get_double(argc,argv,"timesteps",1000);
	integrator_epsilon	= input_get_double(argc,argv,"epsilon",0);
	boxsize = 500;	
	init_box();

	ecc=input_get_double(argc,argv,"e",0);

	// star
	struct particle star;
	star.x = 0.0; 
	star.y = 0.0; 
	star.z = 0; 
	star.vx = 0; 
	star.vy = 0.; 
	star.vz = 0; 
	star.m = 1;
	particles_add(star);
	
	{
		struct particle p2;
		p2.x = 1.+ecc; 
		p2.y = 0; 
		p2.z = 0; 
		p2.vx = 0; 
		p2.vy = sqrt((1.-ecc)/(1.+ecc)*1./p2.x); 
		p2.vz = 0; 
		p2.m = 0;
		particles_add(p2);
 	}

#ifndef INTEGRATOR_WH
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH

	energy_init = energy();	
	velocity_init = velocity();	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_start = tim.tv_sec+(tim.tv_usec/1000000.0);

}


void problem_finish(){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_stop = tim.tv_sec+(tim.tv_usec/1000000.0);
	
#ifdef INTEGRATOR_LEAPFROG
	FILE* of = fopen("energy_leapfrog.txt","a+"); 
#endif
#ifdef INTEGRATOR_IAS15
	FILE* of;
	if (integrator_epsilon>0){
		of = fopen("energy_ias15variable.txt","a+"); 
	}else{
		of = fopen("energy_ias15.txt","a+"); 
	}
#endif
#ifdef INTEGRATOR_WH
	FILE* of = fopen("energy_wh.txt","a+"); 
#endif
	fprintf(of,"%e\t",dt/2./M_PI);					// 1 steps per orbit
	fprintf(of,"%.10e\t",t);					// 2 t (tmax)
	fprintf(of,"%e\t",integrator_epsilon);				// 3 epsilon
	fprintf(of,"%e\t",timing_stop - timing_start);			// 4 timing
	struct particle p = particles[1];
	double r   = sqrt(p.x*p.x + p.y*p.y);
	double x = cos(t);
	double y = sin(t);
	double phi2 = atan2(y,x); 
	double phi = atan2(p.y,p.x); 
	fprintf(of,"%.10e\t",r-1.-ecc);					// 5 error off-track
	fprintf(of,"%.10e\t",(phi-phi2)/(2.*M_PI));			// 6 error on-track
	fprintf(of,"%e\t",tmax/dt);					// 7 Number of timesteps
	double energy_final = energy();	
	double velocity_final = velocity();	
	fprintf(of,"%e\t",(energy_final-energy_init)/energy_init);	// 8 Relagtive energy error (not for wh)
	fprintf(of,"%e\t",(velocity_final-velocity_init)/velocity_init);	// 9 relative velocity error (not for wh)
	fprintf(of,"\n");
	fclose(of);
	printf("\n\nt=%.20f\n\n",t/2./M_PI);
}

void problem_inloop(){
}

void problem_output(){
	if(output_check(tmax/10000)){
	}
//	problem_finish();

}

