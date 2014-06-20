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

double kinetic_energy(){ // specific kinetic energy
	double v = velocity();
	return 0.5 * v*v;
}

double potential_energy(){ // specific potential energy (i.e., potential)
	struct particle s = particles[0];
	struct particle p = particles[1];
	double rx = p.x-s.x;
	double ry = p.y-s.y;
	double rz = p.z-s.z;
	double r2 = rx*rx + ry*ry + rz*rz;
	return G*s.m/sqrt(r2);
}

double energy(){
	double K = kinetic_energy();
	double U = kinetic_energy();
	return K + U;
}

double energy_init;
double K_init;
double U_init;
double velocity_init;
double orbitalperiod;
void problem_init(int argc, char* argv[]){
	G 		= 1.;		

	integrator_epsilon	= input_get_double(argc,argv,"epsilon",3.31e-3);
	double scale_a		= input_get_double(argc,argv,"scale_a",1.);
	double scale_m		= input_get_double(argc,argv,"scale_m",1.);
	double tmax_factor	= input_get_double(argc,argv,"tmax_factor",100.);
	orbitalperiod 		= 2.*M_PI*sqrt(scale_a*scale_a*scale_a/(G*scale_m));
	dt			= orbitalperiod/input_get_double(argc,argv,"timesteps",10);
	boxsize 		= 50.*scale_a;	
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
	star.m = 1.*scale_m;
	particles_add(star);
	
	{
		struct particle p2;
		double xx = (1.+ecc)*scale_a;
		double yy = 0.;
		double zz = 0.;
		double phi = tools_uniform(0.,2*M_PI); // pick random orientation
		p2.x = xx*cos(phi) - yy*sin(phi);
		p2.y = yy*cos(phi) + xx*sin(phi);
		p2.z = zz;
		double vvx = 0.;
		double vvy = sqrt(G*scale_m*(1.-ecc)/(1.+ecc)/scale_a);
		double vvz = 0.;
		p2.vx = vvx*cos(phi) - vvy*sin(phi);
		p2.vy = vvy*cos(phi) + vvx*sin(phi);
		p2.vz = vvz;
		p2.m = 0;
		particles_add(p2);
 	}

#ifndef INTEGRATOR_WH
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH

	energy_init = energy();	
	K_init = kinetic_energy();	
	U_init = potential_energy();	
	velocity_init = velocity();	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_start = tim.tv_sec+(tim.tv_usec/1000000.0);

	tmax = tmax_factor*orbitalperiod;
}


void problem_finish(){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_stop = tim.tv_sec+(tim.tv_usec/1000000.0);
	
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
#ifdef INTEGRATOR_IAS15ORIGINAL
	FILE* of = fopen("energy_ias15original.txt","a+"); 
#endif
	fprintf(of,"%e ",dt/2./M_PI);					// 1 steps per orbit
	fprintf(of,"%.10e ",t);					        // 2 t (tmax)
	fprintf(of,"%e ",integrator_epsilon);				// 3 epsilon
	fprintf(of,"%e ",timing_stop - timing_start);			// 4 timing
	struct particle p = particles[1];
	double r   = sqrt(p.x*p.x + p.y*p.y);
	double phi2 = t; 
	double phi = atan2(p.y,p.x); 
	fprintf(of,"%.10e ",r-1.-ecc);					// 5 error off-track
	double dphi = phi-phi2+M_PI;
	while (dphi<0.){
		dphi+= M_PI*2.;
	}
	while (dphi>2.*M_PI){
		dphi-= M_PI*2.;
	}
	dphi -=M_PI;
	fprintf(of,"%.10e ",dphi);			                // 6 error on-track
	fprintf(of,"%e ",tmax/dt);					// 7 Number of timesteps
	double energy_final = energy();	
	double K_final = kinetic_energy();	
	double U_final = potential_energy();	
	double velocity_final = velocity();	
	fprintf(of,"%e ",(energy_final-energy_init)/energy_init);	// 8 relative energy error (not for wh)
	fprintf(of,"%e ",(K_final-K_init)/K_init);	                // 9 relative kinetic energy error (not for wh)
	fprintf(of,"%e ",(U_final-U_init)/U_init);	                //10 relative potential energy error (not for wh)
	fprintf(of,"%e ",(velocity_final-velocity_init)/velocity_init); //11 relative velocity error (not for wh)
	fprintf(of,"\n");
	fclose(of);
}

void problem_inloop(){
}

void problem_output(){
	if(output_check(tmax/10000)){
	}
//	problem_finish();
//	printf("%e\t%e\n",dt/orbitalperiod,dt);
}

