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

extern int Nmax;
int writeBest; 
double calculate_energy(struct particle* _particles, int _N);
double t0, r0;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	dt		= input_get_double(argc,argv,"dt",0.0001);
	writeBest	= input_get_int(argc,argv,"write",0);

	
	int _N = 100;
	double M = 1;
	double R = 1;
	double E = 3./64.*M_PI*M*M/R;
	r0 = 16./(3.*M_PI)*R;
	t0 = G*pow(M,5./2.)*pow(4.*E,-3./2.)*(double)_N/log(0.4*(double)_N); //Rellaxation time
	tmax		= 2.*t0;
	boxsize 	= 100.*r0;
	softening 	= 0.0001*r0;		

	init_box();
	
	
	
	// http://adsabs.harvard.edu//abs/1974A%26A....37..183A
	for (int i=0;i<_N;i++){
		struct particle star;
		double r = pow(pow(tools_uniform(0,1),-2./3.)-1.,-1./2.);
		double x2 = tools_uniform(0,1);
		double x3 = tools_uniform(0,2.*M_PI);
		star.z = (1.-2.*x2)*r;
		star.x = sqrt(r*r-star.z*star.z)*cos(x3);
		star.y = sqrt(r*r-star.z*star.z)*sin(x3);
		double x5,g,q;
		do{
			x5 = tools_uniform(0.,1.);
			q = tools_uniform(0.,1.);
			g = q*q*pow(1.-q*q,7./2.);
		}while(0.1*x5>g);
		double ve = pow(2.,1./2.)*pow(1.+r*r,-1./4.);
		double v = q*ve;
		double x6 = tools_uniform(0.,1.);
		double x7 = tools_uniform(0.,2.*M_PI);
		star.vz = (1.-2.*x6)*v;
		star.vx = sqrt(v*v-star.vz*star.vz)*cos(x7);
		star.vy = sqrt(v*v-star.vz*star.vz)*sin(x7);
		
		star.x *= 3.*M_PI/64.*M*M/E;
		star.y *= 3.*M_PI/64.*M*M/E;
		star.z *= 3.*M_PI/64.*M*M/E;
		
		star.vx *= sqrt(E*64./3./M_PI/M);
		star.vy *= sqrt(E*64./3./M_PI/M);
		star.vz *= sqrt(E*64./3./M_PI/M);

		star.m = M/(double)_N;


		particles_add(star);
	}
	tools_move_to_center_of_momentum();

	printf("Characteristic size:              %f\n", r0);
	printf("Characteristic time (relaxation): %f\n", t0);
}

void problem_inloop(){
}

double calculate_energy(struct particle* _particles, int _N){
	double energy_kinetic = 0;
	double energy_potential = 0;
	for (int i=0;i<_N;i++){
		for (int j=0;j<i;j++){
			double dx = _particles[i].x - particles[j].x;
			double dy = _particles[i].y - particles[j].y;
			double dz = _particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy_potential += -G*_particles[i].m*_particles[j].m/r;
		}
		double dvx = _particles[i].vx;
		double dvy = _particles[i].vy;
		double dvz = _particles[i].vz;
		energy_kinetic += 1./2. * _particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	return energy_potential + energy_kinetic;
}
int compare (const void * a, const void * b){
	if (*(double*)a < *(double*)b ) return -1;
	if (*(double*)a > *(double*)b ) return 1;
	return 0;
}

void output_radii(){
	tools_move_to_center_of_momentum();
	double* radii = malloc(sizeof(double)*N);
	for (int i=0;i<N;i++){
		double r = sqrt(particles[i].x*particles[i].x + particles[i].y*particles[i].y + particles[i].z*particles[i].z);
		radii[i] = r;
	}
	qsort (radii, N, sizeof(double), compare);
	
	FILE* of = fopen("radii.txt","a+"); 
	fprintf(of,"%e\t",t/t0);
	int index10 = (int)ceil(0.1*(double)N);
	if (index10>=N) index10 = N-1;
	int index50 = (int)ceil(0.5*(double)N);
	if (index50>=N) index50 = N-1;
	int index90 = (int)ceil(0.9*(double)N);
	if (index90>=N) index90 = N-1;
	fprintf(of,"%e\t",radii[index10]/r0);
	fprintf(of,"%e\t",radii[index50]/r0);
	fprintf(of,"%e\t",radii[index90]/r0);
	
//	for (int i=0;i<N;i++){
//		fprintf(of,"%e\t%d\t%e\n",t/t0,i,radii[i]/r0);
//	}

	fprintf(of,"\n");
	fclose(of);
}


void calculate_error(){
	FILE* inf = fopen("best.bin","rb"); 
	int _N;
	double _t;
	fread(&_N,sizeof(int),1,inf);
	fread(&_t,sizeof(double),1,inf);
	struct particle* _particles = malloc(_N*sizeof(struct particle));
	for (int i=0;i<_N;i++){
		fread(&(_particles[i]),sizeof(struct particle),1,inf);
	}
	fclose(inf);
	double dif_pos = 0;
	double dif_posabs = 0;
	double dif_vel = 0;
	for (int i=0;i<N;i++){
		double dx = _particles[i].x - particles[i].x;
		double dy = _particles[i].y - particles[i].y;
		double dz = _particles[i].z - particles[i].z;
		dif_pos += sqrt(dx*dx + dy*dy + dz*dz )/(double)N;
		dif_posabs += fabs(dx)/(double)N;
		dif_posabs += fabs(dy)/(double)N;
		dif_posabs += fabs(dz)/(double)N;
		double dvx = _particles[i].vx - particles[i].vx;
		double dvy = _particles[i].vy - particles[i].vy;
		double dvz = _particles[i].vz - particles[i].vz;
		dif_vel += sqrt(dvx*dvx + dvy*dvy + dvz*dvz )/(double)N;
	}
	double energy_best = calculate_energy(_particles,_N);
	double energy      = calculate_energy(particles,N);
	double dif_energy = fabs((energy-energy_best)/energy_best);
#ifdef INTEGRATOR_LEAPFROG
	FILE* of = fopen("error_leapfrog.txt","a+"); 
#endif
#ifdef INTEGRATOR_EULER
	FILE* of = fopen("error_euler.txt","a+"); 
#endif
#ifdef INTEGRATOR_MODIFIEDEULER
	FILE* of = fopen("error_modifiedeuler.txt","a+"); 
#endif
#ifdef INTEGRATOR_RADAU15
	FILE* of = fopen("error_radau15.txt","a+"); 
#endif
#ifdef INTEGRATOR_IMPULSE
	FILE* of = fopen("error_impulse.txt","a+"); 
#endif
#ifdef INTEGRATOR_IMPULSE2NDORDER
	FILE* of = fopen("error_impulse2ndorder.txt","a+"); 
#endif
	fprintf(of,"%e\t",dt);
	fprintf(of,"%e\t",dif_pos);
	fprintf(of,"%e\t",dif_vel);
	fprintf(of,"%e\t",dif_energy);
	fprintf(of,"%e\t",dif_posabs);
	fprintf(of,"\n");
	fclose(of);
	FILE* ofp = fopen("position.txt","a+"); 
	fprintf(ofp,"%e\t%e\t%e\t%e\n",dt,particles[0].x,particles[0].y,particles[0].z);
	fclose(ofp);
	
}

struct vec3o {
	double x;
	double y;
	double z;
};


void problem_finish(){
	printf("\nFinished. Time: %e. Difference to exact time: %e. Timestep: %e\n",t,t-tmax,dt);
	if (writeBest){
		output_binary("best.bin");
	}
	calculate_error();
}
void problem_output(){
	if (output_check(10.0*dt)) output_timing();
	if (output_check(tmax/1000.)) output_radii();
}

