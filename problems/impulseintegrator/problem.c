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
void output_statistics();
double t0, r0;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	dt		= input_get_double(argc,argv,"dt",0.0001);
	writeBest	= input_get_int(argc,argv,"write",0);

	
	int _N 		= input_get_int(argc,argv,"N",4096);
	double M = 1;
	double R = 64./3./M_PI;
	double E = 1;
	r0 = 16./(3.*M_PI)*R;	// 1.6976527 * R
	t0 = G*pow(M,5./2.)*pow(4.*E,-3./2.)*(double)_N/log(0.4*(double)_N); //Rellaxation time
	tmax		= 2.*t0;
	boxsize 	= 100.*r0;
	softening 	= 0.025;		

	init_box();
	
	
	
	tools_init_plummer(_N, M, R);
	tools_move_to_center_of_momentum();

	output_statistics();
	exit(1);
}

void problem_inloop(){
}
double calculate_mininteractiontime(){
	double mininteractiontime = 1e12;
	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			double dvx = particles[i].vx - particles[j].vx;
			double dvy = particles[i].vy - particles[j].vy;
			double dvz = particles[i].vz - particles[j].vz;
			double v = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
			double interactiontime = fabs(r/v);
			if (interactiontime<mininteractiontime){
				mininteractiontime = interactiontime;
			}
		}
	}
	return mininteractiontime;
}

struct vec3 {
	double x;
	double y;
	double z;
};

void output_statistics(){
	printf("\n\n");
	printf("\t Number of particles: \t%d\n",N);
	double mass = 0;
	for (int i=0;i<N;i++){
		mass += particles[i].m;
	}
	printf("\t Total mass:           \t%e\n",mass);
	// Algorithm with reduced roundoff errors (see wikipedia)
	struct vec3 A = {.x=0, .y=0, .z=0};
	struct vec3 Q = {.x=0, .y=0, .z=0};
	for (int i=0;i<N;i++){
		struct vec3 Aim1 = A;
		struct particle p = particles[i];
		A.x = A.x + (p.vx-A.x)/(double)(i+1);
		A.y = A.y + (p.vy-A.y)/(double)(i+1);
		A.z = A.z + (p.vz-A.z)/(double)(i+1);
		Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
		Q.y = Q.y + (p.vy-Aim1.y)*(p.vy-A.y);
		Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
	}
	Q.x = sqrt(Q.x/(double)N);
	Q.y = sqrt(Q.y/(double)N);
	Q.z = sqrt(Q.z/(double)N);
	double totQ = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
	printf("\t Velocity dispersion: \t%e\n",totQ);

	double energy_kinetic = 0;
	double energy_potential = 0;
	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			energy_potential += -G*particles[i].m*particles[j].m/r;
		}
		double dvx = particles[i].vx;
		double dvy = particles[i].vy;
		double dvz = particles[i].vz;
		energy_kinetic += 1./2. * particles[i].m* (dvx*dvx + dvy*dvy + dvz*dvz);
	}
	
	printf("\t Energy (kinetic): \t%e\n",energy_kinetic);
	printf("\t Energy (potential): \t%e\n",energy_potential);
	printf("\t Energy (total): \t%e\n",energy_potential+energy_kinetic);
	printf("\n\n");
}

double calculate_energy(struct particle* _particles, int _N){
	double energy_kinetic = 0;
	double energy_potential = 0;
	for (int i=0;i<_N;i++){
		for (int j=0;j<i;j++){
			double dx = _particles[i].x - _particles[j].x;
			double dy = _particles[i].y - _particles[j].y;
			double dz = _particles[i].z - _particles[j].z;
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
	//printf("Min interaction time: %e\n",calculate_mininteractiontime());
	if (output_check(10.0*dt)) output_timing();
	if (output_check(tmax/1000.)) output_radii();
}

