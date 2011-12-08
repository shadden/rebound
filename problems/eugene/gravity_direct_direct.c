/**
 * @file 	gravity.c
 * @brief 	Double Direct gravity calculation, O(N).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * 
 * @section LICENSE
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
#include "boundaries.h"
#include "communication_mpi.h"

#ifdef MPI
#warning GRAVITY_DIRECT_DIRECT may not work as expected for your problem.
#warning Make sure you understand what the code is doing.
#endif

void gravity_calculate_acceleration(){
#pragma omp parallel for schedule(guided)
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	if (N_active == -1){
		printf("Error: N_active needs to be set for GRAVITY_DIRECT_DIRECT module.\n");
		exit(-1);
	}
	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
	for (int gbz=-nghostz; gbz<=nghostz; gbz++){
		struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
		// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
		for (int i=0; i<N_active; i++){
		for (int j=N_active; j<N; j++){
			double dx = (gb.shiftx+particles[i].x) - particles[j].x;
			double dy = (gb.shifty+particles[i].y) - particles[j].y;
			double dz = (gb.shiftz+particles[i].z) - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			double prefact = -G/(r*r*r)*particles[j].m;
			particles[i].ax += prefact*dx; 
			particles[i].ay += prefact*dy; 
			particles[i].az += prefact*dz; 
		}
		}
#pragma omp parallel for schedule(guided)
		for (int i=N_active; i<N; i++){
		for (int j=0; j<N_active; j++){
			double dx = (gb.shiftx+particles[i].x) - particles[j].x;
			double dy = (gb.shifty+particles[i].y) - particles[j].y;
			double dz = (gb.shiftz+particles[i].z) - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			double prefact = -G/(r*r*r)*particles[j].m;
			particles[i].ax += prefact*dx; 
			particles[i].ay += prefact*dy; 
			particles[i].az += prefact*dz; 
		}
		}
	}
	}
	}
#ifdef MPI
	// Distribute active particles from root to all other nodes.
	// This assures that round-off errors do not accumulate and 
	// the copies of active particles diverge. 
	for(int i=0;i<N_active;i++){
		struct particle p = particles[i];
		// Sum acceleration
		double a_send[3];
		a_send[0] = p.ax; a_send[1] = p.ay; a_send[2] = p.az; 
		double a_recv[3] = {0,0,0};
		MPI_Reduce(a_send,a_recv,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		particles[i].ax = a_recv[0];
		particles[i].ay = a_recv[1];
		particles[i].az = a_recv[2];
		
		// Average velocity (due to collisions)
		double v_send[3];
		v_send[0] = p.vx; v_send[1] = p.vy; v_send[2] = p.vz; 
		double v_recv[3] = {0,0,0};
		MPI_Reduce(v_send,v_recv,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		particles[i].vx = v_recv[0]/(double)mpi_num;
		particles[i].vy = v_recv[1]/(double)mpi_num;
		particles[i].vz = v_recv[2]/(double)mpi_num;
	}
	MPI_Bcast(particles, N_active, mpi_particle, 0, MPI_COMM_WORLD); 
#endif
}
