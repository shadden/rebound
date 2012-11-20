/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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

struct xyz {
	double x;
	double y;
	double z;
};

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.

struct xyz* x; 
struct xyz* v; 
struct xyz* xp; 
struct xyz* xpp;
struct xyz* vpp;
long _arraymax = 0;
void integrator_part1(){
	if (_arraymax<N){
		x   = realloc(x,  sizeof(struct xyz)*N);
		v   = realloc(v,  sizeof(struct xyz)*N);
		xp  = realloc(xp, sizeof(struct xyz)*N);
		xpp = realloc(xpp,sizeof(struct xyz)*N);
		vpp = realloc(vpp,sizeof(struct xyz)*N);
		_arraymax = N;
	}
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		// Initial conditions
		x[i].x = particles[i].x;
		x[i].y = particles[i].y;
		x[i].z = particles[i].z;
		v[i].x = particles[i].vx;
		v[i].y = particles[i].vy;
		v[i].z = particles[i].vz;
		// Initial guess
		xp[i].x = x[i].x+dt*v[i].x;
		xp[i].y = x[i].y+dt*v[i].y;
		xp[i].z = x[i].z+dt*v[i].z;
	}
	// Do 5 iterations. Let's hope we're converged.
	for (int iterations=0;iterations<5;iterations++){
		for (int i=0;i<N;i++){
			xpp[i].x = x[i].x+dt*v[i].x;
			xpp[i].y = x[i].y+dt*v[i].y;
			xpp[i].z = x[i].z+dt*v[i].z;
			vpp[i].x = v[i].x;
			vpp[i].y = v[i].y;
			vpp[i].z = v[i].z;
			for (int j=0;j<N;j++){
				if (i!=j){
					struct xyz _xp;
					_xp.x = xp[i].x - xp[j].x;
					_xp.y = xp[i].y - xp[j].y;
					_xp.z = xp[i].z - xp[j].z;
					struct xyz _v;
					_v.x = v[i].x - v[j].x;
					_v.y = v[i].y - v[j].y;
					_v.z = v[i].z - v[j].z;
					double _v2 = _v.x*_v.x + _v.y*_v.y + _v.z*_v.z; 
					double _v1  = sqrt(_v2); 
					double _v3 = _v1*_v1*_v1;
					double y2 = _xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening;
					double _xp_v = _xp.x*_v.x+_xp.y*_v.y+_xp.z*_v.z;
					{	// tau = dt
						struct xyz _xp_v_tau;
						_xp_v_tau.x = _xp.x-_v.x*dt;
						_xp_v_tau.y = _xp.y-_v.y*dt;
						_xp_v_tau.z = _xp.z-_v.z*dt;
						double s = sqrt(_xp_v_tau.x*_xp_v_tau.x + _xp_v_tau.y*_xp_v_tau.y + _xp_v_tau.z*_xp_v_tau.z + softening*softening);
						double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));
						
						// Calculating the new v prime
						vpp[i].x += -prefactor*(_v.x*(y2-_xp_v*dt) + _xp.x*(_v2*dt-_xp_v));
						vpp[i].y += -prefactor*(_v.y*(y2-_xp_v*dt) + _xp.y*(_v2*dt-_xp_v));
						vpp[i].z += -prefactor*(_v.z*(y2-_xp_v*dt) + _xp.z*(_v2*dt-_xp_v));


						// Calculating the new x prime
						// First term in bracket
						xpp[i].x += prefactor * _xp.x * (y2-_xp_v*dt);
						xpp[i].y += prefactor * _xp.y * (y2-_xp_v*dt);
						xpp[i].z += prefactor * _xp.z * (y2-_xp_v*dt);
						
						// Second term in bracket
						double numerator = y2*_xp_v + _v2*y2*dt-2.*dt*_xp_v*_xp_v;
						xpp[i].x -= prefactor * _v.x / _v2 * numerator;
						xpp[i].y -= prefactor * _v.y / _v2 * numerator;
						xpp[i].z -= prefactor * _v.z / _v2 * numerator;
						
						// Third term in bracket
						double _log = log(_v2*dt-_xp_v+_v1*s);
						xpp[i].x -= G*particles[j].m*_v.x/_v3 * _log;
						xpp[i].y -= G*particles[j].m*_v.y/_v3 * _log;
						xpp[i].z -= G*particles[j].m*_v.z/_v3 * _log;
					}
					{	// tau = 0
						double s = sqrt(_xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening);
						double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));

						// Calculating the new v prime
						vpp[i].x -= -prefactor*(_v.x*(y2) + _xp.x*(-_xp_v));
						vpp[i].y -= -prefactor*(_v.y*(y2) + _xp.y*(-_xp_v));
						vpp[i].z -= -prefactor*(_v.z*(y2) + _xp.z*(-_xp_v));
						

						// Calculating the new x prime
						// First term in bracket
						xpp[i].x -= prefactor * _xp.x * (y2);
						xpp[i].y -= prefactor * _xp.y * (y2);
						xpp[i].z -= prefactor * _xp.z * (y2);
						
						// Second term in bracket
						double numerator = y2*_xp_v;
						xpp[i].x += prefactor * _v.x / _v2 * numerator;
						xpp[i].y += prefactor * _v.y / _v2 * numerator;
						xpp[i].z += prefactor * _v.z / _v2 * numerator;
						
						// Third term in bracket
						double _log = log(-_xp_v+_v1*s);
						xpp[i].x += G*particles[j].m*_v.x/_v3 * _log;
						xpp[i].y += G*particles[j].m*_v.y/_v3 * _log;
						xpp[i].z += G*particles[j].m*_v.z/_v3 * _log;
					}
				}
			}
		}
		for (int i=0;i<N;i++){
			// Ready for next iteration;
			xp[i].x = xpp[i].x;
			xp[i].y = xpp[i].y;
			xp[i].z = xpp[i].z;
		}
	}
	
	// Set variables for next iteration;
	for (int i=0;i<N;i++){
		particles[i].x  = xp[i].x;
		particles[i].y  = xp[i].y;
		particles[i].z  = xp[i].z;
		particles[i].vx = vpp[i].x;
		particles[i].vy = vpp[i].y;
		particles[i].vz = vpp[i].z;
	}

	t+=dt;
}
void integrator_part2(){
}
	

