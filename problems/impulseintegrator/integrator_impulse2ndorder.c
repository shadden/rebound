/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	 
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


struct xyz* xp; 
struct xyz* xpp;
struct xyz* vp;
struct xyz* vpp;
long _arraymax = 0;
bool oddeven = 0;
void integrator_part1(){
	if (_arraymax<N){
		xp  = realloc(xp, sizeof(struct xyz)*N);
		xpp = realloc(xpp,sizeof(struct xyz)*N);
		vp  = realloc(vp, sizeof(struct xyz)*N);
		vpp  = realloc(vp, sizeof(struct xyz)*N);
		_arraymax = N;
	}
	oddeven = !oddeven;
	if (oddeven){
		// x prime
		for (int i=0;i<N;i++){
			// Initial guess
			xp[i].x = particles[i].x+dt*particles[i].vx;
			xp[i].y = particles[i].y+dt*particles[i].vy;
			xp[i].z = particles[i].z+dt*particles[i].vz;
		}
		// Do 5 iterations. Let's hope we're converged.
		int iterations_N = 5;
		for (int iterations=0; iterations<iterations_N; iterations++){
			for (int i=0;i<N;i++){
				xpp[i].x = particles[i].x+dt*particles[i].vx;
				xpp[i].y = particles[i].y+dt*particles[i].vy;
				xpp[i].z = particles[i].z+dt*particles[i].vz;
				for (int j=0;j<N;j++){
					if (i!=j){
						struct xyz _xp;
						_xp.x = xp[i].x - xp[j].x;
						_xp.y = xp[i].y - xp[j].y;
						_xp.z = xp[i].z - xp[j].z;
						struct xyz _v;
						_v.x = particles[i].vx - particles[j].vx;
						_v.y = particles[i].vy - particles[j].vy;
						_v.z = particles[i].vz - particles[j].vz;
						double _v2 	= _v.x*_v.x + _v.y*_v.y + _v.z*_v.z; 
						double _v1  	= sqrt(_v2); 
						double _v3 	= _v2*_v1;
						double y2 	= _xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening;
						double y1 	= sqrt(y2);
						double _xp_v 	= _xp.x*_v.x  + _xp.y*_v.y  + _xp.z*_v.z;
						{	// tau = dt
							struct xyz _xp_v_tau;
							_xp_v_tau.x = _xp.x - _v.x*dt;
							_xp_v_tau.y = _xp.y - _v.y*dt;
							_xp_v_tau.z = _xp.z - _v.z*dt;
							double s = sqrt(_xp_v_tau.x*_xp_v_tau.x + _xp_v_tau.y*_xp_v_tau.y + _xp_v_tau.z*_xp_v_tau.z + softening*softening);
							double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));
							

							// Calculating the new x prime
							// First term in bracket
							xpp[i].x += prefactor * _xp.x * (y2 - _xp_v*dt);
							xpp[i].y += prefactor * _xp.y * (y2 - _xp_v*dt);
							xpp[i].z += prefactor * _xp.z * (y2 - _xp_v*dt);
							
							// Second term in bracket
							double numerator = y2*_xp_v + _v2*y2*dt - 2.*dt*_xp_v*_xp_v;
							xpp[i].x -= prefactor * _v.x / _v2 * numerator;  
							xpp[i].y -= prefactor * _v.y / _v2 * numerator;
							xpp[i].z -= prefactor * _v.z / _v2 * numerator;
							
							// Third term in bracket
							double _log = log((_v2*dt-_xp_v+_v1*s) / (_v1*y1-_xp_v));
							xpp[i].x += G*particles[j].m*_v.x/_v3 * _log;
							xpp[i].y += G*particles[j].m*_v.y/_v3 * _log;
							xpp[i].z += G*particles[j].m*_v.z/_v3 * _log;
						}
						{	// tau = 0
							double s = sqrt(_xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening);
							double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));


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
							double _log = log((-_xp_v+_v1*s) / (_v1*y1-_xp_v));
							xpp[i].x -= G*particles[j].m*_v.x/_v3 * _log;
							xpp[i].y -= G*particles[j].m*_v.y/_v3 * _log;
							xpp[i].z -= G*particles[j].m*_v.z/_v3 * _log;
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

		// v prime
		for (int i=0;i<N;i++){
			vp[i].x = particles[i].vx;
			vp[i].y = particles[i].vy;
			vp[i].z = particles[i].vz;
			for (int j=0;j<N;j++){
				if (i!=j){
					struct xyz _xp;
					_xp.x = xp[i].x - xp[j].x;
					_xp.y = xp[i].y - xp[j].y;
					_xp.z = xp[i].z - xp[j].z;
					struct xyz _v;
					_v.x = particles[i].vx - particles[j].vx;
					_v.y = particles[i].vy - particles[j].vy;
					_v.z = particles[i].vz - particles[j].vz;
					double _v2 	= _v.x*_v.x + _v.y*_v.y + _v.z*_v.z; 
					double _v1  	= sqrt(_v2); 
					double _v3 	= _v2*_v1;
					double y2 	= _xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening;
					double _xp_v 	= _xp.x*_v.x  + _xp.y*_v.y  + _xp.z*_v.z;
					{	// tau = dt
						struct xyz _xp_v_tau;
						_xp_v_tau.x = _xp.x - _v.x*dt;
						_xp_v_tau.y = _xp.y - _v.y*dt;
						_xp_v_tau.z = _xp.z - _v.z*dt;
						double s = sqrt(_xp_v_tau.x*_xp_v_tau.x + _xp_v_tau.y*_xp_v_tau.y + _xp_v_tau.z*_xp_v_tau.z + softening*softening);
						double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));
						vp[i].x -= prefactor*(_v.x*(y2-_xp_v*dt) + _xp.x*(_v2*dt-_xp_v));
						vp[i].y -= prefactor*(_v.y*(y2-_xp_v*dt) + _xp.y*(_v2*dt-_xp_v));
						vp[i].z -= prefactor*(_v.z*(y2-_xp_v*dt) + _xp.z*(_v2*dt-_xp_v));
					}

					{	// tau = 0
						double s = sqrt(_xp.x*_xp.x + _xp.y*_xp.y + _xp.z*_xp.z + softening*softening);
						double prefactor = G*particles[j].m/(s*(_v2*y2-(_xp_v*_xp_v)));
						vp[i].x += prefactor*(_v.x*(y2) + _xp.x*(-_xp_v));
						vp[i].y += prefactor*(_v.y*(y2) + _xp.y*(-_xp_v));
						vp[i].z += prefactor*(_v.z*(y2) + _xp.z*(-_xp_v));
					}
				}
			}
		}
							
		
		// Set variables for next iteration;
		for (int i=0;i<N;i++){
			particles[i].x   = xp[i].x;
			particles[i].y   = xp[i].y;
			particles[i].z   = xp[i].z;
			particles[i].vx  = vp[i].x;
			particles[i].vy  = vp[i].y;
			particles[i].vz  = vp[i].z;
		}
	}else{
		for (int i=0;i<N;i++){
			// Initial guess
			vp[i].x = particles[i].vx;
			vp[i].y = particles[i].vy;
			vp[i].z = particles[i].vz;
		}
		// Do 5 iterations. Let's hope we're converged.
		int iterations_N = 5;
		for (int iterations=0; iterations<iterations_N; iterations++){
			for (int i=0;i<N;i++){
				vpp[i].x = particles[i].vx;
				vpp[i].y = particles[i].vy;
				vpp[i].z = particles[i].vz;
				for (int j=0;j<N;j++){
					if (i!=j){
						struct xyz _x;
						_x.x = particles[i].x - particles[j].x;
						_x.y = particles[i].y - particles[j].y;
						_x.z = particles[i].z - particles[j].z;
						struct xyz _vp;
						_vp.x = vp[i].vx - vp[j].vx;
						_vp.y = vp[i].vy - vp[j].vy;
						_vp.z = vp[i].vz - vp[j].vz;
						double _vp2 	= _vp.x*_vp.x + _vp.y*_vp.y + _vp.z*_vp.z; 
						double _vp1  	= sqrt(_vp2); 
						double _vp3 	= _vp2*_vp1;
						double yb2 	= _x.x*_x.x + _x.y*_x.y + _x.z*_x.z + softening*softening;
						double _x_vp 	= _x.x*_vp.x  + _x.y*_vp.y  + _x.z*_vp.z;
						{	// tau = dt
							struct xyz _x_vp_tau;
							_x_vp_tau.x = _x.x + _vp.x*dt;
							_x_vp_tau.y = _x.y + _vp.y*dt;
							_x_vp_tau.z = _x.z + _vp.z*dt;
							double sb = sqrt(_x_vp_tau.x*_x_vp_tau.x + _x_vp_tau.y*_x_vp_tau.y + _x_vp_tau.z*_x_vp_tau.z + softening*softening);
							double prefactor = G*particles[j].m/(sb*(_vp2*yb2-(_x_vp*_x_vp)));
							
							// Calculating the new v prime
							vpp[i].x += prefactor*(_vp.x*(yb2+_x_vp*dt) - _x.x*(_vp2*dt+_x_vp));
							vpp[i].y += prefactor*(_vp.y*(yb2+_x_vp*dt) - _x.y*(_vp2*dt+_x_vp));
							vpp[i].z += prefactor*(_vp.z*(yb2+_x_vp*dt) - _x.z*(_vp2*dt+_x_vp));


						}
						{	// tau = 0
							double sb = sqrt(_x.x*_x.x + _x.y*_x.y + _x.z*_x.z + softening*softening);
							double prefactor = G*particles[j].m/(sb*(_vp2*yb2-(_x_vp*_x_vp)));

							// Calculating the new v prime
							vpp[i].x -= prefactor*(_vp.x*(yb2) + _x.x*(_x_vp));
							vpp[i].y -= prefactor*(_vp.y*(yb2) + _x.y*(_x_vp));
							vpp[i].z -= prefactor*(_vp.z*(yb2) + _x.z*(_x_vp));
						}
					}
				}
			}
			for (int i=0;i<N;i++){
				// Ready for next iteration;
				vp[i].x = vpp[i].x;
				vp[i].y = vpp[i].y;
				vp[i].z = vpp[i].z;
			}
		}
		// x prime
		for (int i=0;i<N;i++){
			xp[i].x = particles[i].x + dt*vp[i].x;
			xp[i].y = particles[i].y + dt*vp[i].y;
			xp[i].z = particles[i].z + dt*vp[i].z;
			for (int j=0;j<N;j++){
				if (i!=j){
					struct xyz _x;
					_x.x = particles[i].x - particles[j].x;
					_x.y = particles[i].y - particles[j].y;
					_x.z = particles[i].z - particles[j].z;
					struct xyz _vp;
					_vp.x = vp[i].vx - vp[j].vx;
					_vp.y = vp[i].vy - vp[j].vy;
					_vp.z = vp[i].vz - vp[j].vz;
					double _vp2 	= _vp.x*_vp.x + _vp.y*_vp.y + _vp.z*_vp.z; 
					double _vp1  	= sqrt(_vp2); 
					double _vp3 	= _vp2*_vp1;
					double yb2 	= _x.x*_x.x + _x.y*_x.y + _x.z*_x.z + softening*softening;
					double yb1 	= sqrt(yb2);
					double _x_vp 	= _x.x*_vp.x  + _x.y*_vp.y  + _x.z*_vp.z;
					{	// tau = dt
						struct xyz _x_vp_tau;
						_x_vp_tau.x = _x.x + _vp.x*dt;
						_x_vp_tau.y = _x.y + _vp.y*dt;
						_x_vp_tau.z = _x.z + _vp.z*dt;
						double sb = sqrt(_x_vp_tau.x*_x_vp_tau.x + _x_vp_tau.y*_x_vp_tau.y + _x_vp_tau.z*_x_vp_tau.z + softening*softening);
						double prefactor = -G*particles[j].m/(sb*(_vp2*yb2-(_x_vp*_x_vp)));
						
						// First term in bracket
						xp[i].x += prefactor * _x.x * (yb2 + _x_vp*dt);
						xp[i].y += prefactor * _x.y * (yb2 + _x_vp*dt);
						xp[i].z += prefactor * _x.z * (yb2 + _x_vp*dt);
						
						// Second term in bracket
						double numerator = yb2*_x_vp - _vp2*yb2*dt + 2.*dt*_x_vp*_x_vp;
						xp[i].x -= prefactor * _vp.x / _vp2 * numerator; 
						xp[i].y -= prefactor * _vp.y / _vp2 * numerator;
						xp[i].z -= prefactor * _vp.z / _vp2 * numerator;
						
						// Third term in bracket
						double _log = log((_vp2*dt-_x_vp+_vp1*sb) / (_vp1*yb1-_x_vp));
						xp[i].x += G*particles[j].m*_vp.x/_vp3 * _log;
						xp[i].y += G*particles[j].m*_vp.y/_vp3 * _log;
						xp[i].z += G*particles[j].m*_vp.z/_vp3 * _log;
					}
					{	// tau = 0
						double sb = sqrt(_x.x*_x.x + _x.y*_x.y + _x.z*_x.z + softening*softening);
						double prefactor = -G*particles[j].m/(sb*(_vp2*yb2-(_x_vp*_x_vp)));

						// First term in bracket
						xp[i].x -= prefactor * _x.x * (yb2);
						xp[i].y -= prefactor * _x.y * (yb2);
						xp[i].z -= prefactor * _x.z * (yb2);
						
						// Second term in bracket
						double numerator = yb2*_x_vp;
						xp[i].x += prefactor * _vp.x / _vp2 * numerator; 
						xp[i].y += prefactor * _vp.y / _vp2 * numerator;
						xp[i].z += prefactor * _vp.z / _vp2 * numerator;
						
						// Third term in bracket
						double _log = log((-_x_vp+_vp1*sb) / (_vp1*yb1-_x_vp));
						xp[i].x -= G*particles[j].m*_v.x/_v3 * _log;
						xp[i].y -= G*particles[j].m*_v.y/_v3 * _log;
						xp[i].z -= G*particles[j].m*_v.z/_v3 * _log;
					}
				}
			}
		}

		// Set variables for next iteration;
		for (int i=0;i<N;i++){
			particles[i].x   = xp[i].x;
			particles[i].y   = xp[i].y;
			particles[i].z   = xp[i].z;
			particles[i].vx  = vp[i].x;
			particles[i].vy  = vp[i].y;
			particles[i].vz  = vp[i].z;
		}
	}

	t+=dt;
}
void integrator_part2(){
}
	

