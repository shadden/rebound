/**
 * @file 	problem.c
 * @brief 	Example problem: Cometary Outgassing + general relativity correction
 * @author 	Ridlo W. Wibowo <ridlo.w.wibowo@gmail.com>
 * 
 * @detail 	This example applies outgassing parameter of comet
 * and general relativity correction as additional forces
 * Special thanks goes to Dr. Budi Dermawan
 *
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
#include "main.h"
#include "tools.h"
#include "problem.h"
#include "output.h"
#include "particle.h"
#include "collisions.h"
#include "collision_resolve.h"
#include "boundaries.h"

// Input ephemeris
// JD 2456765.5
// April 18th 2014, HORIZON JPL NASA
// mass (M_sun), equatorial radius(AU), x(AU), y, z, vx(AU/day), vy, vz (barycentric)
double obj[14][8] = 
{
    {1.0, 0.005, 1.605553422646323E-03, -1.968808432193491E-03, -1.074011886583918E-04, 5.463278467466615E-06, 3.343531001090643E-06, -1.278227982325224E-07}, // Sun
    {1.660611985452673e-07, 1.63103926e-5, 3.605757382610923E-01, -4.226381789171472E-02, -3.633487676150426E-02, -2.296900630182340E-03,  2.922170498903202E-02,  2.598469713284350E-03}, // Mercury
    {2.4482737118213123e-06, 4.04537843e-5, -4.362607855867268E-02, -7.272949739012508E-01, -7.437644417360898E-03, 2.005685837037687E-02, -1.331643224431570E-03, -1.175587369983879E-03}, // Venus
    {3.003297890315729e-06, 4.26352325e-5, -8.873685319075689E-01, -4.684852415256254E-01, -9.151929903703799E-05, 7.714218411139476E-03, -1.529183485754389E-02, -2.099000065353641E-07}, // Earth
    {3.695668790833897e-08, 1.16146707e-5, -8.885826964904660E-01, -4.706836245078943E-01,  3.474194046391472E-05, 8.251454771541669E-03, -1.556573935321515E-02,  4.559225559678975E-05}, // Moon
    {3.2277384860480837e-07, 2.26602156e-5, -1.483812590735944E+00, -6.315488622947584E-01,  2.316158310228655E-02, 5.990331667260308E-03, -1.168454224213042E-02, -3.919289618281470E-04}, // Mars
    {0.0009545325625181037, 0.000477894503, -2.100997440426695E+00,  4.788335671772308E+00,  2.704703479214598E-02, -7.001353924993028E-03, -2.675649628445124E-03,  1.677745431632356E-04}, // Jupiter
    {0.00028579654259598984, 0.000402866697, -6.476977890157666E+00, -7.482849654255215E+00,  3.878787094103098E-01, 3.913238296648806E-03, -3.665932303961103E-03, -9.190127966157570E-05}, // Saturn
    {4.365520702584404e-05, 0.000170851362, 1.955642937590881E+01,  4.312556734330099E+00, -2.373436786435929E-01, -8.756791019987000E-04,  3.657502374146099E-03,  2.501797114914699E-05}, // Uranus
    {5.149999195391201e-05, 0.000165550485, 2.720496360818385E+01, -1.258810690688079E+01, -3.677377161902303E-01, 1.297241662899291E-03,  2.867871651305701E-03, -8.855065105036448E-05}, // Neptune
    {4.762176023591137e-10, 3.18320039e-6, -2.341996894144388E+00, -1.126774810951099E+00,  3.966527787360209E-01, 3.961215501996468E-03, -1.006503132698686E-02, -1.045852992933018E-03}, // Ceres
    {1.3412457787962378e-10, 1.77141559e-6, -1.998762960534205E+00, -9.206966932468410E-01,  2.707538615296253E-01, 5.621078135077726E-03, -1.047820653128300E-02, -3.688099589606132E-04}, // Vesta    
    {1.0775176762239439e-10, 1.82154999e-6, -2.200952591440287E+00,  4.472974388461693E-01, -1.262924415683855E-01, -4.098052054240473E-03, -9.542041961277441E-03,  6.939213710195355E-03}, // Pallas
    {0.0, 0.0, -3.572790037832271E+00,  3.528915544836181E+00,  7.857477933955778E-02, -2.019119936704946E-03,  5.414338544696132E-04, -8.703395238496269E-05} // SOHO (C/1999 X3)
};

void outgassing_correction_force();
void relativity_correction_force();
void additional_force();

// const double c  = 173.144632674; // speed of light in AU/day
const double _c2 = 3.335661199685737e-05; // 1/c^2
const double k   = 0.01720209895; // Gaussian constant

#ifdef OPENGL
    extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
    // Setup integration parameters
    dt          = 0.05; // days
    tmax        = 3.652422e6; // 10^4 yr
    N_active    = 13;  // "massive" object
    G           = k*k;

#ifdef OPENGL
    display_wire = 1; // show orbits
#endif

    init_boxwidth(200); // init box with width 200 AU (200 x 200 x 200)

    // Initial conditions
	for (int i=0;i<14;i++){
		struct particle p;
        p.m = obj[i][0];
        // p.r = obj[i][1]; // not needed for collisions_none.c
		p.x  = obj[i][2]; p.y = obj[i][3]; p.z = obj[i][4];
        p.vx = obj[i][5]; p.vy = obj[i][6]; p.vz = obj[i][7];
        p.ax = 0.0; p.ay = 0.0; p.az = 0.0;
		particles_add(p); 
	}

    problem_additional_forces = additional_force;

	system("rm -f orbits.txt");
}

void outgassing_correction_force(){
    // for comet
    struct particle sun = particles[0]; // sun -> copy of parameter
    struct particle* com = &(particles[13]); // comet -> real parameter
    
    // Cometary outgassing componenets is not available for  C/1999 X3    
    // we use parameters from: 96P/Machholz 1 [2013]
    const double A1 = 1.069808937609E-10; // in AU/day^2
    const double A2 = -9.383180440636E-13;
    const double A3 = 0.0;

    // Yeomans et al. 2004
    // alpha = 0.111262;
    // ro = 2.808;
    // m = 2.15;
    // n = 5.093;
    // k = 4.6142;     
    
    double dx = com->x - sun.x; // ! calculate 2x for same object (in relativity and outgassing)
    double dy = com->y - sun.y;
    double dz = com->z - sun.z;
    
    double dvx = com->vx - sun.vx;
    double dvy = com->vy - sun.vy;
    double dvz = com->vz - sun.vz;
    
    double r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 < 88.){ // if body close to the Sun (R < 9.36 AU)
        double r = sqrt(dx*dx + dy*dy + dz*dz); 
        double q = r*0.3561253561253561; // Q = R/R0, R0 = 2.808 AU
        //gr = alpha * pow(Q, -m) * pow((1 + pow(Q, n)), -k);
        double gr = 0.111262 * pow(q, -2.15) * pow((1.0 + pow(q, 5.093)), -4.6142);
        
        /* tangensial term */
        // tangensial vector = (r x v) x r = v(r . r) - r(r . v)
        double rv = dx*dvx + dy*dvy + dz*dvz;
        double tx = r2*dvx - rv*dx;
        double ty = r2*dvy - rv*dy;
        double tz = r2*dvz - rv*dz;
        
        /* normal term */
        double nx = dy*dvz - dz*dvy;
        double ny = dz*dvx - dx*dvz;
        double nz = dx*dvy - dy*dvx;

        double prefact1 = A1*gr/r;
        double prefact2 = A2*gr/sqrt(tx*tx + ty*ty + tz*tz);
        double prefact3 = A3*gr/sqrt(nx*nx + ny*ny + nz*nz);

        // add acceleration to object
        com->ax += (prefact1*dx + prefact2*tx + prefact3*nx);
        com->ay += (prefact1*dy + prefact2*ty + prefact3*ny);
        com->az += (prefact1*dz + prefact2*tz + prefact3*nz);
    }
}

void relativity_correction_force(){
    struct particle* com = &(particles[13]);
    
    // relative to the Sun
    struct particle sun = particles[0];
    double dx = com->x - sun.x;
    double dy = com->y - sun.y;
    double dz = com->z - sun.z;
    
    double dvx = com->vx - sun.vx;
    double dvy = com->vy - sun.vy;
    double dvz = com->vz - sun.vz;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double _r = 1./r; // reciprocal
    
    // Benitez and Gallardo 2008
    // alpha = k^2/(r^3*c^2); beta  = 4k^2/r - v.v; gamma = 4*(v.r)
    // delta_a = alpha*(beta*r + gamma*v)
    double alpha = G*_r*_r*_r*_c2;
    double beta  = 4.*G*_r - (dvx*dvx + dvy*dvy + dvz*dvz);
    double gamma = 4.*(dvx*dx + dvy*dy + dvz*dz);
    
    // add acceleration to comet
    com->ax = alpha * (beta*dx + gamma*dvx);
    com->ay = alpha * (beta*dy + gamma*dvy);
    com->az = alpha * (beta*dz + gamma*dvz);
}

void additional_force(){
    if (N > N_active){ // "comet only" and if still alive, haha -_-
        outgassing_correction_force();
        relativity_correction_force();
    }
}

void problem_inloop(){ // simple direct collision-search only for the central body
	for (int i=N_active; i<N; i++){
        struct particle sun = particles[0];
        struct particle p = particles[i];
        double dx = p.x - sun.x;
        double dy = p.y - sun.y;
        double dz = p.z - sun.z;
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        // remove particles falling into the star (sun)
        if (r < 0.005){ 
            particles[i] = particles[N-1]; // -> rebound scheme to delete particle
            i--;
            N--;
        }
    }

    if (output_check(36.52422)){ // output heliocentric orbital elements every [input] 
		output_append_orbits("orbits.txt"); // output orbit from 1 to N (ignore the central object)
        // t, a, e, i, Omega, omega, l(mean longitude), P (period), f(true anomaly)
	}
}

void problem_output(){
	if (output_check(10.0*dt)){
		output_timing(); // output on terminal
	}
}

void problem_finish(){
}
