/**
 * @file 	boundaries.c
 * @brief 	Implementation of open boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * This file implements open boundary conditions. Every particle that leaves 
 * the box is removed from the simulation 
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
#include "integrator.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

int nghostx = 0;
int nghosty = 0;
int nghostz = 0;

void boundaries_check(){
}

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx  = 0;
	gb.shifty  = 0; 
	gb.shiftz  = 0;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}

/**
 * Checks if a given particle is within the computational domain.
 * @param p Particle to be checked.
 * @return Return value is 1 if particle is inside the box and 0 otherwise.
 */
int boundaries_particle_is_in_box(struct particle p){
	return 1;
}

