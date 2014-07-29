#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void parse_cartesian_data(FILE *fi,int npl, double masses[npl], double pos[npl][3], double vel[npl][3]){
	double mass;
	double x,y,z;
	double vx,vy,vz;
	
	int i =0;
	while( fscanf(fi,"%lf %lf %lf %lf %lf %lf %lf",&mass,&x,&y,&z,&vx,&vy,&vz)==7){
		masses[i] = mass;
		pos[i][0] = x;	pos[i][1] = y;	pos[i][2] = z;
		vel[i][0] = vx;	vel[i][1] = vy;	vel[i][2] = vz;
		i++;
	}
}

int main()
{
	FILE * fi;
	double pos[5][3],vel[5][3], masses[5];

	fi = fopen("test.dat","r");
	
	parse_cartesian_data(fi,5,masses,pos,vel);
	int i;
	for (i=0;i<5;i++){
		printf("%d , %lf\n",i,masses[i]);
		printf("%lf %lf %lf\n",pos[i][0],pos[i][1],pos[i][2]);
		printf("%lf %lf %lf\n",vel[i][0],vel[i][1],vel[i][2]);
	}
	return(0);
}
