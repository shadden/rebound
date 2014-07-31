int main()
{
	FILE * fi;
	int npl, Ndt,NdtT,i,j;
	double dt;
	
	fi = fopen("test.par","r");
	parse_param_data(fi,&npl,&dt,&Ndt,&NdtT);
	fclose(fi);

	double  **pos,**vel,*masses;
	pos = (double **)malloc(sizeof(double *)*npl);
	vel = (double **)malloc(sizeof(double *)*npl);

	for(i=0; i < npl; i++) {
	  pos[i] = (double *)malloc(sizeof(double)*3);
	  vel[i] = (double *)malloc(sizeof(double)*3);
	}
	masses = (double *)malloc(sizeof(double)*npl);
	
	fi = fopen("test.dat","r");
	parse_cartesian_data(fi,npl,masses,pos,vel);
	fclose(fi);

	for (i=0;i<npl;i++){
		printf("%d , %lf\n",i,masses[i]);
		printf("%lf %lf %lf\n",pos[i][0],pos[i][1],pos[i][2]);
		printf("%lf %lf %lf\n",vel[i][0],vel[i][1],vel[i][2]);
	}
	for(i=0; i < Npl; i++) {
	  free(pos[i]);
	  free(vel[i]);
	}
	free(pos);
	free(vel);
	free(masses);

	return(0);
}
