#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "gaussQuad.h"


#define PI 3.14159265358979323846

double L_x = 1.0;
double L_y = 1.0;

double f(double x, double y){
	return -2*(x-x*x+y-y*y);
}

double addContribution(size_t m, size_t n, double x, double y){
	// integrate with Gaussian quadratures
	
	unsigned int M = 64, N = 64;
	unsigned int MN = M*N;

	double points_xit[MN], points_etat[MN], weights[MN];
	gaussQuad2D(points_xit, points_etat, weights, M, N);

	double xit, etat, wt, xt, yt, a_mn = 0, J_2 = 0.25*L_x*L_y;

	double fact_x = (m+1)*PI/L_x;
	double fact_y = (n+1)*PI/L_y;
	for(unsigned int gp = 0; gp < MN; gp++){
		xit = points_xit[gp];
		etat = points_etat[gp];
		wt = weights[gp];

		xt = (xit+1)*L_x/2;
		yt = (etat+1)*L_y/2;

		a_mn +=  f(xt,yt)*sin(fact_x*xt)*sin(fact_y*yt)*J_2*wt;
		printf("a_mn = %e\n", a_mn);

	}
	printf("a_mn = %e\n", a_mn);
	a_mn *= -4/L_x/L_y/(fact_x*fact_x+fact_y*fact_y);
	printf("a_mn = %e\n", a_mn);
	
	exit(EXIT_FAILURE);	
	return a_mn*sin(fact_x*x)*sin(fact_y*y);
}
	

int main(int argc, char **argv){

		

	double uTerm, u = 0;

	double maxTerm, x = 0.2, y = 0.2;
	size_t N = 1, n, m;
	do{
		maxTerm = 0;
		n = N;
		for(size_t m = 0; m <= N; m++){
			uTerm = addContribution(m, n, x, y);
			if(fabs(uTerm) > maxTerm)
				maxTerm = fabs(uTerm);
			u += uTerm;

//			printf("uTerm = %.15f\n", uTerm);
//			exit(EXIT_FAILURE);
			printf("uTerm = %e\n", uTerm);
		}
		m = N;
		for(size_t n = 0; n < N; n++){
			uTerm = addContribution(m, n, x, y);
			if(fabs(uTerm) < maxTerm)
				maxTerm = fabs(uTerm);
			u += uTerm;
		}
		N++;

	}while(maxTerm > 1e-4);
	
	printf("u = %e\n\n", u);

	return 0;
}


