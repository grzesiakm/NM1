#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#define N 5
#define IT_MAX 30

double licz_r(double * a, double * b, int n, double xj) {
	b[n] = 0;
	for(int i = n-1; i >= 0; i--)
		b[i] = a[i+1] + xj*b[i+1];
	return a[0] + xj*b[0];
}

int main() {
	double x0, x1, x2, R0, R1, R2;
	double* a = malloc((N+1)*sizeof(double));
	double* b = malloc((N+1)*sizeof(double));
	FILE* f = fopen("out.dat", "w");
	fprintf(f, "%s | %s | %s | %s |\n", "L", "it", "x_it+1", "R_it+1");
	a[0] = 240.0;
	a[1] = -196.0;
	a[2] = -92.0;
	a[3] = 33.0;
	a[4] = 14.0;
	a[5] = 1.0;
	for(int L = 1; L <= N; L++) {
		int n = N-L+1;
		x0 = 0.0;
		x1 = 0.1;
		R0 = licz_r(a, b, n, x0);
		R1 = licz_r(a, b, n, x1);
		for(int IT = 1;  IT <= IT_MAX; IT++) {
			x2 = x1-R1*(x1-x0)/(R1-R0);
			R2 = licz_r(a, b, n, x2);
			R0 = R1;
			R1 = R2;
			x0 = x1;
			x1 = x2;
			fprintf(f, "%d | %d | %g | %g |\n", L, IT, x2, R2);
			if(fabs(x1-x0) < 10e-8) break;
		}
		for(int i = 0; i <= n; i++)
			a[i] = b[i];
	}
	fclose(f);	
	free(a);
	free(b);
	return 0;
}
