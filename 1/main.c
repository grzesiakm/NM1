#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

#define N 400	 // rozmiar macierzy M: NxN
#define STEP 0.1

int main(void)
{
	float **M, **b;
	// Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);

	// Wypelnienie macierzy M i wektora b
	for (int i = 1; i <= N; ++i) {
		b[i][1] = 0.0;
		for (int j = 1; j <= N; ++j)
			M[i][j] = 0.0;
	}
	b[1][1] = 1.0;

	M[1][1] = 1.0;
	M[2][1] = -1.0;
	M[2][2] = 1.0;

	for (int i = 3; i <= N; i++) {
		M[i][i-2] = 1.0;
		M[i][i-1] = STEP*STEP - 2;
		M[i][i] = 1.0;
	}

	// Rozwiazanie ukladu rownan Mx=b - wywolanie procedury
	gaussj(M, N, b, 1);
	
	// Zapis wynikow do pliku
	FILE * fp;
	fp = fopen("out.dat", "w");

	// Wypisanie rozwiazania zapisanego w wektorze b
	for (int i = 1; i <= N; ++i)
		printf("%f %g\n", (i-1)*STEP, b[i][1]);

	fclose(fp);

	// Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}
