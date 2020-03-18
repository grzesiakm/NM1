#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/ludcmp.c"
#include "/opt/NR/numerical_recipes.c/lubksb.c"
#define N 3

//Wypisywanie macierzy
void printMatrix(char * txt, float ** matrix, int size) {
	printf("%s", txt);
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j ++) {
			printf("%12g ", matrix[i][j]);
		}
		printf("\n");
	}
}
//Zapisywanie macierzy do pliku
void fprintMatrix(FILE * file, char * txt, float ** matrix, int size) {
	fprintf(file, "%s", txt);
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j ++) {
			fprintf(file, "%12g ", matrix[i][j]);
		}
		fprintf(file, "\n");
	}
}

int main(void) {

	FILE * file = fopen("out.dat", "w");
	float ** A  = matrix(1, N, 1, N); //Alokacja macierzy n x n
	float ** B = matrix(1, N, 1, N);
	int iterr = 1;
	//Uzupelnianie macierzy
	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			A[i][j] = iterr;
			B[i][j] = iterr;
			iterr++;
		}
	}
	B[1][1] = 1.1;
	//Wypisywanie macierzy
	printMatrix("Macierz A:\n", A, N);
	fprintMatrix(file, "Macierz A:\n", A, N);
	printMatrix("\nMacierz B:\n", B, N);
	fprintMatrix(file, "\nMacierz B:\n", B, N);
	//Alokacja wektorow permutacji oraz wypelnienie ich procedura ludcmp i inicjalizacja znaku permutacji
	int * indxA = ivector(1, N);
	int * indxB = ivector(1, N);
	float dA;
	float dB;
	ludcmp(A, N, indxA, &dA);
	ludcmp(B, N, indxB, &dB);
    	//Alokacja macierzy L i U dla macierzy A i B
	float ** LA  = matrix(1, N, 1, N);
	float ** LB = matrix(1, N, 1, N);
	float ** UA  = matrix(1, N, 1, N);
	float ** UB = matrix(1, N, 1, N);
	//Uzupelnienie macierzy L i U
	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			if (i == j) {
				LA[i][j] = 1;
				LB[i][j] = 1;
				UA[i][j] = A[i][j];
				UB[i][j] = B[i][j];
			}
			else if(i < j) {
				LA[i][j] = 0;
				LB[i][j] = 0;
				UA[i][j] = A[i][j];
				UB[i][j] = B[i][j];
			}
			else {
				LA[i][j] = A[i][j];
				LB[i][j] = B[i][j];
				UA[i][j] = 0;
				UB[i][j] = 0;
			}
			iterr++;
		}
	}
	//Wypisywanie macierzy L i U
	printMatrix("\nMacierz L dla A:\n", LA, N);
	fprintMatrix(file, "\nMacierz L dla A:\n", LA, N);
	printMatrix("\nMacierz L dla B:\n", LB, N);
	fprintMatrix(file, "\nMacierz L dla B:\n", LB, N);
	printMatrix("\nMacierz U dla A:\n", UA, N);
	fprintMatrix(file, "\nMacierz U dla A:\n", UA, N);
	printMatrix("\nMacierz U dla B:\n", UB, N);
	fprintMatrix(file, "\nMacierz U dla B:\n", UB, N);
	//Alokacja i uzupelnienie wektorow wyrazow wolnych, na ktorych wykonujemy procedure lubksb
	float * a1 = vector(1, N);
	float * a2 = vector(1, N);
	float * a3 = vector(1, N);
	float * b1 = vector(1, N);
	float * b2 = vector(1, N);
	float * b3 = vector(1, N);

	a1[1] = 1; a1[2] = 0; a1[3] = 0;
	a2[1] = 0; a2[2] = 1; a2[3] = 0;
	a3[1] = 0; a3[2] = 0; a3[3] = 1;
	b1[1] = 1; b1[2] = 0; b1[3] = 0;
	b2[1] = 0; b2[2] = 1; b2[3] = 0;
	b3[1] = 0; b3[2] = 0; b3[3] = 1;
	lubksb(A, N, indxA, a1);
	lubksb(A, N, indxA, a2);
	lubksb(A, N, indxA, a3);
	lubksb(B, N, indxB, b1);
	lubksb(B, N, indxB, b2);
	lubksb(B, N, indxB, b3);
	//Alokacja, uzupelnienie i wypisanie macierzy odwrotnych
	float ** A_inv = matrix(1, N, 1, N);
	float ** B_inv = matrix(1, N, 1, N);

	for(int i = 1; i <= N; i++) {
        	A_inv[i][1] = a1[i];
        	A_inv[i][2] = a2[i];
        	A_inv[i][3] = a3[i];
        	B_inv[i][1] = b1[i];
        	B_inv[i][2] = b2[i];
        	B_inv[i][3] = b3[i];
    	}
	printMatrix("\nMacierz A^-1:\n", A_inv, N);
	fprintMatrix(file, "\nMacierz A^-1:\n", A_inv, N);
	printMatrix("\nMacierz B^-1:\n", B_inv, N);
	fprintMatrix(file, "\nMacierz B^-1:\n", B_inv, N);
	iterr = 1;

	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
  			A[i][j] = iterr;
  			B[i][j] = iterr;
  			iterr++;
  		}
  	}
  	B[1][1] = 1.1;
	//Liczenie norm macierzy i macierzy odwrotnych
	float maxA = 0;
	float maxB = 0;
	float maxA_inv = 0;
	float maxB_inv = 0;

	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			if(fabs(A[i][j]) > maxA)
				maxA = fabs(A[i][j]);
			if(fabs(B[i][j]) > maxB)
				maxB = fabs(B[i][j]);
			if(fabs(A_inv[i][j]) > maxA_inv)
				maxA_inv= fabs(A_inv[i][j]);
			if(fabs(B_inv[i][j]) > maxB_inv)
				maxB_inv= fabs(B_inv[i][j]);
		}
	}

	float kappaA = maxA * maxA_inv;
	float kappaB = maxB * maxB_inv;

	printf("\n||A|| = %g\n||A^-1|| = %f = %g\nKappaA = %f = %g\n", maxA, maxA_inv, maxA_inv, kappaA, kappaA);
	fprintf(file, "\n||A|| = %g\n||A^-1|| = %f = %g\nKappaA = %f = %g\n", maxA, maxA_inv, maxA_inv, kappaA, kappaA);
	printf("\n||B|| = %g\n||B^-1|| = %f\nKappaB = %f\n", maxB, maxB_inv, kappaB);	
	fprintf(file, "\n||B|| = %g\n||B^-1|| = %f\nKappaB = %f\n", maxB, maxB_inv, kappaB);		
	iterr = 1;

	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			A[i][j] = iterr;
			B[i][j] = iterr;
			iterr++;
		}
	}
	B[1][1] = 1.1;
	//Alokowanie macierzy dla iloczynu macierzy oraz uzupelnianie i wypisywanie ich
	float ** A_product = matrix(1, N, 1, N);
	float ** B_product = matrix(1, N, 1, N);
	
	for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
			A_product[i][j] = B_product[i][j] = 0;
	int k;
	float sumA = 0, sumB = 0;

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			sumA = sumB = 0;
			for (int k = 1; k <= N; k++) {
				sumA = sumA + A[i][k] * A_inv[k][j];
				sumB = sumB + B[i][k] * B_inv[k][j];
			}
			A_product[i][j] = sumA;
			B_product[i][j] = sumB;
		}
	}
	printMatrix("\nA x A^-1:\n", A_product, N); 
	fprintMatrix(file, "\nA x A^-1:\n", A_product, N); 
	printMatrix("\nB x B^-1:\n", B_product, N);
	fprintMatrix(file, "\nB x B^-1:\n", B_product, N);

	fclose(file);
	//Zwalnianie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);
	free_matrix(LA, 1, N, 1, N);
	free_matrix(LB, 1, N, 1, N);
	free_matrix(UA, 1, N, 1, N);
	free_matrix(UB, 1, N, 1, N);
	free_matrix(A_inv, 1, N, 1, N);
	free_matrix(B_inv, 1, N, 1, N);
	free_matrix(A_product, 1, N, 1, N);
	free_matrix(B_product, 1, N, 1, N);
	free_vector(a1, 1, N);
	free_vector(a2, 1, N);
	free_vector(a3, 1, N);
	free_vector(b1, 1, N);
	free_vector(b2, 1, N);
	free_vector(b3, 1, N);
	free_ivector(indxA, 1, N);
	free_ivector(indxB, 1, N);
	
	return 0;
}
