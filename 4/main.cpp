#include <iostream>
#include <cmath>
#include</usr/include/gsl/gsl_eigen.h>

int delta(int i, int j) {
	if (i == j)
		return 1;
	return 0;
}

double ro(int alfa, double x) {
	return 1 + 4 * alfa * x * x;
}


int main(void) {
	FILE * file1 = fopen("result.dat", "w");
	FILE * file2 = fopen("vectors0.dat", "w");
	FILE * file3 = fopen("vectors100.dat", "w");

	int n = 200, N = 1;
	double L = 10, x;
	double delta_x = L / (n+1);

	gsl_matrix * A = gsl_matrix_calloc(n, n);
	gsl_matrix * B = gsl_matrix_calloc(n, n);
	gsl_vector *eval = gsl_vector_calloc(n);
	gsl_matrix *evec = gsl_matrix_calloc(n,n);
	gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);


	for (int alfa = 0; alfa <= 100; alfa +=2) {

		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				x = -L/2.0 + delta_x* (i + 1);
				gsl_matrix_set(A, i, j, (-1 * delta(i, j+1) + 2 * delta(i,j) - delta(i, j-1)) / (delta_x * delta_x));
				gsl_matrix_set(B, i, j, (ro(alfa, x) * delta(i, j)) / N);
			}
		}

		gsl_eigen_gensymmv(A, B, eval, evec, w);
		gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		fprintf(file1, "%3d ", alfa);
		for (int k = 0; k < 6; k++) 
			fprintf(file1, "%12g ", sqrt(gsl_vector_get(eval, k)));
		fprintf(file1, "\n");

		if (alfa == 0) {
			for (int i = 0; i < n; i++) {
				x = -L/2.0 + delta_x * (i + 1);
				fprintf(file2, "%12g ", x);
				for (int k = 0; k < 6; k++) 
					fprintf(file2, "%12g ", (gsl_matrix_get(evec,i,k)));
				fprintf(file2, "\n");
			}
		}

		if (alfa == 100) {
			for (int i = 0; i < n; i++) {
				x = -L/2.0 + delta_x * (i + 1);
				fprintf(file3, "%12g ", x);
				for (int k = 0; k < 6; k++) 
					fprintf(file3, "%12g ", (gsl_matrix_get(evec,i,k)));
				fprintf(file3, "\n");
			}
		}
	}

	fclose(file1);
	fclose(file2);
	fclose(file3);

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_gensymmv_free(w);


	return 0;
}
