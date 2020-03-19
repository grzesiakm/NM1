#include <cmath>
#include <cstdio>
#include <iostream>
#include <chrono>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

#define N 1000
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))

void multiplyAvec (double** A, double* b, double* result, int n, int m){
	for (int i = 0; i < n; i++){
		result[i] = 0.0;
		for (int j = max(0, i - m); j <= min(n - 1, i + m); j++)
			result[i] += A[i][j] * b[j];
	}
}

double multiplyScalar(double* a, double* b, int n){
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return sum;
}

int main(void) {
	double** A = new double* [N];
	for(int i = 0; i < N; i++)
		A[i] = new double [N];
	double* b = new double [N];
	double* x = new double [N];
	int m = 5;

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if (abs(i - j) <=  m)
				A[i][j] =  1./(1 + (abs(i - j)));
			else
				A[i][j] =  0;
		}
		b[i] = i + 1;
		x[i] = 0;
	}

	double* r = new double[N];
	double* v = new double[N];
	double* y = new double[N];
	double* tmp = new double[N];

	double alfa, beta;

	auto timeStart = std::chrono::steady_clock::now();
	multiplyAvec(A, x, y, N, m);

	for (int i = 0; i < N; i++){
		v[i] = r[i] = b[i] - y[i];
	}

	FILE * file = fopen("result.dat", "w");
	int k = 0;
	while(sqrt(multiplyScalar(r, r, N)) > pow(10, -6)) {

		double tmpScalar = multiplyScalar(r, r, N);
		double tmpScalar2 = multiplyScalar(x, x, N);
		multiplyAvec(A, v, tmp, N, m);
		alfa = tmpScalar / multiplyScalar(v, tmp, N);

		for(int i = 0; i < N; i++) {
			x[i] += alfa * v[i];
			r[i] -= alfa * tmp[i];
		}

		beta = multiplyScalar(r, r, N) / tmpScalar;

		for(int i = 0; i < N; i++) {
			v[i] = v[i] * beta + r[i];
		}

		fprintf(file, "%d %g %g %g %g\n", k, sqrt(tmpScalar), alfa, beta, sqrt(tmpScalar2));
		k++;
	}

	auto timeStop = std::chrono::steady_clock::now();

	std::cout << "Part one execution time: " << std::chrono::duration<double, std::milli>(timeStop - timeStart).count() << "ms" << std::endl;

	fclose(file);


	float** A2 = matrix(1, N, 1, N);
	float** B2 = matrix(1, N, 1, 1);

	for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j++) {
			if (abs(i - j) <=  m)
				A2[i][j] =  1./(1 + (abs(i - j)));
			else
				A2[i][j] =  0;
		}
		B2[i][1] = i + 1;
	}

	timeStart = std::chrono::steady_clock::now();

	gaussj(A2,N,B2,1);

	timeStop = std::chrono::steady_clock::now();

	std::cout << "Part two execution time: " << std::chrono::duration<double, std::milli>(timeStop - timeStart).count() << "ms" << std::endl;

	for(int i = 0; i < N; i++)
		delete[] A[i];

	delete[] A;
	delete[] b;
	delete[] x;
	delete[] v;
	delete[] r;
	delete[] y;
	delete[] tmp;

	free_matrix(A2,1,N,1,N);
	free_matrix(B2,1,N,1,1);

	return 0;
}
