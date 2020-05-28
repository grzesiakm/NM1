#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../nrutil.h"
#include "../nrutil.c"


double f1(double x) {
	if(fabs(x) < 0.0000001) {
		x += 0.0000001;
    }

	return sin(x) / x;
}

double f2(double x) {
	if(fabs(x) < 0.0000001) {
		x += 0.0000001;
    }

	return (cos(x) - exp(x)) / sin(x);
}

double f3(double t) {
	if(fabs(t) < 0.0000001) {
		t += 0.0000001;
    }

	return exp(-1 / t) / t;
}

void Romberg(double** D, int n, double a, double b, double (*f)(double)) {
    D[0][0] = 0.5 * (b - a) * (f(a) + f(b));

    for (int w = 1; w <= n; w++) {
        double sum = 0.0;
        double hw = (b - a) / pow(2, w);

        for (int i = 1; i <= pow(2, w - 1); i++) {
            sum += f(a + (2.0 * i - 1) * hw);
        }

        D[w][0] = 0.5 * D[w -1][0] + hw * sum;
    }

    for (int i = 1; i <= n; i++) {
        for (int w = i; w <= n; w++) {
            D[w][i] = (pow(4, i) * D[w][i - i] - D[w - 1][i - 1]) / (pow(4, i) - 1);
        }
    }
}


int main() {
    FILE* file1 = fopen("out1.dat", "w");
    FILE* file2 = fopen("out2.dat", "w");
    FILE* file3 = fopen("out3.dat", "w");
    const int n[] = {7, 15, 7};
    const double a[] = {0.0, -1.0, 0.0};
    double (*f[])(double) = {f1, f2, f3};
    FILE* files[] = {file1, file2, file3};

    for(int i = 0; i < 3; i++) {
        double** D = dmatrix(0, n[i], 0, n[i]);
        Romberg(D, n[i], a[i], 1.0, f[i]);

        fprintf(files[i], "%s %15s %15s\n", "w", "D[w][0]", "D[w][w]");
        for (int j = 0; j <= n[i]; j++) {
            fprintf(files[i], "%d %15f %15f\n", j, D[j][0], D[j][j]);
        }

        fclose(files[i]);
        free_dmatrix(D, 0, n[i], 0, n[i]);
    }

    return 0;
}