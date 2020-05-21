#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../four1.c"
#include "../nrutil.h"
#include "../nrutil.c"
#define T 1.0
#define MAX(X,Y) ((X)>(Y) ? (X):(Y))


float randDelta() {
    return (float)rand()/RAND_MAX - 0.5;
}

int main() {
    srand(time(NULL));
    FILE *f1 = fopen("k8.dat", "w");
    FILE *f2 = fopen("k10.dat", "w");
    FILE *f3 = fopen("k12.dat", "w");
    FILE* files[] = {f1, f2, f3};
    const int k[] = {8, 10, 12};

    for(int i = 0; i < 3; i++){
        FILE* current = files[i];
        const int N = pow(2, k[i]);
        const float dt = (float)(3 * T / N);
        const float sigma = T / 20.0;
        float* f = vector(1, 2 * N);
        float* y = vector(1, 2 * N);
        float* g1 = vector(1, 2 * N);
        float* g2 = vector(1, 2 * N);

        for(int i = 1; i <= N; i++) {
            float ti = dt * (i - 1);
            f[2 * i - 1] = y[2 * i - 1] = sin((2 * M_PI / T) * ti) +  sin(2 * (2 * M_PI / T) * ti) + sin(3 * (2 * M_PI / T) * ti) + randDelta();
			g1[2 * i - 1] = g2[2 * i - 1] = (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-(ti * ti) / (2 * sigma * sigma));
			f[2 * i] = y[2 * i] = g1[2 * i] = g2[2 * i] = 0;

		}

		four1(f, N, 1);
		four1(g1, N, 1);
		four1(g2, N, -1);

        for (int i = 1; i <= N; i++) {
			float a1 = f[2 * i - 1];
			float b1 = f[2 * i];
			float a2 = g1[2 * i - 1] + g2[2 * i - 1];
			float b2 = g1[2 * i] + g2[2 * i];
			f[2 * i - 1] = a1 * a2 - b1 * b2;
			f[2 * i] = a1 * b2 + a2 * b1;
		}

		four1(f, N, -1);
        float max = fabs(f[1]);

        for(int i = 2; i <= N; i++) {
            max = MAX(fabs(max), fabs(f[2 * i - 1]));
        }

        for(int j = 1; j <= N; j++) {
            float ti = dt * (j - 1);
			fprintf(current, "%10f %10f\n", ti, y[2 * j - 1]);
        }

        fprintf(current, "\n\n");

        for (int j = 1; j <= N; j++) {
			float ti = dt * (j - 1);
            fprintf(current, "%10f %10f\n", ti, f[2 * j - 1] * 2.5 / max);
        }

        free_vector(f, 1, 2 * N);
		free_vector(y, 1, 2 * N);
		free_vector(g1, 1, 2 * N);
		free_vector(g2, 1, 2 * N);
        fclose(current);
    }

    return 0;
}