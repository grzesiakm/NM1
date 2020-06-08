#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../nrutil.h"
#include "../nrutil.c"
#include "../gauleg.c"
#include "../gaulag.c"
#include "../gammln.c"
#include "../gauher.c"

#define A 0.0
#define B 2.0
#define c3_dok 0.1919832644


float c1(float x) {
    return x / (4 * x * x + 1);
}

float c2(float x, int k) {
    return pow(x, k);
}

float c3(float x, float y) {
    return pow(sin(x), 2) * pow(sin(y), 4);
}

float c1_a(float a, float b) {
    return 1.0 / 8.0 * (log(fabs(4.0 * b * b + 1.0)) - log(fabs(4.0 * a * a + 1.0)));
}

float c2_a(int k) {
    int result = 1;
    for (int i = 1; i <= k; i++) {
        result *= i;
    } 
    return result;
}

void Gauss_Legandre1(int n1, int n2, float(*fun)(float), float(*c)(float, float), FILE* f) {
    for (int n = n1; n <= n2; n++) {
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauleg(A, B, x, w, n);
        float sum = 0.0;
        for (int i = 1; i <= n; i++) {
            sum += fun(x[i]) * w[i];
        }
        fprintf(f, "%d %f\n", n, fabs(sum - c(A, B)));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }
}

void Gauss_Legandre2(int n1, int n2, int k1, int k2, float(*fun)(float, int), float(*c)(int), FILE* f) {
    for (int k = k1; k <= k2; k += 5) {
        for (int n = n1; n <= n2; n++) {
            float* x = vector(1, n);
            float* w = vector(1, n);
            
            gaulag(x, w, n, 0);
            float sum = 0.0;
            for (int i = 1; i <= n; i++) {
                sum += fun(x[i], k) * w[i];
            }
            fprintf(f, "%d %f\n", n, fabs(sum - c(k)));
        
            free_vector(x, 1, n);
            free_vector(w, 1, n);
        }
        fprintf(f, "\n\n");
    }
}

void Gauss_Hermite(int n1, int n2, float(*fun)(float, float), FILE* f) {
    for (int n = n1; n <= n2; n++) {
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauher(x, w, n);
        float sum = 0.0;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                sum += fun(x[i], x[j]) * w[i] * w[j];
            }
        }
        fprintf(f, "%d %f\n", n, fabs(sum - c3_dok));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }
}


int main(){
    FILE* f1 = fopen("c1.dat", "w");
    FILE* f2 = fopen("c2.dat", "w");
    FILE* f3 = fopen("c3.dat", "w");

    Gauss_Legandre1(2, 20, c1, c1_a, f1);
    Gauss_Legandre2(2, 20, 5, 10, c2, c2_a, f2);
    Gauss_Hermite(2, 15, c3, f3);

    fclose(f1);
    fclose(f2);
    fclose(f3);
    return 0;
}