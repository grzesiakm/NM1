#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../gaussj.c"
#include "../nrutil.h"
#include "../nrutil.c"
#define _USE_MATH_DEFINES
#define n 100
#define xmax 2*M_PI
#define xmin 0.0
#define step 0.01


float fun1(float x) {
	return 2.0 * sin(x) + sin(2.0 * x) + 2.0 * sin(3.0 * x);
}

float fun14(float x) {
    float alpha = rand() / (RAND_MAX + 1.0) - 0.5;
	return 2.0 * sin(x) + sin(2.0 * x) + 2.0 * sin(3.0 * x) + alpha;
}

float fun2(float x) {
    return 2.0 * sin(x) + sin(2.0 * x) + 2.0 * cos(x) + cos(2.0 * x);
}

float fun3(float x) {
    return 2.0 * sin(1.1 * x) + sin(2.1 * x) + 2.0 * sin(3.1 * x);
}

float approx(float* ak, float* bj, int Ms, int Mc, float x) {
    float sum = 0.0;
    for(int k = 1; k <= Ms; k++) 
        sum += ak[k] * sin(k * x); 
    for (int j = 1; j <= Mc; j++)  
        sum += bj[j] * cos(j * x); 
    return sum + bj[0]/2.0;
}

int main() {

    srand(time(NULL));
    FILE* f1 = fopen("pkt.dat", "w");
    FILE* f2 = fopen("approx1.dat", "w");
    FILE* f3 = fopen("approx2.dat", "w");
    FILE* f4 = fopen("approx3.dat", "w");
    FILE* f5 = fopen("approx4.dat", "w");
    FILE* f6 = fopen("ab4.dat", "w");

    const float h = (xmax - xmin) / (n - 1);
    float* xi = vector(0, n - 1);

    for(int i = 0; i <= n - 1; i++)
	xi[i] = i * h;

    //fun1
    float* y1 = vector(0, n - 1);
    for(int i = 0; i <= n - 1; i++)
        y1[i] = fun1(xi[i]);

    for(int i = 0; i < n; i++)
        fprintf(f1, "%f %f\n", xi[i], y1[i]);
        
    fprintf(f1, "\n\n"); 

    unsigned int Ms = 5;
    unsigned int Mc = 5;
    float* a1 = calloc(sizeof(float), Ms + 1);
    float* b1 = calloc(sizeof(float), Mc + 1);

    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a1[k] += 2.0/((double)n) * y1[i] * sin(k * xi[i]);
            b1[k] += 2.0/((double)n) * y1[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d : a1[i] = %f, b1[i] = %f\n", i, a1[i], b1[i]);

    printf("\n\n");
    
    for(float x = xmin; x <= xmax; x += step)
        fprintf(f2, "%f %f\n", x, approx(a1, b1, Ms, Mc, x));   

    free_vector(y1, 0, n - 1);
    free(a1);
    free(b1);

    //fun2
    float* y2 = vector(0, n - 1);
    for(int i = 0; i <= n - 1; i++)
        y2[i] = fun2(xi[i]);

    for(int i = 0; i < n; i++)
        fprintf(f1, "%f %f\n", xi[i], y2[i]);
        
    fprintf(f1, "\n\n"); 

    float* a2 = calloc(sizeof(float), Ms + 1);
    float* b2 = calloc(sizeof(float), Mc + 1);

    
    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a2[k] += 2.0/((double)n) * y2[i] * sin(k * xi[i]);
            b2[k] += 2.0/((double)n) * y2[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d: a2[i] = %f, b2[i] = %f\n", i, a2[i], b2[i]);

    printf("\n\n");

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f3, "%f %f\n", x, approx(a2, b2, Ms, Mc, x)); 

    free_vector(y2, 0, n - 1);
    free(a2);
    free(b2);

    //fun3
    float* y3 = vector(0, n - 1);
    for(int i = 0; i <= n - 1; i++)
        y3[i] = fun3(xi[i]);

    for(int i = 0; i < n; i++)
        fprintf(f1, "%f %f\n", xi[i], y3[i]);
        
    fprintf(f1, "\n\n"); 

    Ms = 5;
    Mc = 0;
    float* a3 = calloc(sizeof(float), Ms + 1);
    float* b3 = calloc(sizeof(float), Mc + 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= Mc; j++)
            b3[j] += 2.0 / ((double)n) * y3[i] * cos(j * xi[i]);
        for (int k = 0; k <= Ms; k++) 
            a3[k] += 2.0 / ((double)n) * y3[i] * sin(k * xi[i]);
    }

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d : a3[i] = %f, b3[i] = %f\n", i, a3[i], b3[i]);

    printf("\n\n");

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f4, "%f %f\n", x, approx(a3, b3, Ms, Mc, x)); 

    fprintf(f4, "\n\n"); 

    free(a3);
    free(b3);

    Ms = Mc =5;
    float* a33 = calloc(sizeof(float), Ms + 1);
    float* b33 = calloc(sizeof(float), Mc + 1);

    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a33[k] += 2.0 / ((double)n) * y3[i] * sin(k * xi[i]);
            b33[k] += 2.0 / ((double)n) * y3[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i < Ms + 1; i++) 
        printf("i = %d : a33[i] = %f, b33[i] = %f\n", i, a33[i], b33[i]);

    printf("\n\n");

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f4, "%f %f\n", x, approx(a33, b33, Ms, Mc, x)); 

    fprintf(f4, "\n\n"); 

    free(a33);
    free(b33);

    Ms = Mc = 10;
    float* a333 = calloc(sizeof(float), Ms + 1);
    float* b333 = calloc(sizeof(float), Mc + 1);

    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a333[k] += 2.0 / ((double)n) * y3[i] * sin(k * xi[i]);
            b333[k] += 2.0 / ((double)n) * y3[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d : a333[i] = %f, b333[i] = %f\n", i, a333[i], b333[i]);

    printf("\n\n");

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f4, "%f %f\n", x, approx(a333, b333, Ms, Mc, x)); 

    free_vector(y3, 0, n - 1);
    free(a333);
    free(b333);

    //fun4
    float* y4 = vector(0, n - 1);
    for(int i = 0; i <= n - 1; i++)
        y4[i] = fun14(xi[i]);

    for(int i = 0; i < n; i++)
        fprintf(f1, "%f %f\n", xi[i], y4[i]);

    Ms = Mc = 5;
    float* a4 = calloc(sizeof(float), Ms + 1);
    float* b4 = calloc(sizeof(float), Mc + 1);

    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a4[k] += 2.0 / ((double)n) * y4[i] * sin(k * xi[i]);
            b4[k] += 2.0 / ((double)n) * y4[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i <= Ms; i++){
        fprintf(f6, "%d %f %f\n", i, a4[i], b4[i]);
    }

    fprintf(f6, "\n\n"); 

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d : a4[i] = %f, b4[i] = %f\n", i, a4[i], b4[i]);

    printf("\n\n");

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f5, "%f %f\n", x, approx(a4, b4, Ms, Mc, x));

    fprintf(f5, "\n\n"); 
    
    free(a4);
    free(b4);

    Ms = Mc = 30;
    float* a44 = calloc(sizeof(float), Ms + 1);
    float* b44 = calloc(sizeof(float), Mc + 1);

    for(int i = 0; i <= n - 1; i++) {
        for(int k = 0; k <= Ms; k++) {
            a44[k] += 2.0 / ((double)n) * y4[i] * sin(k * xi[i]);
            b44[k] += 2.0 / ((double)n) * y4[i] * cos(k * xi[i]);
        }
    }

    for(int i = 0; i <= Ms; i++){
        fprintf(f6, "%d %f %f\n", i, a44[i], b44[i]);
    }

    for(int i = 0; i <= Ms; i++) 
        printf("i = %d : a44[i] = %f, b44[i] = %f\n", i, a44[i], b44[i]);

    for(float x = xmin; x <= xmax; x += step)
        fprintf(f5, "%f %f\n", x, approx(a44, b44, Ms, Mc, x));

    free_vector(y4, 0, n - 1);
    free(a44);
    free(b44);
    free_vector(xi, 0, n - 1);

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);

    return 0;
}