#include <stdio.h>
#include <stdlib.h>
#define krok 0.01


void liczIloraz(double* fm, double* xm, double* ym, int rzadI) {
	for(int j = 0; j < rzadI; j++) {
		for(int k = 0; k <= j; k++) {
			double el = 1;
			for(int l = 0; l <= j; l++) {
				if (l != k)
					el *= 1.0 / (xm[k] - xm[l]);
			}
			fm[j] += el * ym[k];
		}
	}
}

double f(double x) {
	return 1.0 / (1.0 + x*x);
}

double interpolujNewton(double x, int n, double* xm, double* fm) {
	double wynikN = 0;
	for(int j = 0; j <= n; j++) {
		double el = 1;
		for(int i = 0; i < j; i++) {
			el *= x - xm[i];
		}
		wynikN += el * fm[j];
	}
	return wynikN;
}

int main() {
	int n = 6; //stopien wielomianu interpolacyjnego
	int J = 7; //rzad ilorazu

	double* fm = malloc((n+1) * sizeof(double)); //tablica wartosci ilorazow roznicowych
	double xm[] = {-5.0, -2.0, -0.5, 0.0, 0.5, 2.0, 5.0}; //tablica polozen wezlow
	double* ym = malloc((n+1) * sizeof(double)); //tablica wartosci funkcji f w wezlach
	double* fmRowne = malloc((n+1) * sizeof(double));
	double* xmRowne = malloc((n+1) * sizeof(double));
	double* ymRowne = malloc((n+1) * sizeof(double));

	FILE* f1 = fopen("zad_1.dat", "w");
	FILE* f2 = fopen("zad_2.dat", "w");
	FILE* f3 = fopen("ilorazy.dat", "w");

	for(int h = 0; h < J; h++) {
		fm[h] = 0.0;
		fmRowne[h] = 0.0;

		ym[h] = f(xm[h]);

		xmRowne[h] = -5.0 + h * (10.0 / (n));
		ymRowne[h] = f(xmRowne[h]);
	}

	liczIloraz(fm, xm, ym, J);
	liczIloraz(fmRowne, xmRowne, ymRowne, J);

	for(double x = -5.0; x <= 5.0; x += krok) {
		fprintf(f1, "%.2f %f\n", x, interpolujNewton(x, n, xm, fm));
		fprintf(f2, "%.2f %f\n", x, interpolujNewton(x, n, xmRowne, fmRowne));
	}

	fprintf(f3, "%s | %39s | %39s \n", "j", "ilorazy roznicowe dla nierownych wezlow", "ilorazy roznicowe dla rownych wezlow");	
	for(int j = 0; j < J; j++)
		fprintf(f3, "%d | %39g | %39g \n", j, fm[j], fmRowne[j]);

	fclose(f1);
	fclose(f2);
	fclose(f3);

	free(fm);
	free(ym);
	free(fmRowne);
	free(xmRowne);
	free(ymRowne);

	return 0;	
}
