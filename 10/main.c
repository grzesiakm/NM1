#include <stdio.h>
#include <math.h>
#define h 0.1
#define max 1000


double f(double x, double y) {
    return 5.0/2.0 * (x * x - y) * (x * x - y) + (1 - x) * (1 - x);
}

double dfdx(double xi, double yi, double dx) {
    return (f(xi + dx, yi) - f(xi - dx, yi)) / (2.0 * dx);
}

double dfdy(double xi, double yi, double dy) {
    return (f(xi, yi + dy) - f(xi, yi - dy)) / (2.0 * dy);
}


int main() {
    FILE* f1 = fopen("eps1.dat", "w");
    FILE* f2 = fopen("eps2.dat", "w");

	const double step = 1.0e-4;
	const double e1 = 1.0e-2;
	const double e2 = 1.0e-3;
    double x0 = -0.75;
    double y0 = 1.75;

	fprintf(f1, "%f %f\n", x0, y0);
	fprintf(f2, "%f %f\n", x0, y0);

	//epsilon1
    for (int i = 0; i < max; i++) {
    	double x1 = x0 - h * dfdx(x0, y0, step);
		double y1 = y0 - h * dfdy(x0, y0, step);
        
		fprintf(f1, "%f %f\n", x1, y1);

	    if (sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)) < e1) {
	        break;
	    }
	
		x0 = x1;
		y0 = y1;
    }

	//epsilon2
    x0 = -0.75;
    y0 = 1.75;

    for (int i = 0; i < max; i++) {
        double x1 = x0 - h * dfdx(x0, y0, step);
        double y1 = y0 - h * dfdy(x0, y0, step);

        fprintf(f2, "%f %f\n", x1, y1);

	    if (sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)) < e2) {
            break;
        }

        x0 = x1;
        y0 = y1;
    }

    fclose(f1);
    fclose(f2);

    return 0;
}
