#include <stdio.h>
#include <math.h>
#include <omp.h>

double rect_integral(double, double, double);

int main(void) {
    double x1 = 0;
    double x2 = 2 * M_PI;
    double dx = 0.00001;

    double integral = rect_integral(x1, x2, dx);

    printf("openmp will be using %d threads\n", omp_get_max_threads());

    printf(
        "numerical integral of sin(x) over %.5f to %.5f with the step of %.5f = %.5f\n",
        x1, x2, dx,
        integral
    );

    return EXIT_SUCCESS;
}

double rect_integral(double x1, double x2, double dx) {
    int steps = (x2 - x1) / dx;
    double area = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:area)
    for (i = 1; i <= steps; i++)
        area += sin(x1 + i * dx);

    area *= dx;

    return area;
}
