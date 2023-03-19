/*#include <stdio.h>
#include <math.h>

#define N 2 // Number of dimensions

// Structure containing the coordinates of a point
typedef struct {
    double x[N];
} Point;

// Rosenbrock function
double func(Point P) {
    return 100 * pow(P.x[1] - pow(P.x[0], 2), 2) + pow(1 - P.x[0], 2);
}

//Save function values for a range of input values to a file
void file_data(char* textfile, double xlow, double xhigh, int values) {
    double x0 = xlow;
    double step = (xhigh - xlow) / (values - 1);

    FILE* fp = fopen(textfile, "w");

    for (int i = 0; i < 100; i++, x0 += step) {
        Point P = { x0, 1 };
        fprintf(fp, "%.3lf, %.3lf\n", x0, func(P));
    }
    fclose(fp);
}

int main() {
    // Generate and display values in text file
    file_data("data.txt", -2., 2., 100);

    return 0;
}*/