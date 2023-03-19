/*****************************************************
 *                                                   *
 * Candidate number: 24563                           *
 *                                                   *
 * Downhill Simplex algorithm for finding the        *
 * minimum of the Rosenbrock function                *
 *                                                   *
 * Starting vertices:                                *
 * P_0 = (0,0), P_1 = (2,0), P_2 = (0,2)             *
 *                                                   *
 * This code can by amended for higher dimensions by *
 * changing N to the number of variables in the new  *
 * function and changing the initial vertices        *
 *                                                   *
 *****************************************************/

/*#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define N 2 // Number of dimensions

// Constants for reflection, contraction, expansion, and shrinking
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
#define RHO 0.5

#define TOL 1e-8 // Tolerance for convergence
#define MAX_ITER 1000 // Maximum number of iterations

// Structure containing the coordinates of a point
typedef struct {
    double x[N];
} Point;

// Rosenbrock function
double func(Point P) {
    return 100 * pow(P.x[1] - pow(P.x[0], 2), 2) + pow(1 - P.x[0], 2);
}

// Sort points in ascending order based on function values
void sortPoints(double Y[], Point P[]) {
    double temp;
    Point pTemp;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N + 1; j++) {
            if (Y[i] > Y[j]) {
                temp = Y[i];
                Y[i] = Y[j];
                Y[j] = temp;

                pTemp = P[i];
                P[i] = P[j];
                P[j] = pTemp;
            }
        }
    }
}

// Calculate centroid (Pbar) for the first N points
void centroid(Point P[], Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pbar->x[i] = 0;
        for (int j = 0; j < N; j++) {
            Pbar->x[i] += P[j].x[i] / N;
        }
    }
}


// Reflect the worst (highest value) point (Ph) with respect to the centroid (Pbar)
void reflect(Point* Ps, Point* Pbar, Point P[]) {
    for (int i = 0; i < N; i++) {
        Ps->x[i] = Pbar->x[i] + ALPHA * (Pbar->x[i] - P[N].x[i]);
    }
}

// Contract the worst point (Ph) towards or reflected point (P*) away from the centroid (Pbar)
void contract(Point* Pss, Point P[], Point* Pbar, Point* Ps, bool inside) {
    for (int i = 0; i < N; i++) {
        if (inside) {
            Pss->x[i] = Pbar->x[i] + BETA * (P[N].x[i] - Pbar->x[i]);
        }
        else {
            Pss->x[i] = Pbar->x[i] + BETA * (Ps->x[i] - Pbar->x[i]);
        }
    }
}

// Expand the reflected point (P*) further away from the centroid (Pbar)
void expand(Point* Pss, Point* Ps, Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pss->x[i] = Pbar->x[i] + GAMMA * (Ps->x[i] - Pbar->x[i]);
    }
}

// Shrink all points towards the best (lowest value) point (Pl)
void shrink(Point P[]) {
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N; j++) {
            P[i].x[j] = P[0].x[j] + RHO * (P[i].x[j] - P[0].x[j]);
        }
    }
}

// Replace an original point (Ph) with a new point (P* or P**)
void replacePoint(Point* new, Point* orig) {
    *new = *orig;
}

// Test if the simplex has reached convergence (standard deviation < 10^-8)
bool minCon(double Y[]) {
    double Ybar = 0, sum = 0;
    for (int i = 0; i < N + 1; i++) {
        Ybar += Y[i] / (N + 1);
    }
    for (int i = 0; i < N + 1; i++) {
        sum += pow(Y[i] - Ybar, 2) / N;
    }
    return (sqrt(sum) < TOL);
}

// Downhill simplex algorithm implementation
void simplex(Point P[]) {
    Point Pbar, Ps, Pss;
    double Ys, Yss, Y[N + 1];

    int a;

    for (a = 0; a < MAX_ITER; a++) {

        // Evaluate function values for each point
        for (int i = 0; i < N + 1; i++) {
            Y[i] = func(P[i]);
        }

        // Sort points and find centroid
        sortPoints(Y, P);
        centroid(P, &Pbar);

        // Reflect Ph
        reflect(&Ps, &Pbar, P);
        Ys = func(Ps);

        // If Y* is greater than or equal to Yl but less than the second worst value
        if (Ys >= Y[0] && Ys < Y[N - 1]) {
            replacePoint(&P[N], &Ps); // Replace Ph with P*
        }

        // If Y* is less than Yl
        else if (Ys < Y[0]) {
            expand(&Pss, &Ps, &Pbar); // Calculate the expanded point P**
            Yss = func(Pss);

            // If Y** is less than Y*
            if (Yss < Ys) {
                replacePoint(&P[N], &Pss); // Replace Ph with P**
            }

            else {
                replacePoint(&P[N], &Ps); // Replace Ph with P*
            }
        }

        else {
            
            // If Y* is less than Yh
            if (Ys < Y[N]) {
                contract(&Pss, P, &Pbar, &Ps, 0); // Calculate P**; contract P* away from Pbar
                Yss = func(Pss);

                // If Y** is less than Y*
                if (Yss < Ys) {
                    replacePoint(&P[N], &Pss); // Replace Ph with P**
                }

                else {
                    shrink(P); // Shrink all points towards Pl
                }
            }

            else {
                contract(&Pss, P, &Pbar, NULL, 1); // Calculate P**; contract Ph towards Pbar
                Yss = func(Pss);

                // If Y** is less than Yh
                if (Yss < Y[N]) {
                    replacePoint(&P[N], &Pss); // Replace Ph with P**
                }

                else {
                    shrink(P); // Shrink all points towards Pl
                }
            }
        }

        // Check for convergence
        if (minCon(Y)) {
            break;
        }
    }

    // Print results
    printf("Final vertices:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("(");
        for (int j = 0; j < N; j++) {
            printf("%lf", P[i].x[j]);
            if (j < N - 1) {
                printf(", ");
            }
        }
        printf(")\n");
    }

    printf("\nEvaluations:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("%lf\n", func(P[i]));
    }

    printf("\nIterations: %d\n", a + 1);
}

int main() {
    // Define initial points
    Point P[N + 1] = { {0, 0}, {2, 0}, {0, 2} };

    // Run the Downhill Simplex method
    simplex(P);

    return 0;
}*/