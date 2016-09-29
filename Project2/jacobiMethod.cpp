#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;

double maxOffDiagonals(arma::mat &A, int &k, int &l, int n) {

    double max=0;
    for (int i = 0; i < n; ++i) {

        // Indexing will only include upper triangular part
        // This assumes that matrix is symmetric
        // Advantage: excludes diagonal, which is what we want
        for (int j = i+1; j < n; ++j) {

            double aij = fabs(A(i,j));

            if (aij > max) {
                // Maximum non-diagonal value
                max = aij;
                // Position of this value
                k = i;
                l = j;
            }
        }

    }

    cout << "Largest off-diagonal element is: " << max << endl;

    return max;
}

void jacobiRotation(arma::mat &A, arma::mat &R, int &k, int &l, int n) {
    // A is the input matrix and S is the matrix with eigenvectors after
    // enough rotations

    /* -----------------------------------------------
     * COMPUTING c = cosine(angle) AND s = sine(angle)
     * -----------------------------------------------
     *
     * I don't want to compute cosine or sine of anything.
     * (Important steps to the right)
     * tan(angle) = t = sin / cos = s / c
     * cot(2*angle) = tau = (a_ll - a_kk) / 2*a_kl 	->	(No 1: tau)
     * using: cot(2*angle) = 1/2 * (cot(angle) - tan(angle))
     * -> t**2 + 2*tau*t - 1 = 0
     * which gives t = +/- 1/ (tau + sqrt(1 + tau**2))->(No. 2: t)
     * FINALLY: c = 1 / (sqrt(1 + t**2)) AND s = c*t -> (No.3: c and s)
     */

    // No. 1: tau
    double tau = (A(l,l) - A(k,k)) / (2*A(k,l));

    // No. 2: t
    // t can be either + or - depending on t**2 + 2*tau*t
    double t;
    if ( tau >= 0 ) {
        t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
        t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }

    // No. 3:
    double c = 1.0 / sqrt(1 + t*t);				// cos(angle)
    double s = c*t;								// sin(angle)
    if(A(k,l) == 0) {
        c = 1.0;
        s = 0.0;
    }

    /* ----------------------------------------------
     * ROTATING
     * ----------------------------------------------
     * for i neq k, i neq l:
     * 1) b_ii = a_ii
     * 2) b_ik = a_ik*c - a_il*s
     * 3) b_il = a_il*c - a_ik*s
     * for all i:
     * 4) b_kk = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s
     * 5) b_ll = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s
     * 6) b_kl = (a_kk - a_ll)*c*s + a_kl*(c*c - s*s)
     * The last term should be zero by definition of c and s
     * can be hardcoded to zero
     * b_?? are not used, instead a_?? is overwritten
     * this means that 1) is unnecesary
     */

    // Declaring variables
    double a_kk = A(k, k);
    double a_ll = A(l, l);
    A(k ,k) = a_kk*c*c - 2*A(k, l)*c*s + a_ll*s*s; // 4)
    A(l, l) = a_ll*c*c + 2*A(k, l)*c*s + a_kk*s*s; // 5)
    A(l, k) = 0.0; // Hard-coding non-diagonal elements,
    A(k, l) = 0.0; // these should be equal to zero

    for (int i = 0; i < n; i++) {
        // i neq k, i neq l
        if (i != k && i != l) {
            double a_ik = A(i, k);
            double a_il = A(i, l);
            A(i, k) = a_ik*c - a_il*s; // 2)
            A(k, i) = A(i, k); // Symmetric matrix
            A(i, l) = a_il*c - a_ik*s; // 3)
            A(l, i) = A(i, l);
        }

        // The new eigenvectors
        double r_ik = R(i, k);
        double r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
    }
}


void jacobiMethod(arma::mat &A, arma::mat &R, int n) {

    // Tolerance for the non-diagonals
    double eps = 1.0E-8;

    // Defining maximum number of iterations
    int maxiter=n*n*n;
    int iterations = 0;
    int k = 0;
    int l = 0;
    double max_nondiagonal = maxOffDiagonals(A, k, l, n);
    // The Jacobi rotation algorithm
    while (max_nondiagonal > eps && iterations <=maxiter ) {
        jacobiRotation(A, R, k, l, n);

        cout << "With max diagonal: "<< max_nondiagonal << endl;
        cout << "A" << endl;
        cout << A << endl;
        cout << "R" << endl;
        cout << R << endl;
        max_nondiagonal = maxOffDiagonals(A, k, l, n);
        iterations++;
    }
}
