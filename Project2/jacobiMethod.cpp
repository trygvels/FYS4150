#include <iostream>
#include <math.h>
#include <armadillo>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <iostream>

using namespace std;

arma::mat constructA(double &rho_min, double &rho_max, int n,int interacting,double Omega_r){
    //Step length
    double h = (rho_max-rho_min)/(n+1);

    //Defining the constants in diagonal elements
    arma::mat A = arma::zeros<arma::mat>(n,n);
    double d = 2/(h*h);
    double e = -1/(h*h);

    //This is 0 so not included
    //int l = 0;
    //double OrbitalFactor = l*(l+1);

    //Empty potential
    double V;

    //Creating the tridiagonal matrix
    for (int i = 0; i < (n); i++) {
        double rho = (i+1)*h;

        //Choosing potential
        if (interacting==true){
            V=Omega_r*Omega_r*rho*rho + 1/(rho);
        }else{
            V=rho*rho;
        }

        //Setting diagonal elements
        A(i,i)   = d + V;
        if(i<n-1){
            A(i,i+1) = A(i+1,i) = e;
        }
    }

    //Return constructed matrix
    return A;
}

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

    //cout << "Largest off-diagonal element is: " << max << endl;

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
            A(i, l) = a_il*c + a_ik*s; // 3)
            A(l, i) = A(i, l);
        }

        // The new eigenvectors
        double r_ik = R(i, k);
        double r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
    }

}

ofstream ofile;
void output(double rho_max , double rho_min, int n, arma::mat &R, double Omega_r,int lowestvalueindex){
  ofile.open("dataproj2e.txt");
  ofile <<"Omega_r = " << setw(15) << setprecision(8) << Omega_r << endl;
  ofile <<"rho_min = " << setw(15) << setprecision(8) << rho_min << endl;
  ofile <<"rho_max = " << setw(15) << setprecision(8) << rho_max << endl;
  ofile <<"Number of steps = " << setw(15) << n << endl;
  ofile << "Eigenvector corresponding to lowest eigenvalue:" << endl;
  for (int i = 0; i < n; i++) {
          ofile << R(lowestvalueindex,i) << endl;
      }
      ofile <<endl;

  ofile.close();
}


void jacobiMethod(arma::mat &A, arma::mat &R, int n) {
    double rho_min = 0.0;
    double rho_max = 5.0;

    // If it is interacting then true
    bool interacting = true;

    //Choosing angular freq
    double Omega_r = 0.01;

    // Constructing matrix A
    A = constructA(rho_min, rho_max,n,interacting,Omega_r);

    // Tolerance for the non-diagonals
    double eps = 1.0e-8;

    // Defining maximum number of iterations
    int maxiter=n*n*n;
    int iterations = 0;
    int k = 0;
    int l = 0;
    double max_nondiagonal = maxOffDiagonals(A, k, l, n);
    // The Jacobi rotation algorithm
    while (max_nondiagonal > eps && iterations <=maxiter ) {
        jacobiRotation(A, R, k, l, n);


//        cout << "With max diagonal: "<< max_nondiagonal << endl;
//        cout << "A" << endl;
//        cout << A << endl;
//        cout << "R" << endl;
//        cout << R << endl;

        max_nondiagonal = maxOffDiagonals(A, k, l, n);
        iterations++;
    }

    //sorting eigenvalues and printing
    arma::vec lambda = A.diag();

    int lowestvalueindex = lambda.index_min();

    lambda=sort(lambda);

    //Printing 5 lowest eigenvalues
    int i=0;
    while(i<3){
        cout <<"Eigenvalue #" << i << ": "<< lambda(i)<<endl;
        i++;
    }

    output(rho_max, rho_min, n, R,Omega_r,lowestvalueindex);

}



/*
The following functions are tests
  */



void jacobiEigTest(arma::mat &A, arma::mat &R, int n){

    //Defining the exact eigenvalues for a chosen matrix
    double * lambda_exact = new double[n];
    lambda_exact[0] = 2.0;
    lambda_exact[1] = 2.0;
    lambda_exact[2] = 3.0;
    lambda_exact[3] = 4.0;
    lambda_exact[4] = 6.0;

    // Defining tolerance
    double eps = 1.0E-8;

    // fill A with eigenvalues
    jacobiMethod(A, R, n);

    bool s1=false;
    bool s2=false;
    bool s3=false;
    bool s4=false;
    bool s5=false;

    //Checking difference between exact and numerical to tolerance
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            s1 = s1 or ( abs( A(i,j) - lambda_exact[0]) < eps);
            s2 = s2 or ( abs( A(i,j) - lambda_exact[1]) < eps);
            s3 = s3 or ( abs( A(i,j) - lambda_exact[2]) < eps);
            s4 = s4 or ( abs( A(i,j) - lambda_exact[3]) < eps);
            s5 = s5 or ( abs( A(i,j) - lambda_exact[4]) < eps);
        }
    }

    // Assertion
    assert(s1);
    assert(s2);
    assert(s3);
    cout<<"successfully calculated eigenvalues"<<endl;
}

void jacobiOrthogTest(arma::mat &A,arma::mat &R,int n){
    jacobiMethod(A, R, n);
    double eps=1e-14;
    arma::mat R_trans = arma::trans(R);
    arma::mat I; I.eye(n,n);
    bool s4=false;
    arma::mat compare = abs(R_trans*R-I);

    if(arma::all(all(compare<eps,1))){
            s4=true;
    }

    assert(s4);
    cout << "Orthogonality conserved"<<endl;
}

void jacobiMaxoffTest(arma::mat &A,int n){
    int k=0;
    int l=0;
    double max = maxOffDiagonals(A, k, l, n);

    bool s5=false;
    if(max==1){
        s5=true;
    }
    assert(s5);
    cout << "Max value found"<<endl;
}

void tests(arma::mat &A,arma::mat &R,int n){
    jacobiMaxoffTest(A,n);
    jacobiEigTest(A, R, n);
    jacobiOrthogTest(A,R,n);
    cout << "All tests completed" << endl;

}
