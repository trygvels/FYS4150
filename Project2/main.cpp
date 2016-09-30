#include <iostream>
#include <math.h>
#include <armadillo>

#include <jacobimethod.h>
#include <assert.h>

using namespace std;

int main()
{
    int n = 300;
    // If it is interacting. interacting=1
    int interacting = 1;

    arma::mat A(n,n);
    arma::mat R; R.eye(n,n);
    A = constructA(n,interacting);
    //cout <<A<<endl;

    jacobiMethod(A, R, n);

//    A   << 4 << 0 << 0 << arma::endr
//        << 0 << 4 << -2 << arma::endr
//        << 0 << -2 << 1 << arma::endr;
//    This matrix has eigvals  0,4,5

//    A << 4 << 1 << 1 << 0 <<-1 << arma::endr
//      << 1 << 4 << 0 << 1 <<-1 << arma::endr
//      << 1 << 0 << 3 << 0 << 0 << arma::endr
//      << 0 << 1 << 0 << 3 << 0 << arma::endr
//      <<-1 <<-1 << 0 << 0 << 3 << arma::endr;

//    checking eigenvalues for simple matrix
//    jacobiEigTest(A, R, n);
//    cout << "success" << endl;

}

