#include <iostream>
#include <math.h>
#include <armadillo>

#include <jacobimethod.h>

using namespace std;

int main()
{
    int n = 3;

    arma::mat A(n,n);
    arma::mat R; R.eye(n,n);

    A 	<< 4 << 1 << 0 << arma::endr
        << 1 << -1 << 0 << arma::endr
        << 0 << 0 << 5 << arma::endr;
    //maxOffDiagonals(A, k, l, n);
    cout << A << endl;
//    exit(1);
    jacobiMethod(A, R, n);

}
