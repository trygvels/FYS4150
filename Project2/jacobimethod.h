#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

double maxOffDiagonals(arma::mat &A, int &k, int &l, int n);
void jacobiRotation(arma::mat &A, arma::mat &R, int &k, int &l, int n);
void jacobiMethod(arma::mat &A, arma::mat &R, int n);
void jacobiEigTest(arma::mat &A, arma::mat &R, int n);
arma::mat constructA(int n,int interacting);

#endif // JACOBIMETHOD_H
