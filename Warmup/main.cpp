/*
* Warmup Exercise
*/
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <iomanip>
using namespace std;

double f(double x){

    return atan(x);
}

double f2c(double x, double h){

    return ( f(x+h)-f(x) ) /h;
}

double f3c(double x, double h){

    return ( f(x+h)-f(x-h) )/(2*h);
}

int main() {

    int N = 10;
    double x = sqrt(2);

    double *h = new double[N];
    for (int i=0; i<N; i++) {
        h[i]=pow(10,-i);
    }
    double f2cs[N];
    double f3cs[N];

    ofstream ofs ("../Warmup/data.txt", ofstream::out);

    for (int i=0; i<N; i++){
        f2cs[i] = f2c(x,h[i]);
        f3cs[i] = f3c(x,h[i]);

        ofs << setprecision(16) << f2cs[i] << " ";
        ofs << setprecision(16)  << f3cs[i] << " ";
        ofs << setprecision(16)  << h[i] << " "<< endl;


    }


    ofs.close();

    //    double h_best = h[min_element(abs(x_exact-f2cs))];
    //    double f2c_best = f2c(x,h_best);
    //    double f3c_best = f3c(x,h_best);

    //    cout << "Exact: " << x_exact << endl;
    //    cout << "2c:" << f2c_best << "with error:" << min(abs(x_exact-f2cs)) << endl;
    //    cout << "3c:" << f3c_best << "with error:" << min(abs(x_exact-f3cs)) << endl;
    //    cout << "Best h:" << h_best << endl;


}
