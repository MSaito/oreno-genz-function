#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "mt19937_64.hpp"
#include "make_parameters.h"
#include "make_wafomc_parameters.h"
#include "cvmean.h"

using namespace std;

void makeWafomParameter(int func_index, int dim, int seed,
                        double a[], double b[], double alpha[], double beta[],
                        bool verbose, double mag)
{
    double c = calc_c_for_cvmean(dim, 64);
    double twoc = pow(2.0, c);
    mt19937_64 mt(seed);
    for (int i = 0; i < dim; i++) {
        a[i] = 0.0;
        b[i] = 1.0;
    }
    double maxa = 1.0;
    switch (func_index) {
    case 2:
        maxa = 3.0 * mag;
        break;
    case 1:
    case 5:
    case 6:
        maxa = twoc * mag;
        break;
        //case 2:
    case 4:
    default:
        maxa = sqrt(twoc) * mag;
    }
    double total = 0;
    for (int i = 0; i < dim; i++) {
        alpha[i] = maxa * mt.getDouble01();
        beta[i] = mt.getDouble01();
        total += alpha[i] / 2;
    }
    if (func_index == 6) {
        for (int i = 2; i < dim; i++) {
            beta[i] = 1.0;
        }
    }
    if (func_index == 1) {
        beta[0] = 1.0 - total / (2 * M_PI);
    }
    if (verbose) {
        cout << "#wafom mag = " << mag << endl;
        cout << "#wafom c = " << c << endl;
        cout << "#2^c = " << twoc << endl;
        cout << "#max a = " << maxa << endl;
        //printArray("a", a, dim);
        //printArray("b", b, dim);
        printArray("alpha", alpha, dim);
        printArray("beta", beta, dim);
    }
}
