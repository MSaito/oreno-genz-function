#pragma once
#ifndef MAKE_PARAMETERS
#define MAKE_PARAMETERS

#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "mt19937_64.hpp"

void printArray(std::string str, double array[], int len)
{
    using namespace std;
    cout << "#" << str << " = {";
    for (int i = 0; i < len; i++) {
        cout << scientific << setprecision(18) << array[i] << " ";
    }
    cout << "};" << endl;
}

void makeParameter(int func_index, int dim, int seed, int original,
                   double a[], double b[], double alpha[], double beta[],
                   bool verbose)
{
    mt19937_64 mt(seed);

    if (original ==  1) {
        double exn;
        double dfclt;
        double expnts[] = {1.5, 2.0, 2.0, 1.0, 2.0, 2.0};
        //double difclt[] = {110.0, 600.0, 600.0, 100.0, 150.0, 100.0};
        double difclt[] = {0.9, 0.725, 0.185, 0.703, 2.04, 0.43};

        exn = expnts[func_index - 1];
        dfclt = difclt[func_index - 1] * dim;

        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
        }
        double total = 0;
        for (int i = 0; i < dim; i++) {
            if (func_index == 1 || func_index == 2) {
                alpha[i] = 1.5 - mt.getDouble01();
            } else {
                alpha[i] = mt.getDouble01();
            }
            beta[i] = mt.getDouble01();
            total += alpha[i];
        }
        double dfact = total * pow(dim, exn) / dfclt;
        //double dfact = total / dfclt;
        //cout << "#dfact = " << dfact << endl;
        if (func_index != 1 && func_index != 2) {
            for (int i = 0; i < dim; i++) {
                alpha[i] = alpha[i] / dfact;
            }
        }
#if 1
        if ((func_index == 1) || (func_index == 3)) {
            for (int i = 0; i < dim; i++) {
                b[i] = alpha[i];
            }
        }
#endif
        if (func_index == 6) {
            for (int i = 2; i < dim; i++) {
                beta[i] = 1.0;
                //beta[i] = 2 / M_PI;
            }
        }
    } else if (original == -1) {
        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
            alpha[i] = 1.0;
            if (func_index == 3) {
                beta[i] = mt.getDouble01();
            } else if (func_index == 6) {
                beta[i] = 2 / M_PI;
            } else {
                beta[i] = 1.0;
            }
        }
    } else {
        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
            switch (func_index) {
            case 1:
                alpha[i] = 1.5 - mt.getDouble01();
                if (i == 0) {
                    beta[0] = 1.5 - mt.getDouble01();
                } else {
                    beta[i] = 0;
                }
                break;
            case 2:
                alpha[i] = 1.5 - mt.getDouble01();
                beta[i] = mt.getDouble01();
                break;
            case 3:
                alpha[i] = 1.5 - mt.getDouble01();
                beta[i] = 0;
                break;
            case 4:
                alpha[i] = 1.5 - mt.getDouble01();
                beta[i] = mt.getDouble01();
                break;
            case 5:
                alpha[i] = mt.getDouble01();
                beta[i] = mt.getDouble01();
                break;
            case 6:
            default:
                alpha[i] = mt.getDouble01() / dim;
                if (i < 2) {
                    beta[i] = mt.getDouble01();
                } else {
                    beta[i] = 1.0;
                }
            }
        }
    }
    if (verbose) {
        std::cout << "original = " << std::dec << original << std::endl;
        printArray("a", a, dim);
        printArray("b", b, dim);
        printArray("alpha", alpha, dim);
        printArray("beta", beta, dim);
    }
}


#endif // MAKE_PARAMETERS
