#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "adjust_parameters.h"
#include "testpack.h"
#include "make_parameters.h"

#define DEBUG

using namespace std;

class param_value {
public:
    param_value(int dim) {
        this->dim = dim;
        alpha = new double[dim];
        beta = new double[dim];
    }
    void set(double valuep, const double alphap[], const double betap[]) {
        value = valuep;
        for (int i = 0; i < dim; i++) {
            alpha[i] = alphap[i];
            beta[i] = betap[i];
        }
    }
    double get(double alphap[], double betap[]) {
        for (int i = 0; i < dim; i++) {
            alphap[i] = alpha[i];
            betap[i] = beta[i];
        }
        return value;
    }

    ~param_value() {
        delete[] alpha;
        delete[] beta;
    }
    double value;
    double * alpha;
    double * beta;
private:
    int dim;
};

void print_parameters(bool verbose, int dim, double a[], double b[],
                      double alpha[], double beta[]);

double adjustParameter(int func_index, int dim,
                       double a[], double b[],
                       double alpha[], double beta[], bool verbose)
{
    param_value pre(dim);
    double expected = genz_integral(func_index, dim, a, b, alpha, beta);
    double ratio = 0.9;
    bool large = true;
    double alae = abs(1.0 - abs(log2(abs(expected))));
#if defined(DEBUG)
    cout << "expected = " << expected << endl;
    cout << "alae = " << alae << endl;
#endif
#if 0
    if (func_index == 1) {
        print_parameters(verbose, dim, a, b, alpha, beta);
        return expected;
    }
#endif
    if (alae < 1.0) {
        print_parameters(verbose, dim, a, b, alpha, beta);
        return expected;
    }
    pre.set(expected, alpha, beta);
    if (expected > 1.0) {
        large = true;
        ratio = 0.9;
    } else {
        large = false;
        ratio = 1.05;
    }
    for (int i = 0; i < 1000; i++) {
        if (func_index == 1 || func_index == 6) {
            for (int j = 0; j < dim; j++) {
                beta[j] *= ratio;
            }
        } else {
            for (int j = 0; j < dim; j++) {
                alpha[j] *= ratio;
            }
        }
        expected = genz_integral(func_index, dim, a, b, alpha, beta);
        alae = abs(1.0 - abs(log2(abs(expected))));
#if defined(DEBUG)
        cout << "expected = " << expected << endl;
        cout << "alae = " << alae << endl;
#endif
        if (alae < 1.0) {
            break;
        }
        if (alae <= abs(1.0 - abs(log2(abs(pre.value))))) {
            pre.set(expected, alpha, beta);
            continue;
        } else {
            // 前回より悪化した
#if defined(DEBUG)
            cout << "get worse " << endl;
#endif
            break;
#if 0
            if (large) {
                ratio = 1.1;
                large = false;
            } else {
                ratio = 0.9;
                large = true;
            }
#endif
        }
    }
    if (alae > abs(1.0 - abs(log2(abs(pre.value))))) {
        expected = pre.get(alpha, beta);
    }
    print_parameters(verbose, dim, a, b, alpha, beta);
    return expected;
}

void print_parameters(bool verbose, int dim, double a[], double b[],
                      double alpha[], double beta[])
{
    if (verbose) {
        cout << "# adjust parameters" << endl;
        printArray("a", a, dim);
        printArray("b", b, dim);
        printArray("alpha", alpha, dim);
        printArray("beta", beta, dim);
    }

}
