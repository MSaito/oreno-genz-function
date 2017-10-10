#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "adjust_parameters.h"
#include "testpack.h"
#include "make_parameters.h"

//#define DEBUG

using namespace std;

namespace {

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
    double adjustOscillatory(int dim,
                             double a[], double b[],
                             double alpha[], double beta[], bool verbose);

}

double adjustParameter(int func_index, int dim,
                       double a[], double b[],
                       double alpha[], double beta[], bool verbose)
{
    if (func_index == 1) {
        return adjustOscillatory(dim, a, b, alpha, beta, verbose);
    }
    param_value pre(dim);
    double expected = genz_integral(func_index, dim, a, b, alpha, beta);
    double ratio = 0.9;
    double step = 0.0002;

    //double alae = abs(1.0 - abs(log2(abs(expected))));
    double alae = abs(log2(abs(expected)));
#if defined(DEBUG)
    cout << "expected = " << expected << endl;
    cout << "alae = " << alae << endl;
#endif

    if (alae < 1.0) {
        print_parameters(verbose, dim, a, b, alpha, beta);
        return expected;
    }
    pre.set(expected, alpha, beta);
    if (expected > 1.0) {
        ratio = 0.9;
        step = -step;
    } else {
        ratio = 1.05;
    }
    if (dim >= 2) {
        ratio = 1.0;
    } else {
        step = 0;
    }
    for (int i = 0; i < 1000; i++) {
        if (func_index == 1 || func_index == 6) {
            for (int j = 0; j < dim; j++) {
                beta[j] *= ratio + step;
            }
        } else {
            for (int j = 0; j < dim; j++) {
                alpha[j] *= ratio + step;
            }
        }
        expected = genz_integral(func_index, dim, a, b, alpha, beta);
        //alae = abs(1.0 - abs(log2(abs(expected))));
        alae = abs(log2(abs(expected)));
#if defined(DEBUG)
        cout << "expected = " << expected << endl;
        cout << "alae = " << alae << endl;
#endif
        if (alae < 1.0) {
            break;
        }
        //if (alae <= abs(1.0 - abs(log2(abs(pre.value))))) {
        if (alae <= abs(log2(abs(pre.value)))) {
            pre.set(expected, alpha, beta);
            continue;
        } else {
            // 前回より悪化した
#if defined(DEBUG)
            cout << "get worse " << endl;
#endif
            break;
        }
    }
    //if (alae > abs(1.0 - abs(log2(abs(pre.value))))) {
    if (alae > abs(log2(abs(pre.value)))) {
        expected = pre.get(alpha, beta);
    }
    print_parameters(verbose, dim, a, b, alpha, beta);
    return expected;
}

namespace {
    double adjustOscillatory(int dim,
                             double a[], double b[],
                             double alpha[], double beta[], bool verbose)
    {
#if defined(DEBUG)
        cout << "in adjustOscillatory" << endl;
#endif
        param_value minparam(dim);
        double minalae = INFINITY;
        double ratio = 0.9;
        for (int i = 0; i < 10; i++) {
            double sum = 0;
            for (int j = 0; j < dim; j++) {
                sum += alpha[j] / 2;
            }
            beta[0] = 1.0 - sum / (2 * M_PI);
            double expected = genz_integral(1, dim, a, b, alpha, beta);
            double alae = abs(1.0 - abs(log2(abs(expected))));
#if defined(DEBUG)
            cout << "expected = " << expected << endl;
#endif
            if (minalae > alae) {
                minalae = alae;
                minparam.set(expected, alpha, beta);
            }
            if (expected > 0.1) {
                break;
            }
            for (int j = 0; j < dim; j++) {
                alpha[j] = alpha[j] * ratio;
                if (alpha[j] == 0) {
                    alpha[j] = 0.00001;
                    break;
                }
            }
        }
        double expected = minparam.get(alpha, beta);
        print_parameters(verbose, dim, a, b, alpha, beta);
#if defined(DEBUG)
        cout << "expected = " << expected << endl;
#endif
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
}
