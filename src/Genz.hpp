#pragma once
#ifndef GENZ_HPP
#define GENZ_HPP

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace GenzNS {
    class GenzFunction {
    public:
        virtual ~GenzFunction(){};
        virtual double operator()(int ndim,
                                  const double z[],
                                  const double alpha[],
                                  const double beta[]) const = 0;
        virtual double integral(int ndim,
                                const double a[],
                                const double b[],
                                const double alpha[],
                                const double beta[]) const = 0;
        virtual std::string name() const = 0;
        virtual double difficulty() const = 0;
        virtual void setTestParams(int dim, double a[], double b[],
                                   double alpha[], double beta[]) const = 0;
    };

    class Oscillatory: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    class ProductPeak: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    class CornerPeak: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    class Gaussian: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    class C0Function: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    class Discontinuous: public GenzFunction {
    public:
        double operator()(int ndim, const double z[],
                          const double alpha[], const double beta[]) const;
        double integral(int ndim, const double a[], const double b[],
                        const double alpha[], const double beta[]) const;
        std::string name() const;
        double difficulty() const;
        void setTestParams(int dim, double a[], double b[],
                           double alpha[], double beta[]) const;
    };

    template<typename D>
    double integral(const GenzFunction& func, D& digitalNet, int count, int dim,
                    const double alpha[], const double beta[])
    {
        using namespace std;
        double a[dim];
        double b[dim];
        double alpha1[dim];
        double beta1[dim];
        double total = 0;
        for (int i = 0; i < dim; i++) {
            a[i] = 0;
            b[i] = 1.0;
            alpha1[i] = alpha[i];
            beta1[i] = beta[i];
            total += alpha1[i];
        }
        double dfclt = func.difficulty();
        double dfact = total / (dfclt * dim);
        for (int i = 0; i < dim; i++) {
            alpha1[i] = alpha1[i] / dfact;
        }
        func.setTestParams(dim, a, b, alpha1, beta1);
#if defined(DEBUG)
        cout << "after setTestParams" << endl;
        cout << "a = (";
        for (int i = 0; i < dim; i++) {
            cout << a[i] << " ";
        }
        cout << ")" << endl;
        cout << "b = (";
        for (int i = 0; i < dim; i++) {
            cout << b[i] << " ";
        }
        cout << ")" << endl;
        cout << "alpha = (";
        for (int i = 0; i < dim; i++) {
            cout << alpha1[i] << " ";
        }
        cout << ")" << endl;
        cout << "beta = (";
        for (int i = 0; i < dim; i++) {
            cout << beta1[i] << " ";
        }
        cout << ")" << endl;
#endif
        //digitalNet.pointInitialize();
        double sum = 0;
        //for (int i = 1; i < count; i++) {
        for (int i = 0; i < count; i++) {
            sum += func(dim, digitalNet.getPoint(), alpha1, beta1);
            digitalNet.nextPoint();
        }
        sum = sum / count;
        double expected = func.integral(dim, a, b, alpha1, beta1);
#if defined(DEBUG)
        cout << "expected = " << expected << endl;
        cout << "sum = " << sum << endl;
#endif
        return abs(expected - sum);
    }

    template<typename D>
    void integralAll(const GenzFunction& func, D& digitalNet,
                     int start_m, int stop_m, int dim,
                     const double alpha[], const double beta[])
    {
        using namespace std;
        cout << digitalNet.getName() << endl;
        cout << func.name() << endl;
        for (int m = start_m; m <= stop_m; m++) {
            int count = 1 << m;
            double err = integral(func, digitalNet, count, dim, alpha, beta);
            cout << "m = " << dec << m << ": error = " << err << endl;
        }
    }
}

#endif // GENZ_HPP
