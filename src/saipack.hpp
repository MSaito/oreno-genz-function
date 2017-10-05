#pragma once
#ifndef SAIPACK_HPP
#define SAIPACK_HPP

#include <string>
#include <random>
#include <cmath>

class Saipack {
public:
    virtual ~Saipack() {};
    virtual double operator()(const double x[]) = 0;
    virtual void setParam(int dim, double ap[], double bp[]) = 0;
    virtual double expected(int dim, double ap[], double bp[]) = 0;
    virtual const std::string getName() = 0;
    virtual void makeParameter(int type, int dim, std::mt19937_64& mt,
                               double a[], double b[], bool verbose) = 0;
};

static inline void printParameter(int dim, double a[], double b[]) {
    std::cout << "# a = {";
    for (int i = 0; i < dim; i++) {
        std::cout << a[i] << ",";
    }
    std::cout << "}" << std::endl;
    std::cout << "# b = {";
    for (int i = 0; i < dim; i++) {
        std::cout << b[i] << ",";
    }
    std::cout << "}" << std::endl;
}

class AddSai : Saipack {
public:
    AddSai() {
        a = NULL;
        b = NULL;
    }
    ~AddSai() {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
    }
    void setParam(int dim, double ap[], double bp[]) {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
        s = dim;
        a = new double[dim];
        b = new double[dim];
        for (int i = 0; i < dim; i++) {
            a[i] = ap[i];
            b[i] = bp[i];
        }
    }
    double operator()(const double x[]) {
        double sum = 0;
        for (int i = 0; i < s; i++) {
            sum += x[i] * x[i] + a[i] * x[i] + b[i];
        }
        return sum;
    }
    double expected(int dim, double ap[], double bp[]) {
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            sum += 1 / 3 + ap[i] / 2 + bp[i];
        }
        return sum;
    }

    void makeParameter(int type, int dim, std::mt19937_64& mt,
                       double a[], double b[], bool verbose) {
        switch (type) {
        case 0:
            for (int i = 0; i < dim; i++) {
                a[i] = 1 / dim;
                b[i] = 1 / dim;
            }
            break;
        case 1:
        default:
            std::uniform_real_distribution<double> unif01(0, 1);
            for (int i = 0; i < dim; i++) {
                a[i] = unif01(mt) / dim;
                b[i] = unif01(mt) / dim;
            }
            break;
        }
        if (verbose) {
            printParameter(dim, a, b);
        }
        return;
    }

    const std::string getName() {
        return std::string("AddSai");
    }
private:
    int s;
    double * a;
    double * b;
};

class MulSai : Saipack {
public:
    MulSai() {
        a = NULL;
        b = NULL;
    }

    ~MulSai() {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
    }
    void setParam(int dim, double ap[], double bp[]) {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
        s = dim;
        a = new double[dim];
        b = new double[dim];
        for (int i = 0; i < dim; i++) {
            a[i] = ap[i];
            b[i] = bp[i];
        }
    }
    double operator()(const double x[]) {
        double prod = 1.0;
        for (int i = 0; i < s; i++) {
            prod *= x[i] * x[i] * x[i]
                + a[i] * x[i] * x[i] + b[i];
        }
        return prod;
    }
    double expected(int dim, double ap[], double bp[]) {
        double prod = 1.0;
        for (int i = 0; i < dim; i++) {
            prod *= 1 / 4 + ap[i] / 3 + bp[i];
        }
        return prod;
    }

    void makeParameter(int type, int dim, std::mt19937_64& mt,
                       double a[], double b[], bool verbose) {
        switch (type) {
        case 0:
            for (int i = 0; i < dim; i++) {
                a[i] = 1;
                b[i] = 5 / 12;
            }
            break;
        case 1:
        default:
            std::uniform_real_distribution<double> unif01(0, 1);
            for (int i = 0; i < dim; i++) {
                a[i] = unif01(mt);
                b[i] = 1 - 1 / 4 - a[i];
            }
            break;
        }
        if (verbose) {
            printParameter(dim, a, b);
        }
        return;
    }

    const std::string getName() {
        return std::string("MulSai");
    }

private:
    int s;
    double * a;
    double * b;
};

class SinSai : Saipack {
public:
    SinSai() {
        a = NULL;
        b = NULL;
    }

    ~SinSai() {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
    }
    void setParam(int dim, double ap[], double bp[]) {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
        s = dim;
        a = new double[dim];
        b = new double[dim];
        for (int i = 0; i < dim; i++) {
            a[i] = ap[i];
            b[i] = bp[i];
        }
    }
    double operator()(const double x[]) {
        double prod = 1.0;
        for (int i = 0; i < s; i++) {
            prod *= sin(a[i] * x[i] + 2 * M_PI * b[i]);
        }
        return 1 + prod;
    }
    double expected(int dim, double ap[], double bp[]) {
        double prod = 1.0;
        for (int i = 0; i < dim; i++) {
            prod *= cos(2 * M_PI * bp[i]) * cos(ap[i]) / ap[i]
                + sin(2 * M_PI * bp[i]) * sin(ap[i]) / ap[i]
                - cos(2 * M_PI * bp[i]) / ap[i];
        }
        return 1 + prod;
    }

    void makeParameter(int type, int dim, std::mt19937_64& mt,
                       double a[], double b[], bool verbose) {
        switch (type) {
        case 0:
            for (int i = 0; i < dim; i++) {
                a[i] = 1;
                b[i] = 0.1 * i;
            }
            break;
        case 1:
        default:
            std::uniform_real_distribution<double> unif01(0, 1);
            for (int i = 0; i < dim; i++) {
                a[i] = unif01(mt);
                b[i] = unif01(mt);
            }
            break;
        }
        double ex = expected(dim, a, b);
        if (ex > 100 || ex < 0.01) {
            ex = pow(ex, 1 / dim);
        }
        for (int i = 0; i < dim; i++) {
            a[i] = a[i] * ex;
        }
        if (verbose) {
            printParameter(dim, a, b);
        }
        return;
    }

    const std::string getName() {
        return std::string("SinSai");
    }

private:
    int s;
    double * a;
    double * b;
};

/**
 * f = Prod_{i=1}^s(a_i(x_i + 1)(x_i - 2)(x_i - 0.2b_i)(x_i - 0.4b_i)(x_i - 0.6b_i)(x_i - 0.8b_i)
 * 1 + 0という値になりがち
 * a_i は小さいと積分値が小さくなるはず。絶対値が1より大きい方がよいと思われる
 * b_i 負：たぶん意味なし
 * b_i < 1.0 前より
 * b_i > 1.0 山数がへっていく
 * b_i >= 5 谷か山
 */
class PolSai : Saipack {
public:
    PolSai() {
        a = NULL;
        b = NULL;
    }

    ~PolSai() {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
    }
    void setParam(int dim, double ap[], double bp[]) {
        if (a != NULL) {
            delete[] a;
            delete[] b;
        }
        s = dim;
        a = new double[dim];
        b = new double[dim];
        for (int i = 0; i < dim; i++) {
            a[i] = ap[i];
            b[i] = bp[i];
        }
    }
    double operator()(const double x[]) {
        double prod = 1.0;
        for (int i = 0; i < s; i++) {
            prod *= a[i];
            prod *= (x[i] + 1) * (x[i] - 2);
            for (int j = 0; j < 4; j++) {
                prod *= (x[i] - 0.2 * b[i] * (j + 1));
            }
        }
        return 1 + prod;
    }

    double expected(int dim, double ap[], double bp[]) {
        double prod = 1.0;
        for (int i = 0; i < dim; i++) {
            double q = bp[i];
            prod *= ap[i];
            prod *= - 48 * q * q * q * q / 625
                - 2 / 625 * q * q * q * (-125 + 6 * q)
                + 2 / 1875 * q * (-875 + 125 * q + 12 * q * q)
                - 1 / 20 * q * (-20 + 7 * q + 2 * q * q)
                + 1 / 25 * (-10 + 10 * q + 7 * q * q)
                + 1 / 6 * (-1 + 2 * q)
                + 1 / 7;
        }
        return 1 + prod;
    }

    void makeParameter(int type, int dim, std::mt19937_64& mt,
                       double a[], double b[], bool verbose) {
        switch (type) {
        case 0:
            for (int i = 0; i < dim; i++) {
                a[i] = 1;
                b[i] = 0.1 * i;
            }
            break;
        case 1:
        default:
            std::uniform_real_distribution<double> unif01(0, 1);
            for (int i = 0; i < dim; i++) {
                a[i] = unif01(mt);
                b[i] = unif01(mt);
            }
            break;
        }
        double ex = expected(dim, a, b);
        if (ex > 100 || ex < 0.01) {
            ex = pow(ex, 1 / dim);
        }
        for (int i = 0; i < dim; i++) {
            a[i] = a[i] * ex;
        }
        if (verbose) {
            printParameter(dim, a, b);
        }
        return;
    }

    const std::string getName() {
        return std::string("PolSai");
    }

private:
    int s;
    double * a;
    double * b;
};

#endif // SAIPACK_HPP
