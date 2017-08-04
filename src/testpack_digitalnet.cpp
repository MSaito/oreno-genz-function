#define _USE_MATH_DEFINES
#include <math.h>
#include <cerrno>
#include <getopt.h>
#include <cstdlib>
#include <string>
#include <random>
#include <MCQMCIntegration/DigitalNet.h>
#include "testpack.h"
#include "kahan.hpp"

using namespace std;
using namespace MCQMCIntegration;

namespace {
    struct cmd_opt_t {
        uint32_t s_dim;
        uint32_t start_m;
        uint32_t end_m;
        uint32_t seed;
        int genz_no;
        int dn_id;
        int rmse;
    };
    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    template<typename D>
    double integral(int func_index, D& digitalNet, int count, int dim,
                    mt19937_64& mt, int rmse);
}

int main(int argc, char *argv[]) {
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    DigitalNetID dnid = static_cast<DigitalNetID>(opt.dn_id);
    mt19937_64 mt(opt.seed);
    cout << "#" << genz_name(opt.genz_no) << endl;
    cout << "#" << getDigitalNetName(opt.dn_id) << endl;
    if (opt.rmse > 0) {
        cout << "#m, abs err, log2(RMSE[" << dec << opt.rmse << "])" << endl;
    } else {
        cout << "#m, abs err, log2(err)" << endl;
    }
    for (uint32_t m = opt.start_m; m <= opt.end_m; m++) {
        DigitalNet<uint64_t> dn(dnid, opt.s_dim, m);
        int count = 1 << m;
        double error = integral(opt.genz_no, dn, count, opt.s_dim,
                                mt, opt.rmse);
        cout << dec << m << "," << error << "," << log2(error) << endl;
    }
    return 0;
}

namespace {
    void cmd_message(const string& pgm)
    {
        cout << pgm << "-s s_dim -m start_m -M end_m -S seed -g genz_no"
             << " -d digitalnet_id " << endl;
    }

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv)
    {
        int c;
        bool error = false;
        string pgm = argv[0];
        static struct option longopts[] = {
            {"s-dim", required_argument, NULL, 's'},
            {"start-m", required_argument, NULL, 'm'},
            {"end-m", required_argument, NULL, 'M'},
            {"seed", required_argument, NULL, 'S'},
            {"genz-no", required_argument, NULL, 'g'},
            {"digitalnet-id", required_argument, NULL, 'd'},
            {"rmse", optional_argument, NULL, 'r'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.start_m = 0;
        opt.end_m = 0;
        opt.seed = 1;
        opt.genz_no = 0;
        opt.dn_id = 0;
        opt.rmse = 0;
        errno = 0;
        for (;;) {
            c = getopt_long(argc, argv, "s:m:M:S:g:d:r:", longopts, NULL);
            if (error) {
                break;
            }
            if (c == -1) {
                break;
            }
            switch (c) {
            case 's':
                opt.s_dim = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "s_dim should be a number" << endl;
                    error = true;
                }
                break;
            case 'm':
                opt.start_m = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "start_m should be a number" << endl;
                    error = true;
                }
                break;
            case 'M':
                opt.end_m = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "end_m should be a number" << endl;
                    error = true;
                }
                break;
            case 'S':
                opt.seed = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "seed should be a number" << endl;
                    error = true;
                }
                break;
            case 'g':
                opt.genz_no = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "genz_no should be a number" << endl;
                    error = true;
                }
                break;
            case 'd':
                opt.dn_id = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "digitalnet_id should be a number" << endl;
                    error = true;
                }
                break;
            case 'r':
                if (optarg == NULL) {
                    opt.rmse = 100;
                } else {
                    opt.rmse = strtoul(optarg, NULL, 10);
                }
                if (errno) {
                    cout << "digitalnet_id should be a number" << endl;
                    error = true;
                }
                break;
            case '?':
            default:
                error = true;
                break;
            }
        }
        if (error) {
            cmd_message(pgm);
            return false;
        }
        return true;
    }

    template<typename D>
    double integral(int func_index, D& digitalNet, int count, int dim,
                    mt19937_64& mt, int rmse)
    {
        double expected;
        double a[dim];
        double b[dim];
        double alpha[dim];
        double beta[dim];
        uniform_real_distribution<double> unif01(0.0, 1.0);

#if defined(SAME_AS_GENZ) || 0
        double exn;
        double dfclt;
        double expnts[] = {1.5, 2.0, 2.0, 1.0, 2.0, 2.0};
        double difclt[] = {110.0, 600.0, 100.0, 150.0, 100.0};

        exn = expnts[func_index - 1];
        dfclt = difclt[func_index - 1];
        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
        }
        double total = 0;
        for (int i = 0; i < dim; i++) {
            alpha[i] = unif01(mt);
            beta[i] = unif01(mt);
            total += alpha[i];
        }
        double dfact = total * pow(dim, exn) / dfclt;
        for (int i = 0; i < dim; i++) {
            alpha[i] = alpha[i] / dfact;
        }
        if ((func_index == 1) || (func_index == 3)) {
            for (int i = 0; i < dim; i++) {
                b[i] = alpha[i];
            }
        }
        if (func_index == 6) {
            for (int i = 0; i < dim; i++) {
                beta[i] = 2 / M_PI;
            }
        }
#else
        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
            alpha[i] = 1.0;
            if (func_index == 3) {
                beta[i] = unif01(mt);
            } else if (func_index == 6) {
                beta[i] = 2 / M_PI;
            } else {
                beta[i] = 1.0;
            }
        }
#endif
        // RMSE　Root Mean Squared Error
        expected = genz_integral(func_index, dim, a, b, alpha, beta);
        if (rmse > 0) {
            Kahan esum;
            for (int z = 0; z < 100; z++) {
                Kahan sum;
                for (int i = 0; i < count; i++) {
                    const double *tuple = digitalNet.getPoint();
                    sum.add(genz_function(func_index, dim, tuple, alpha, beta));
                    digitalNet.nextPoint();
                }
                double er = expected - sum.get() / count;
                esum.add(er * er);
                digitalNet.setDigitalShift(true);
                digitalNet.pointInitialize();
            }
            return sqrt(esum.get() / 100);
        } else {
            Kahan sum;
            for (int i = 0; i < count; i++) {
                const double *tuple = digitalNet.getPoint();
                sum.add(genz_function(func_index, dim, tuple, alpha, beta));
                digitalNet.nextPoint();
            }
            return abs(expected - sum.get() / count);
        }
    }
}
