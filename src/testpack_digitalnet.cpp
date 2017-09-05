#define _USE_MATH_DEFINES
#include <math.h>
#include <cerrno>
#include <getopt.h>
#include <cstdlib>
#include <string>
#include <MCQMCIntegration/DigitalNet.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "testpack.h"
#include "kahan.hpp"
#include "make_parameters.hpp"
#include "RandomNet.hpp"

using namespace std;
using namespace MCQMCIntegration;

//#define DEBUG

namespace {
    struct cmd_opt_t {
        uint32_t s_dim;
        uint32_t start_m;
        uint32_t end_m;
        uint32_t seed;
        int genz_no;
        int dn_id;
        int rmse;
        int original;
        bool verbose;
        string dnfile;
    };

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    template<typename D>
    double integral(int func_index, D& digitalNet, int count, int dim,
                    double alpha[], double beta[], double expected, int rmse,
                    bool verbose);
    int file_genz(cmd_opt_t& opt);
    int random_genz(cmd_opt_t& opt);
}

int main(int argc, char *argv[]) {
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    if (opt.dn_id < 0) {
        return file_genz(opt);
    }
    if (opt.dn_id >= 100) {
        return random_genz(opt);
    }
    double a[opt.s_dim];
    double b[opt.s_dim];
    double alpha[opt.s_dim];
    double beta[opt.s_dim];
    for (size_t i = 0; i < opt.s_dim; i++) {
        a[i] = 0;
        b[i] = 0;
        alpha[i] = 0;
        beta[i] = 0;
    }
    makeParameter(opt.genz_no, opt.s_dim, opt.seed, opt.original,
                  a, b, alpha, beta, opt.verbose);
    double expected = genz_integral(opt.genz_no, opt.s_dim, a, b, alpha, beta);
    DigitalNetID dnid = static_cast<DigitalNetID>(opt.dn_id);
    cout << "#" << genz_name(opt.genz_no) << endl;
    cout << "#" << getDigitalNetName(opt.dn_id) << endl;
    if (opt.rmse > 0) {
        cout << "#m, abs err, log2(RMSE[" << dec << opt.rmse << "])" << endl;
    } else {
        cout << "#m, abs err, log2(err)" << endl;
    }
    cout << "#expected = " << expected << endl;
    for (uint32_t m = opt.start_m; m <= opt.end_m; m++) {
        DigitalNet<uint64_t> dn(dnid, opt.s_dim, m);
        int count = 1 << m;
        double error = integral(opt.genz_no, dn, count, opt.s_dim,
                                alpha, beta, expected, opt.rmse, opt.verbose);
        cout << dec << m << "," << error << "," << log2(error) << endl;
    }
    return 0;
}

namespace {
    int file_genz(cmd_opt_t& opt)
    {
        ifstream dnstream(opt.dnfile);
        if (!dnstream) {
            cout << "can't open digital_net_file" << endl;
            return -1;
        }
        DigitalNet<uint64_t> dn(dnstream);
        dn.pointInitialize();
#if defined(DEBUG) && 0
        dn.showStatus(cout);
#endif
        int s = dn.getS();
        int m = dn.getM();
        double a[s];
        double b[s];
        double alpha[s];
        double beta[s];
        for (int i = 0; i < s; i++) {
            a[i] = 0;
            b[i] = 0;
            alpha[i] = 0;
            beta[i] = 0;
        }
        makeParameter(opt.genz_no, s, opt.seed, opt.original,
                      a, b, alpha, beta, opt.verbose);
        double expected = genz_integral(opt.genz_no, s, a, b, alpha, beta);
        cout << "#" << genz_name(opt.genz_no) << endl;
        cout << "# filename = " << opt.dnfile << endl;
        cout << "# s = " << dec << s << endl;
        cout << "# m = " << dec << m << endl;
        if (opt.rmse > 0) {
            cout << "#m, abs err, log2(RMSE[" << dec << opt.rmse << "])"
                 << endl;
        } else {
            cout << "#m, abs err, log2(err)" << endl;
        }
        int count = 1 << m;
        double error = integral(opt.genz_no, dn, count, s,
                                alpha, beta, expected, opt.rmse, opt.verbose);
        cout << dec << m << "," << error << "," << log2(error) << endl;
        return 0;
    }

    int random_genz(cmd_opt_t& opt)
    {
        int s = opt.s_dim;
        RandomNet dn(s, 100);
        dn.pointInitialize();
        double a[s];
        double b[s];
        double alpha[s];
        double beta[s];
        for (int i = 0; i < s; i++) {
            a[i] = 0;
            b[i] = 0;
            alpha[i] = 0;
            beta[i] = 0;
        }

        makeParameter(opt.genz_no, s, opt.seed, opt.original,
                      a, b, alpha, beta, opt.verbose);
        double expected = genz_integral(opt.genz_no, s, a, b, alpha, beta);
        cout << "#" << genz_name(opt.genz_no) << endl;
        cout << "# Random" << endl;
        cout << "# s = " << dec << s << endl;
        if (opt.dn_id >= 101) {
            cout << "# mask by m" << endl;
        } else {
            cout << "# not mask by m" << endl;
        }
        if (opt.rmse > 0) {
            cout << "#m, abs err, log2(RMSE[" << dec << opt.rmse << "])" << endl;
        } else {
            cout << "#m, abs err, log2(err)" << endl;
        }
        cout << "#expected = " << expected << endl;
        int mask = 64;
        dn.setMask(mask);
        for (uint32_t m = opt.start_m; m <= opt.end_m; m++) {
            if (opt.dn_id >= 101) {
                dn.setMask(m);
            }
            int count = 1 << m;
            double error = integral(opt.genz_no, dn, count, opt.s_dim,
                                    alpha, beta, expected, opt.rmse,
                                    opt.verbose);
            cout << dec << m << "," << error << "," << log2(error) << endl;
        }
    return 0;
    }

    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -m start_m -M end_m -S seed -g genz_no"
             << " [-d digitalnet_id] [digitalnet_file] [-o] [-v]" << endl;
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
            {"orignal", optional_argument, NULL, 'o'},
            {"bad-parameter", no_argument, NULL, 'x'},
            {"verbose", no_argument, NULL, 'v'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.start_m = 0;
        opt.end_m = 0;
        opt.seed = 1;
        opt.genz_no = 0;
        opt.dn_id = -1;
        opt.rmse = 0;
        opt.original = 0;
        opt.verbose = false;
        errno = 0;
        for (;;) {
            c = getopt_long(argc, argv, "s:m:M:S:g:d:r:o::vx", longopts, NULL);
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
            case 'o':
                if (optarg == NULL) {
                    opt.original = 1;
                } else if (optarg[0] == 'x') {
                    opt.original = -1;
                } else {
                    opt.original = strtol(optarg, NULL, 10);
                    if (errno) {
                        cout << "original shoud be one of {-1, 0, 1}." << endl;
                        error = true;
                    }
                }
                break;
            case 'x':
                opt.original = -1;
                break;
            case 'v':
                opt.verbose = true;
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
        if (opt.dn_id >= 0) {
            return true;
        }
        argc -= optind;
        argv += optind;
        if (argc <= 0) {
            error = true;
        } else {
            opt.dnfile = argv[0];
        }
        if (error) {
            cmd_message(pgm);
            return false;
        }
        return true;
    }

    template<typename D>
    double integral(int func_index, D& digitalNet, int count, int dim,
                    double alpha[], double beta[], double expected, int rmse,
                    bool verbose)
    {
#if defined(DEBUG)
        cout << "func_index = " << dec << func_index << endl;
        cout << "count = " << dec << count << endl;
        cout << "dim = " << dec << dim << endl;
        cout << "rmse = " << dec << rmse << endl;
#endif
        // RMSE　Root Mean Squared Error
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
#if defined(DEBUG)
                cout << "expected = " << expected << endl;
                cout << "calculated = " << (sum.get() / count) << endl;
#endif
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
#if defined(DEBUG)
            cout << "expected = " << expected << endl;
#endif
            if (verbose) {
                cout << "calculated = " << (sum.get() / count) << endl;
            }
            return abs(expected - sum.get() / count);
        }
    }
}
