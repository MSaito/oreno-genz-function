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
#include "saipack.hpp"
#include "kahan.hpp"
#include "RandomNet.hpp"
#include <memory>
#include <random>


using namespace std;
using namespace MCQMCIntegration;

//#define DEBUG

namespace {
    struct cmd_opt_t {
        uint32_t s_dim;
        uint32_t start_m;
        uint32_t end_m;
        uint32_t seed;
        int sai_no;
        int dn_id;
        int rmse;
        int parameter;
        bool verbose;
        string dnfile;
    };
#if 1
    shared_ptr<Saipack> functions[] = {
        shared_ptr<Saipack>(reinterpret_cast<Saipack *>(new AddSai())),
        shared_ptr<Saipack>(reinterpret_cast<Saipack *>(new MulSai())),
        shared_ptr<Saipack>(reinterpret_cast<Saipack *>(new SinSai())),
        shared_ptr<Saipack>(reinterpret_cast<Saipack *>(new PolSai()))
    };
#else
    Saipack functions[] = {
        AddSai(),
        MulSai(),
        SinSai(),
        PolSai()
    };
#endif
    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    template<typename D>
    double integral(Saipack& sai, D& digitalNet, int count,
                    double expected, int rmse,
                    bool verbose);
    int file_sai(cmd_opt_t& opt, Saipack& sai, double expected);
    int random_sai(cmd_opt_t& opt, Saipack& sai, double expected);
//    template<typename D>
//    void loop_integral(D& digitalNet, cmd_opt_t& opt,
//                       int s, int start_m, int end_m);
    void print_header(cmd_opt_t& opt, const string& test_name,
                      const string& dn_name,
                      double expected);
}

int main(int argc, char *argv[]) {
#if defined(DEBUG)
    cout << "main step 1" << endl;
#endif
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    double a[opt.s_dim];
    double b[opt.s_dim];
    for (size_t i = 0; i < opt.s_dim; i++) {
        a[i] = 0;
        b[i] = 0;
    }
#if defined(DEBUG)
    cout << "main step 2" << endl;
#endif
    Saipack& func = *functions[opt.sai_no];
    std::mt19937_64 mt(opt.seed);
    func.makeParameter(opt.parameter, opt.s_dim, mt, a, b, opt.verbose);
    func.setParam(opt.s_dim, a, b);
    double expected = func.expected(opt.s_dim, a, b);
#if defined(DEBUG)
    cout << "main step 3" << endl;
#endif
    if (opt.dn_id < 0) {
        return file_sai(opt, func, expected);
    }
    if (opt.dn_id >= 100) {
        return random_sai(opt, func, expected);
    }
#if defined(DEBUG)
    cout << "main step 4" << endl;
#endif
    DigitalNetID dnid = static_cast<DigitalNetID>(opt.dn_id);
    print_header(opt, func.getName(), getDigitalNetName(opt.dn_id),
                 expected);
    for (uint32_t m = opt.start_m; m <= opt.end_m; m++) {
        DigitalNet<uint64_t> dn(dnid, opt.s_dim, m);
        int count = 1 << m;
        double error = integral(func, dn, count,
                                expected, opt.rmse, opt.verbose);
        cout << dec << m << "," << error << "," << log2(error) << endl;
    }
    return 0;
}

namespace {
    void print_header(cmd_opt_t& opt, const string& test_name,
                      const string& dn_name,
                      double expected)
    {
    cout << "#" << test_name << endl;
    cout << "#" << dn_name << endl;
    if (opt.rmse > 0) {
        cout << "#m, abs err, log2(RMSE[" << dec << opt.rmse << "])" << endl;
    } else {
        cout << "#m, abs err, log2(err)" << endl;
    }
    cout << "#expected = " << expected << endl;
    }

    int file_sai(cmd_opt_t& opt, Saipack& func, double expected)
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
        if (dn.getS() != opt.s_dim) {
            cout << "s_dim != dn.getS()" << endl;
            return -1;
        }
        int m = dn.getM();
        print_header(opt, func.getName(), opt.dnfile, expected);
        cout << "# m = " << dec << m << endl;
        int count = 1 << m;
        double error = integral(func, dn, count,
                                expected, opt.rmse, opt.verbose);
        cout << dec << m << "," << error << "," << log2(error) << endl;
        return 0;
    }

    int random_sai(cmd_opt_t& opt, Saipack& func, double expected)
    {
        int s = opt.s_dim;
        RandomNet dn(s, 100);
        dn.pointInitialize();
        print_header(opt, func.getName(), "Random", expected);
        int mask = 64;
        dn.setMask(mask);
        for (uint32_t m = opt.start_m; m <= opt.end_m; m++) {
            int count = 1 << m;
            double error = integral(func, dn, count,
                                    expected, opt.rmse,
                                    opt.verbose);
            cout << dec << m << "," << error << "," << log2(error) << endl;
        }
        return 0;
    }

    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -m start_m -M end_m -S seed -n sai_no"
             << " [-d digitalnet_id] [-p] [-v] [digitalnet_file]" << endl;
    }

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv)
    {
#if defined(DEBUG)
        cout << "parse_opt step 1" << endl;
#endif
        int c;
        bool error = false;
        string pgm = argv[0];
        static struct option longopts[] = {
            {"s-dim", required_argument, NULL, 's'},
            {"start-m", required_argument, NULL, 'm'},
            {"end-m", required_argument, NULL, 'M'},
            {"seed", required_argument, NULL, 'S'},
            {"sai-no", required_argument, NULL, 'n'},
            {"digitalnet-id", required_argument, NULL, 'd'},
            {"rmse", optional_argument, NULL, 'r'},
            {"parameter", required_argument, NULL, 'p'},
            {"verbose", no_argument, NULL, 'v'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.start_m = 0;
        opt.end_m = 0;
        opt.seed = 1;
        opt.sai_no = 0;
        opt.dn_id = -1;
        opt.rmse = 0;
        opt.parameter = 0;
        opt.verbose = false;
        errno = 0;
#if defined(DEBUG)
        cout << "parse_opt step 2" << endl;
#endif
        for (;;) {
            c = getopt_long(argc, argv, "s:m:M:S:n:d:r:p:v", longopts, NULL);
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
            case 'n':
                opt.sai_no = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "sai_no should be a number" << endl;
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
            case 'p':
                opt.parameter = strtol(optarg, NULL, 10);
                if (errno) {
                    cout << "param shoud be a number." << endl;
                    error = true;
                }
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
    double integral(Saipack& func, D& digitalNet, int count,
                    double expected, int rmse,
                    bool verbose)
    {
#if defined(DEBUG)
        cout << "name = " << func.getName() << endl;
        cout << "count = " << dec << count << endl;
        cout << "rmse = " << dec << rmse << endl;
#endif
        // RMSE　Root Mean Squared Error
        if (rmse > 0) {
            Kahan esum;
            for (int z = 0; z < 100; z++) {
                Kahan sum;
                for (int i = 0; i < count; i++) {
                    const double *tuple = digitalNet.getPoint();
                    sum.add(func(tuple));
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
                sum.add(func(tuple));
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
