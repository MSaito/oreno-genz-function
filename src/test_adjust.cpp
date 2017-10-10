#include <iostream>
#include <iomanip>
#include <cerrno>
#include <getopt.h>
#include <cstdlib>
#include <string>
#include "adjust_parameters.h"
#include "make_parameters.h"

using namespace std;

namespace {
    struct cmd_opt_t {
        uint32_t s_dim;
        int genz_no;
        bool verbose;
    };

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
}

int main(int argc, char * argv[])
{
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    int dim = opt.s_dim;
    int genz_no = opt.genz_no;
    int seed = 1;
    int original = 1;
    double a[dim];
    double b[dim];
    double alpha[dim];
    double beta[dim];
    for (int i = 0; i < dim; i++) {
        a[i] = 0;
        b[i] = 0;
        alpha[i] = 0;
        beta[i] = 0;
    }
    cout << "#dim = " << dim << endl;
    makeParameter(genz_no, dim, seed, original, a, b, alpha, beta, false);
    double expect = adjustParameter(genz_no, dim, a, b,
                                    alpha, beta, opt.verbose);
    cout << "expect = " << dec << expect << endl;
    return 0;
}

namespace {
    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -g genz_no -v" << endl;
        cout << "\t--s-dim, -s\t\tdim" << endl;
        cout << "\t--genz-no, -g\t\tgenz-no" << endl;
        cout << "\t--verbosel, -v\t\tverbose" << endl;
    }

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv)
    {
        int c;
        bool error = false;
        string pgm = argv[0];
        static struct option longopts[] = {
            {"s-dim", required_argument, NULL, 's'},
            {"genz-no", required_argument, NULL, 'g'},
            {"verbose", no_argument, NULL, 'v'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.genz_no = 0;
        opt.verbose = 0;
        errno = 0;
        for (;;) {
            c = getopt_long(argc, argv, "s:g:v", longopts, NULL);
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
            case 'g':
                opt.genz_no = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "genz_no should be a number" << endl;
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
        if (opt.genz_no < 1 || opt.genz_no > 6) {
            cout << "genz_no shoule be 1 <= genz_no <= 6" << endl;
            error = true;
        }
        if (error) {
            cmd_message(pgm);
            return false;
        }
        return true;
    }
}
