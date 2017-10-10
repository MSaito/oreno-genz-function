#define _USE_MATH_DEFINES
#include "config.h"
#include <math.h>
#include <cerrno>
#include <getopt.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "testpack.h"
#include "make_parameters.h"

#if defined(HAVE_MPI_H)
#include <mpi.h>
#include "change_output.h"
#else
#define MPI_Finalize()
#endif

using namespace std;

//#define DEBUG

namespace {
    struct cmd_opt_t {
        uint32_t s_dim;
        uint32_t e_dim;
        uint32_t add;
        uint32_t seed;
        int genz_no;
        int original;
    };
    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    int calc_theoretical(uint32_t dim, uint32_t seed, uint32_t genz_no,
                         int original);
}

int main(int argc, char *argv[]) {
    int rank;
    int num_process;
#if defined(HAVE_MPI_H)
    // MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
#else
    rank = 0;
    num_process = 1;
#endif
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        MPI_Finalize();
        return -1;
    }
#if defined(HAVE_MPI_H)
    if (!change_output(argv[0], rank, opt.s_dim, opt.genz_no, opt.seed)) {
        MPI_Finalize();
        return -1;
    }
#endif
    cout << "#" << genz_name(opt.genz_no) << endl;
    if (opt.original == 1) {
        cout << "#orignal genz parameters" << endl;
    } else if (opt.original == -1) {
        cout << "#saito old parameters" << endl;
    } else {
        cout << "#saito new parameters" << endl;
    }
    cout << "#func_index = " << opt.genz_no << endl;
    cout << "#seed = " << opt.seed << endl;
    cout << endl;
    for (uint32_t dim = opt.s_dim; dim <= opt.e_dim; dim += opt.add) {
        if (static_cast<int>(dim) % num_process == rank) {
            calc_theoretical(dim, opt.seed, opt.genz_no, opt.original);
        }
    }
    MPI_Finalize();
    return 0;
}

namespace {
    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -e e_dim -S seed -g genz_no -o" << endl;
        cout << "\t--s-dim, -s\t\tstart dim" << endl;
        cout << "\t--e-dim, -e\t\tend dim" << endl;
        cout << "\t--add, -a\t\tstep of dim" << endl;
        cout << "\t--seed, -S\t\tseed of random" << endl;
        cout << "\t--genz-no, -g\t\tgenz-no" << endl;
        cout << "\t--orignal, -o\t\torignal genz parameters" << endl;
    }

    bool parse_opt(cmd_opt_t& opt, int argc, char **argv)
    {
        int c;
        bool error = false;
        string pgm = argv[0];
        static struct option longopts[] = {
            {"s-dim", required_argument, NULL, 's'},
            {"e-dim", required_argument, NULL, 'e'},
            {"add", required_argument, NULL, 'a'},
            {"seed", required_argument, NULL, 'S'},
            {"genz-no", required_argument, NULL, 'g'},
            {"orignal", optional_argument, NULL, 'o'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.e_dim = 0;
        opt.add = 1;
        opt.seed = 1;
        opt.genz_no = 0;
        opt.original = 0;
        errno = 0;
        for (;;) {
            c = getopt_long(argc, argv, "s:e:a:S:g:o::", longopts, NULL);
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
            case 'e':
                opt.e_dim = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "e_dim should be a number" << endl;
                    error = true;
                }
                break;
            case 'a':
                opt.add = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "add should be a number" << endl;
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
            case 'o':
                if (optarg == NULL) {
                    opt.original = 1;
                } else if (optarg[0] == 'x') {
                    opt.original = -1;
                } else {
                    opt.original = strtol(optarg, NULL, 10);
                    if (errno) {
                        cout << "original shoud be one of {x(-1), 0, 1}."
                             << endl;
                        error = true;
                    }
                }
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
        if (opt.e_dim < opt.s_dim) {
            opt.e_dim = opt.s_dim;
        }
        return true;
    }

    int calc_theoretical(uint32_t dim, uint32_t seed, uint32_t genz_no,
                         int original)
    {
        double *a = new double[dim];
        double *b = new double[dim];
        double *alpha = new double[dim];
        double *beta = new double[dim];
        for (size_t i = 0; i < dim; i++) {
            a[i] = 0;
            b[i] = 0;
            alpha[i] = 0;
            beta[i] = 0;
        }
        cout << "#dim = " << dim << endl;
        makeParameter(genz_no, dim, seed, original,
                      a, b, alpha, beta, true);
        double value = genz_integral(genz_no, dim, a, b, alpha, beta);
        cout << "value = " << scientific << setprecision(18) << value << endl;
        cout << endl;
        delete[] a;
        delete[] b;
        delete[] alpha;
        delete[] beta;
        return 0;
    }
}
