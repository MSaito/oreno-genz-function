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
//#include "kahan.hpp"
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
        int bits;
        int offset;
        int dn_id;
        char type;
        bool digital_shift;
        bool verbose;
        string dnfile;
    };
    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    template<typename D>
    void output(D& digitalNet, int count, const string& dn_name,
        bool digital_shift);
    int file_sai(cmd_opt_t& opt);
    int random_sai(cmd_opt_t& opt);
}

int main(int argc, char *argv[]) {
#if defined(DEBUG)
    cout << "main step 1" << endl;
#endif
    cmd_opt_t opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
#if defined(DEBUG)
    cout << "main step 2" << endl;
#endif
    if (opt.dn_id < 0) {
        return file_sai(opt);
    }
    if (opt.dn_id >= 100) {
        return random_sai(opt);
    }
#if defined(DEBUG)
    cout << "main step 4" << endl;
#endif
    DigitalNetID dnid = static_cast<DigitalNetID>(opt.dn_id);
    DigitalNet<uint64_t> dn(dnid, opt.s_dim, opt.start_m);
    int count = 1 << opt.start_m;
    output(dn, count, getDigitalNetName(opt.dn_id), opt.digital_shift);
    return 0;
}

namespace {
    int file_sai(cmd_opt_t& opt)
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
        int count = 1 << m;
        output(dn, count, opt.dnfile, opt.digital_shift);
        return 0;
    }

    int random_sai(cmd_opt_t& opt)
    {
        int s = opt.s_dim;
        RandomNet dn(s, 100);
        dn.pointInitialize();
        int mask = 64;
        dn.setMask(mask);
        int count = 1 << opt.start_m;
        output(dn, count, "Random", opt.digital_shift);
        return 0;
    }

    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -m start_m "
             << " [-d digitalnet_id] [-v] [-D] [digitalnet_file]" << endl;
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
            {"type", required_argument, NULL, 't'},
            {"bits", required_argument, NULL, 'b'},
            {"offset", required_argument, NULL, 'o'},
            {"digitalnet-id", required_argument, NULL, 'd'},
            {"digital-shift", no_argument, NULL, 'D'},
            {"verbose", no_argument, NULL, 'v'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.start_m = 0;
        opt.dn_id = -1;
        opt.type = ' ';
        opt.bits = 3;
        opt.offset = 0;
        opt.digital_shift = false;
        opt.verbose = false;
        errno = 0;
#if defined(DEBUG)
        cout << "parse_opt step 2" << endl;
#endif
        for (;;) {
            c = getopt_long(argc, argv, "s:m:t:b:o:d:vD", longopts, NULL);
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
            case 't':
                opt.type = optarg[0];
                break;
            case 'b':
                opt.bits = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "bits should be a number" << endl;
                    error = true;
                }
                break;
            case 'o':
                opt.offset = strtoul(optarg, NULL, 10);
                if (errno) {
                    cout << "offset should be a number" << endl;
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
            case 'v':
                opt.verbose = true;
                break;
            case 'D':
                opt.digital_shift = true;
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

    //単純に出力するだけ
    template<typename D>
    void output(D& digitalNet, int count, const string& dn_name,
                bool digital_shift)
    {
        int s = digitalNet.getS();
        cout << "# " << dn_name << " digital_shift = " << digital_shift << endl;
        cout << "# s = " << dec << s << endl;
        cout << "# count = " << dec << count << endl;
        cout << scientific << setprecision(18);
        digitalNet.setDigitalShift(digital_shift);
        digitalNet.pointInitialize();
        for (int i = 0; i < count; i++) {
            const double *tuple = digitalNet.getPoint();
            for (int j = 0; j < s; j++) {
                cout << tuple[j] << ",";
            }
            cout << endl;
            digitalNet.nextPoint();
        }
    }
}
