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
        int bits;
        int offset;
        int dn_id;
        char type;
        bool verbose;
        string dnfile;
    };
    bool parse_opt(cmd_opt_t& opt, int argc, char **argv);
    void cmd_message(const string& pgm);
    template<typename D>
    void counter1(D& digitalNet, int m, int count, const string& dn_name);
    template<typename D>
    void counterh(D& digitalNet, int count, int bits, int offset,
                  const string& dn_name);
    template<typename D>
    void counterv(D& digitalNet, int m, int count, int bits, int offset,
                  const string& dn_name);
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
    if (opt.type == '1') {
        counter1(dn, opt.start_m, count, getDigitalNetName(opt.dn_id));
    } else if (opt.type == 'v') {
        counterv(dn, opt.start_m, count, opt.bits, opt.offset,
                 getDigitalNetName(opt.dn_id));
    } else {
        counterh(dn, count, opt.bits, opt.offset,
                 getDigitalNetName(opt.dn_id));
    }
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
        if (opt.type == '1') {
            counter1(dn, m, count, opt.dnfile);
        } else if (opt.type == 'v') {
            counterv(dn, m, count, opt.bits, opt.offset, opt.dnfile);
        } else {
            counterh(dn, count, opt.bits, opt.offset, opt.dnfile);
        }
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
        if (opt.type == '1') {
            counter1(dn, opt.start_m, count, "Random");
        } else if (opt.type == 'v') {
            counterv(dn, opt.start_m, count, opt.bits, opt.offset, "Random");
        } else {
            counterh(dn, count, opt.bits, opt.offset, "Random");
        }
        return 0;
    }

    void cmd_message(const string& pgm)
    {
        cout << pgm << " -s s_dim -m start_m -t type [-b bits] [-o offset]"
             << " [-d digitalnet_id] [-v] [digitalnet_file]" << endl;
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
            {"verbose", no_argument, NULL, 'v'},
            {NULL, 0, NULL, 0}};
        opt.s_dim = 0;
        opt.start_m = 0;
        opt.dn_id = -1;
        opt.type = ' ';
        opt.bits = 3;
        opt.offset = 0;
        opt.verbose = false;
        errno = 0;
#if defined(DEBUG)
        cout << "parse_opt step 2" << endl;
#endif
        for (;;) {
            c = getopt_long(argc, argv, "s:m:t:b:o:d:v", longopts, NULL);
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

    //ビットの数を数える
    template<typename D>
    void counter1(D& digitalNet, int m, int count, const string& dn_name)
    {
#if defined(DEBUG)
        cout << "count = " << dec << count << endl;
#endif
        int s = digitalNet.getS();
        cout << "# " << dn_name << endl;
        cout << "# s = " << dec << s << endl;
        cout << "# count = " << dec << count << endl;
        int sum[s][m];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < m; j++) {
                sum[i][j] = 0;
            }
        }
        digitalNet.pointInitialize();
        for (int i = 0; i < count; i++) {
            const double *tuple = digitalNet.getPoint();
            for (int j = 0; j < s; j++) {
                double mag = 2;
                for (int k = 0; k < m; k++) {
                    if (static_cast<uint32_t>(tuple[j] * mag) & 1) {
                        sum[j][k] += 1;
                    }
                    mag = mag * 2;
                }
            }
            digitalNet.nextPoint();
        }
        cout << "#s, ";
        for (int i = 0; i < m; i++) {
            cout << dec << i << ", ";
        }
        cout << endl;
        for (int i = 0; i < s; i++) {
            cout << dec << i << ", ";
            for (int j = 0; j < m; j++) {
                cout << dec << sum[i][j] << ",";
            }
            cout << endl;
        }
    }
    // 横方向にカウントする
    template<typename D>
    void counterh(D& digitalNet, int count, int bits, int offset,
                  const string& dn_name)
    {
#if defined(DEBUG)
        cout << "count = " << dec << count << endl;
#endif
        int s = digitalNet.getS();
        cout << "# " << dn_name << endl;
        cout << "# s = " << dec << s << endl;
        cout << "# offset = " << dec << offset << endl;
        cout << "# count = " << dec << count << endl;
        int k = 1 << bits;
        int o = 1 << offset;
        int omask = k - 1;
        int sum[s][k];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < k; j++) {
                sum[i][j] = 0;
            }
        }
        digitalNet.pointInitialize();
        for (int i = 0; i < count; i++) {
            const double *tuple = digitalNet.getPoint();
            for (int j = 0; j < s; j++) {
                int mask = static_cast<int>(tuple[j] * k * o) & omask;
                sum[j][mask] += 1;
            }
            digitalNet.nextPoint();
        }
        cout << "#s, ";
        for (int i = 0; i < k; i++) {
            cout << dec << i << ", ";
        }
        cout << endl;
        for (int i = 0; i < s; i++) {
            cout << dec << i << ", ";
            for (int j = 0; j < k; j++) {
                cout << dec << sum[i][j] << ",";
            }
            cout << endl;
        }
    }

    // 縦方向にカウントする（全部は無理）
    template<typename D>
    void counterv(D& digitalNet, int m, int count, int bits, int offset,
                  const string& dn_name)
    {
#if defined(DEBUG)
        cout << "count = " << dec << count << endl;
#endif
        int s = digitalNet.getS();
        cout << "# counter_v " << dn_name << endl;
        cout << "# s = " << dec << s << endl;
        cout << "# offset = " << dec << offset << endl;
        cout << "# m = " << dec << m << endl;
        cout << "# count = " << dec << count << endl;
        int c = 1 << bits;
        uint64_t sum[m][c];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < c; j++) {
                sum[i][j] = 0;
            }
        }
        digitalNet.pointInitialize();
#if defined(DEBUG)
        cout << "counter v step 1 c = " << dec << c << endl;
#endif
        // このへんからよく考える
        for (int i = 0; i < count; i++) {
            const double *tuple = digitalNet.getPoint();
            double mag = 2;
            for (int j = 0; j < m; j++) {
#if defined(DEBUG)
                cout << "j, m = " << dec << j << "," << m << endl;
#endif
                int kei = 0;
                //for (int k = 0; k < s && k < bits; k++) {
                for (int k = offset; k < offset + bits && k < s; k++) {
                    kei = kei << 1;
                    if (static_cast<int>(tuple[k] * mag) & 1) {
                        kei += 1;
                    }
                }
#if defined(DEBUG)
                cout << "j, kei = " << dec << j << "," << kei << endl;
#endif
                sum[j][kei] += 1;
                mag = mag * 2;
            }
#if defined(DEBUG)
            cout << "before nextPoint" << endl;
#endif
            digitalNet.nextPoint();
        }
#if defined(DEBUG)
        cout << "counter v step 2 " << endl;
#endif
        cout << "#m, ";
        for (int i = 0; i < c; i++) {
            cout << dec << i << ", ";
        }
        cout << endl;
        for (int i = 0; i < m; i++) {
            cout << dec << i << ", ";
            for (int j = 0; j < c; j++) {
                cout << dec << sum[i][j] << ",";
            }
            cout << endl;
        }
    }
}
