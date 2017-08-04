#include <MCQMCIntegration/DigitalNet.h>
#include <stdlib.h>
#include <errno.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>

#define DEBUG
#include "Genz.hpp"

using namespace std;
using namespace MCQMCIntegration;
using namespace GenzNS;

int main(int argc, char * argv[])
{
    if (argc < 7) {
        cout << argv[0] << " digital_net_id s start_m stop_m genz_id seed"
             << endl;
        return -1;
    }
    errno = 0;
    int digital_net_id = strtol(argv[1], NULL, 10);
    int s = strtol(argv[2], NULL, 10);
    int start_m = strtol(argv[3], NULL, 10);
    int stop_m = strtol(argv[4], NULL, 10);
    int genz_id = strtol(argv[5], NULL, 10);
    uint32_t seed = strtoull(argv[6], NULL, 10);
    if (errno) {
        cout << argv[0] << " digital_net_id s start_m stop_m genz_id seed"
             << endl;
        return -1;
    }
    cout << "#digital_net_id = " << dec << digital_net_id << endl;
    cout << "#digital_net_name = " << getDigitalNetName(digital_net_id)
    << endl;
    cout << "#s = " << dec << s << endl;
    cout << "#start_m = " << dec << start_m << endl;
    cout << "#stop_m = " << dec << stop_m << endl;
    cout << "#genz_id = " << dec << genz_id << endl;
    cout << "#seed = " << dec << seed << endl;
    mt19937_64 mt(seed);
    uniform_real_distribution<double> unif01(0.0, 1.0);
    double alpha[s];
    double beta[s];
    for (int i = 0; i < s; i++) {
        alpha[i] = unif01(mt);
        beta[i] = unif01(mt);
    }
    //int genz_size = 6;
    shared_ptr<GenzFunction> genz[] = {
        shared_ptr<GenzFunction>(new Oscillatory()),
        shared_ptr<GenzFunction>(new ProductPeak()),
        shared_ptr<GenzFunction>(new CornerPeak()),
        shared_ptr<GenzFunction>(new Gaussian()),
        shared_ptr<GenzFunction>(new C0Function()),
        shared_ptr<GenzFunction>(new Discontinuous())
    };
    if (genz_id < 1 || genz_id > 6) {
        cout << "genz_id should be 1 ... 6" << endl;
        return -1;
    }
    GenzFunction& func = *genz[genz_id-1];
    cout << "#genz_name = " << func.name() << endl;
    cout << "#m, abs err, log2(err)" << endl;
    for (int m = start_m; m <= stop_m; m++) {
        DigitalNetID dnid = static_cast<DigitalNetID>(digital_net_id);
        //cout << "dnid = " << dnid << endl;
        DigitalNet<uint64_t> dn(dnid, s, m);
        int count = 1 << m;
#if defined(DEBUG)
#endif
        double err = integral(func, dn, count, s, alpha, beta);
        double l = log2(err);
        cout << dec << m << "," << err << "," << l << endl;
    }
    return 0;
}
