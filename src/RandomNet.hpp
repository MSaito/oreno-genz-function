#pragma once
#ifndef RANDOMNET_HPP
#define RANDOMNET_HPP

#include "mt19937_64.hpp"

namespace MCQMCIntegration {

    class RandomNet {
    public:
        RandomNet(int s, uint32_t seed) {
            this->s = s;
            mt.seed(seed);
            mask = 0;
            mask = ~mask;
            point = new double[s];
        }
        ~RandomNet() {
            delete[] point;
        }
        void setMask(int m) {
            mask = 0;
            mask = ~mask >> (64 - m);
            mask = mask << (64 - m);
        }
        const double * getPoint() const {
            return point;
        }
        void pointInitialize() {
            nextPoint();
        }
        void nextPoint() {
            for (int i = 0; i < s; i++) {
                point[i] = ((mt.getUint64() & mask) >> 11)
                    * (1.0/9007199254740992.0)
                    + pow(2.0, -54);
            }
        }
        void setDigitalShift(bool) {
        }
        int getS() {
            return s;
        }
    private:
        int s;
        uint64_t mask;
        double * point;
        mt19937_64 mt;
    };
}

#endif // RANDOMNET_HPP
