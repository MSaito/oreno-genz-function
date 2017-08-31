#pragma once
#ifndef MT19937_64_HPP
#define MT19937_64_HPP
/*
  C++ class by Mutsuo Saito. 2017-8-23.
  This is C++ code for machines which does not support C++11.
  bug report: saito@manieth.com
*/

/*
  A C-program for MT19937-64 (2014/2/23 version).
  Coded by Takuji Nishimura and Makoto Matsumoto.

  This is a 64-bit version of Mersenne Twister pseudorandom number
  generator.

  Before using, initialize the state by using init_genrand64(seed)
  or init_by_array64(init_key, key_length).

  Copyright (C) 2004, 2014, Makoto Matsumoto and Takuji Nishimura,
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  3. The names of its contributors may not be used to endorse or promote
  products derived from this software without specific prior written
  permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  References:
  T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
  ACM Transactions on Modeling and
  Computer Simulation 10. (2000) 348--357.
  M. Matsumoto and T. Nishimura,
  ``Mersenne Twister: a 623-dimensionally equidistributed
  uniform pseudorandom number generator''
  ACM Transactions on Modeling and
  Computer Simulation 8. (Jan. 1998) 3--30.

  Any feedback is very welcome.
  http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

#include <inttypes.h>

class mt19937_64 {
public:
    mt19937_64() {
        MATRIX_A = UINT64_C(0xB5026F5AA96619E9);
        UM = UINT64_C(0xFFFFFFFF80000000);
        LM = UINT64_C(0x7FFFFFFF);
        array = new uint64_t[NN];
        seed(UINT64_C(19650218));
    }

    mt19937_64(uint64_t seed_value) {
        MATRIX_A = UINT64_C(0xB5026F5AA96619E9);
        UM = UINT64_C(0xFFFFFFFF80000000);
        LM = UINT64_C(0x7FFFFFFF);
        array = new uint64_t[NN];
        seed(seed_value);
    }

    ~mt19937_64() {
        delete[] array;
    }

    void seed(uint64_t value) {
        array[0] = value;
        for (int i = 1; i < NN; i++) {
            array[i] =  (UINT64_C(6364136223846793005)
                         * (array[i - 1] ^ (array[i - 1] >> 62)) + i);
        }
        index = NN;
    }

    void seed(uint64_t init_key[], uint64_t key_length) {
        unsigned int i, j;
        uint64_t k;
        seed(UINT64_C(19650218));
        i=1;
        j=0;
        k = (NN>key_length ? NN : key_length);
        for (; k; k--) {
            array[i] = (array[i] ^ ((array[i-1] ^ (array[i-1] >> 62))
                                    * UINT64_C(3935559000370003845)))
                + init_key[j] + j; /* non linear */
            i++; j++;
            if (i>=NN) { array[0] = array[NN-1]; i=1; }
            if (j>=key_length) j=0;
        }
        for (k=NN-1; k; k--) {
            array[i] = (array[i] ^ ((array[i-1] ^ (array[i-1] >> 62))
                                    * UINT64_C(2862933555777941757)))
                - i; /* non linear */
            i++;
            if (i>=NN) { array[0] = array[NN-1]; i=1; }
        }
        /* MSB is 1; assuring non-zero initial array */
        array[0] = UINT64_C(1) << 63;
        index = NN;
    }

    uint64_t getUint64() {
        if (index >= NN) {
            genall();
        }
        uint64_t x = array[index++];
        return temper(x);
    }

    uint32_t getUint32() {
        return getUint64() >> 32;
    }

    /* generates a random number on [0,1)-real-interval */
    double getDouble01() {
        return (getUint64() >> 11) * (1.0/9007199254740992.0);
    }

    /* generates a random number on [0,1]-real-interval */
    double getDouble01Close() {
        return (getUint64() >> 11) * (1.0/9007199254740991.0);
    }

    /* generates a random number on (0,1)-real-interval */
    double getDouble01Open() {
        return ((getUint64() >> 12) + 0.5) * (1.0/4503599627370496.0);
    }

private:
    enum {NN = 312, MM = 156};
    int index;
    uint64_t *array;
    uint64_t MATRIX_A;
    /* Most significant 33 bits */
    uint64_t UM;
    /* Least significant 31 bits */
    uint64_t LM;

    void genall() {
        int i;
        uint64_t x;
        const uint64_t mag01[2]={UINT64_C(0), MATRIX_A};
        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        for (i=0;i<NN-MM;i++) {
            x = (array[i]&UM)|(array[i+1]&LM);
            array[i] = array[i+MM] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        for (;i<NN-1;i++) {
            x = (array[i]&UM)|(array[i+1]&LM);
            array[i] = array[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        x = (array[NN-1]&UM)|(array[0]&LM);
        array[NN-1] = array[MM-1] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        index = 0;
    }

    uint64_t temper(uint64_t y) {
        y ^= (y >> 29) & UINT64_C(0x5555555555555555);
        y ^= (y << 17) & UINT64_C(0x71D67FFFEDA60000);
        y ^= (y << 37) & UINT64_C(0xFFF7EEE000000000);
        y ^= (y >> 43);
        return y;
    }
};

#endif // mt19937_64
