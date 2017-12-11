#pragma once
#ifndef MAKE_WAFOMC_PARAMETERS_H
#define MAKE_WAFOMC_PARAMETERS_H

#include <inttypes.h>

void makeWafomParameter(int func_index, int dim, int seed,
                        double a[], double b[], double alpha[], double beta[],
                        bool verbose, double mag);

#endif // MAKE_WAFOMC_PARAMETERS_H
