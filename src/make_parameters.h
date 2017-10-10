#pragma once
#ifndef MAKE_PARAMETERS_H
#define MAKE_PARAMETERS_H

#include <inttypes.h>

void printArray(std::string str, double array[], int len);
void makeParameter(int func_index, int dim, int seed, int original,
                   double a[], double b[], double alpha[], double beta[],
                   bool verbose, double difficulty = -1);

#endif // MAKE_PARAMETERS_H
