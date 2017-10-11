noinst_bindir=./
genz_files = Genz.hpp
testpack_file = testpack.h kahan.hpp make_parameters.h mt19937_64.hpp \
change_output.h config.h RandomNet.hpp adjust_parameters.h cvmean.h

noinst_bin_PROGRAMS = genz_test testpack_digitalnet calc_theoretical \
saipack_digitalnet count_highbit simpleout test_adjust \
test_calc_cvmean

genz_test_SOURCES = genz_test.cpp Genz.cpp $(genz_files)
testpack_digitalnet_SOURCES = testpack_digitalnet.cpp testpack.cpp \
make_parameters.cpp adjust_parameters.cpp $(testpack_files)
calc_theoretical_SOURCES = calc_theoretical.cpp testpack.cpp \
make_parameters.cpp $(testpack_files)
saipack_digitalnet_SOURCES = saipack_digitalnet.cpp saipackpack.hpp \
make_parameters.cpp $(testpack_files)
count_highbit_SOURCES = count_highbit.cpp $(testpack_files)
simpleout_SOURCES = simpleout.cpp $(testpack_files)
test_adjust_SOURCES = test_adjust.cpp adjust_parameters.cpp testpack.cpp \
make_parameters.cpp $(testpack_files)
test_calc_cvmean_SOURCES = test_calc_cvmean.cpp cvmean.cpp $(testpack_files)

AM_CXXFLAGS = -I../src -O3 -Wall -Wextra -D__STDC_CONSTANT_MACROS
