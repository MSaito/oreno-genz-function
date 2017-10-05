noinst_bindir=./
genz_files = Genz.hpp
testpack_file = testpack.h kahan.hpp make_parameters.hpp mt19937_64.hpp \
change_output.h config.h RandomNet.hpp

noinst_bin_PROGRAMS = genz_test testpack_digitalnet calc_theoretical \
saipack_digitalnet count_highbit simpleout

genz_test_SOURCES = genz_test.cpp Genz.cpp $(genz_files)
testpack_digitalnet_SOURCES = testpack_digitalnet.cpp testpack.cpp \
$(testpack_files)
calc_theoretical_SOURCES = calc_theoretical.cpp testpack.cpp \
$(testpack_files)
saipack_digitalnet_SOURCES = saipack_digitalnet.cpp saipackpack.hpp \
$(testpack_files)
count_highbit_SOURCES = count_highbit.cpp $(testpack_files)
simpleout_SOURCES = simpleout.cpp $(testpack_files)

AM_CXXFLAGS = -I../src -O3 -Wall -Wextra -D__STDC_CONSTANT_MACROS
