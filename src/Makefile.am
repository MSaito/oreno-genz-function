noinst_bindir=./
genz_files = Genz.hpp
testpack_file = testpack.h kahan.hpp

noinst_bin_PROGRAMS = genz_test testpack_digitalnet

genz_test_SOURCES = genz_test.cpp Genz.cpp $(genz_files)
testpack_digitalnet_SOURCES = testpack_digitalnet.cpp testpack.cpp \
$(testpack_files)

AM_CXXFLAGS = -I../src -O3 -Wall -Wextra -D__STDC_CONSTANT_MACROS
