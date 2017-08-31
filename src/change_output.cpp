#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <cerrno>
#include <cstdio>

#include "change_output.h"

using namespace std;

bool change_output(const char * pgm, int rank,
                   int s, int g, int seed)
{
    char fname[500];
    errno = 0;
    sprintf(fname, "%s.rank%d.s%d.g%d.S%d.txt", pgm, rank, s, g, seed);
    int fd;
    fd = open(fname, O_WRONLY | O_CREAT, 0660);
    if (fd < 0) {
        cerr << pgm << ": open file error\n" << endl;
        return false;
    }
    dup2(fd, 1);
    if (errno) {
        cerr << pgm << ": dup error.\n" << endl;
        close(fd);
        return false;
    }
    return true;
}
