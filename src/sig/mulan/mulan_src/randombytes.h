#ifndef MULAN_RANDOMBYTES_H
#define MULAN_RANDOMBYTES_H

#define _GNU_SOURCE

#include <unistd.h>

void mulan_randombytes(unsigned char *x, size_t xlen);

#endif
