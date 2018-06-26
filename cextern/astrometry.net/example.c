// Simple hello world example
#include <stdio.h>
#include "healpix.h"

int main(int argc, char const *argv[])
{
    double side_length = healpix_side_length_arcmin(4);
    printf("side_length: %f\n", side_length);
    return 0;
}
