

/* This is auto-generated code. Edit at your own peril. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "CuTest.h"


extern void test_side_length(CuTest*);
extern void test_make_map(CuTest*);
extern void test_healpix_distance_to_radec(CuTest*);
extern void test_healpix_neighbours(CuTest*);
extern void test_big_nside(CuTest*);
extern void test_distortion_at_pole(CuTest*);


void RunAllTests(void) 
{
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_side_length);
    SUITE_ADD_TEST(suite, test_make_map);
    SUITE_ADD_TEST(suite, test_healpix_distance_to_radec);
    SUITE_ADD_TEST(suite, test_healpix_neighbours);
    SUITE_ADD_TEST(suite, test_big_nside);
    SUITE_ADD_TEST(suite, test_distortion_at_pole);

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);

    printf("%s\n", output->buffer);

}

int main(int argc, char** args)
{
    RunAllTests();
}
