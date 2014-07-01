#include <stdio.h>
#include <check.h>

/**
 * ccs, this is an example grabbed from here: http://isti.bitbucket.org/2012/06/01/cmake-check.html
 * it currently won't compile
 */

#include <gsl/gsl_histogram.h>
#include "flow-fns.h"



START_TEST (test_update)
{
  double testX = 1.0;
  fail_if(updateMean(0, 0.0, testX) != testX, "failed to update zero element mean");
  fail_if(updateMean(4, 2.5, testX) != 2.20, "updated mean is wrong");
}
END_TEST

START_TEST (test_sub)
{
  double testX = 1.0;
  fail_if(subMean(2, 1.0, 1.0) != 1.0, "failed to sub mean of single item");
  fail_if(subMean(5, 2.2, testX) != 2.5, "sub mean is wrong");
}
END_TEST


Suite* mean_suite (void) {
        Suite *suite = suite_create("Core");
        TCase *tcase = tcase_create("RunningMean");
        tcase_add_test(tcase, test_update);
        tcase_add_test(tcase, test_sub);        
        suite_add_tcase(suite, tcase);
        return suite;
}

int main (int argc, char *argv[]) {
        int number_failed;
        Suite *suite = mean_suite();
        SRunner *runner = srunner_create(suite);
        srunner_run_all(runner, CK_NORMAL);
        number_failed = srunner_ntests_failed(runner);
        srunner_free(runner);
        return number_failed;
}
