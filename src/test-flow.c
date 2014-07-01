#include <stdio.h>
#include <check.h>

#include "isti.h"
#include "isti_str.h"


START_TEST (test_str)
{
        isti_str *s = NULL;
        fail_if(isti_str_alloc(&s), "Could not allocate s");
        fail_if(!s, "s null after allocation");
        fail_if(s->free(&s, ISTI_OK), "Could not free s");
        fail_if(s, "s not null after free");
}
END_TEST

Suite* str_suite (void) {
        Suite *suite = suite_create("isti_str");
        TCase *tcase = tcase_create("case");
        tcase_add_test(tcase, test_str);
        suite_add_tcase(suite, tcase);
        return suite;
}

int main (int argc, char *argv[]) {
        int number_failed;
        Suite *suite = str_suite();
        SRunner *runner = srunner_create(suite);
        srunner_run_all(runner, CK_NORMAL);
        number_failed = srunner_ntests_failed(runner);
        srunner_free(runner);
        return number_failed;
}
