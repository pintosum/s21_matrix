#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "s21_matrix.h"

#define RED "\033[31m"
#define GREEN "\033[32m"
#define NOCOLOR "\033[0m"

int s21_equal_double(double, double);
int s21_is_valid_matrix_t(matrix_t *a);

void s21_fill_matrix(matrix_t *a, const char *numbers) {
  char *copy = calloc(strlen(numbers), 1);
  strncpy(copy, numbers, strlen(numbers) - 1);
  char *number = strtok(copy, " ");
  int i = 0;
  while (number && i < a->rows * a->columns) {
    a->matrix[i / a->columns][i % a->columns] = atof(number);
    i++;
    number = strtok(NULL, " ");
  }
  free(copy);
}

START_TEST(create_mtrx) {
  matrix_t result;
  int err = s21_create_matrix(-1, 4, &result);
  ck_assert_int_eq(err, ERROR);

  err = s21_create_matrix(3, 3, NULL);
  ck_assert_int_eq(err, ERROR);

  err = s21_create_matrix(1, 1, &result);
  ck_assert_int_eq(err, OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);
  ck_assert_double_eq(result.matrix[0][0], 0);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(remove_mtrx) {
  matrix_t mtrx;
  s21_create_matrix(2, 2, &mtrx);
  s21_remove_matrix(&mtrx);
  ck_assert_ptr_null(mtrx.matrix);
  ck_assert_int_eq(mtrx.rows, 0);
  ck_assert_int_eq(mtrx.columns, 0);
}
END_TEST

START_TEST(eq_mtrx) {
  matrix_t a, b;
  s21_create_matrix(3, 2, &a);
  s21_create_matrix(2, 2, &b);
  ck_assert_int_eq(FALSE, s21_eq_matrix(&a, &b));
  s21_remove_matrix(&b);
  s21_create_matrix(3, 2, &b);
  s21_fill_matrix(&a,
                  " 1 2.0 "
                  " 2 2.0 "
                  "-2 0.0 ");
  s21_fill_matrix(&b,
                  " 1 2.0 "
                  " 2 2.0 "
                  "-2 0.0 ");
  ck_assert_int_eq(TRUE, s21_eq_matrix(&b, &a));
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
}
END_TEST

START_TEST(not_eq_matr) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[1][1] = 1;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FALSE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_mtrx) {
  matrix_t a, b, result, zero, c;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &b);
  s21_create_matrix(3, 3, &zero);
  s21_create_matrix(1, 1, &c);
  s21_fill_matrix(&a,
                  "13 4  2 "
                  "-4 2 -1 "
                  " 0 1  1 ");
  s21_fill_matrix(&b,
                  "-13 -4 -2 "
                  " 4  -2  1 "
                  " 0  -1 -1 ");

  int err = s21_sum_matrix(&a, &b, &result);
  ck_assert_int_eq(err, OK);
  ck_assert_int_eq(s21_eq_matrix(&zero, &result), TRUE);
  err = s21_sum_matrix(&a, &c, &result);
  ck_assert_int_eq(err, CALCULATION_ERROR);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&c);
  s21_remove_matrix(&result);
  s21_remove_matrix(&zero);
  err = s21_sum_matrix(&a, &b, &result);
  ck_assert_int_eq(err, ERROR);
}

START_TEST(sub_mtrx) {
  matrix_t a, b, result, zero;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 2, &b);
  s21_create_matrix(3, 3, &zero);
  s21_fill_matrix(&a,
                  "13 4  2 "
                  "-4 2 -1 "
                  " 0 1  1 ");
  s21_sub_matrix(&a, &a, &result);
  ck_assert_int_eq(s21_eq_matrix(&zero, &result), TRUE);
  int err = s21_sub_matrix(&a, &b, &result);
  ck_assert_int_eq(err, CALCULATION_ERROR);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&result);
  s21_remove_matrix(&zero);
  err = s21_sub_matrix(&a, &a, &result);
  ck_assert_int_eq(ERROR, err);
}

START_TEST(mult_num) {
  matrix_t a, result, test;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &test);
  s21_fill_matrix(&a,
                  " 10  20  30 "
                  "-10  25 -30 "
                  " 40 -60  -5 ");
  s21_fill_matrix(&test,
                  " 1    2     3 "
                  "-1  2.5    -3 "
                  " 4   -6  -0.5 ");

  s21_mult_number(&a, 0.1, &result);
  ck_assert_int_eq(s21_eq_matrix(&result, &test), TRUE);
  s21_remove_matrix(&a);
  s21_remove_matrix(&result);
  s21_remove_matrix(&test);
  int err = s21_mult_number(&a, 1, &result);
  ck_assert_int_eq(ERROR, err);
}
END_TEST

START_TEST(mult_mtrx) {
  matrix_t a, b, result, test;
  s21_create_matrix(4, 2, &a);
  s21_create_matrix(2, 2, &b);
  s21_fill_matrix(&a,
                  "14 23 "
                  " 2  5 "
                  "-1 -4 "
                  " 8  3 ");
  s21_fill_matrix(&b,
                  " 3 -3 "
                  " 2  6 ");

  s21_create_matrix(4, 2, &test);
  s21_fill_matrix(&test,
                  " 88  96 "
                  " 16  24 "
                  "-11 -21 "
                  " 30  -6 ");
  s21_mult_matrix(&a, &b, &result);
  ck_assert_int_eq(s21_eq_matrix(&result, &test), TRUE);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&result);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(mult_err_mtrx) {
  matrix_t a, b, result;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(2, 2, &b);
  int err = s21_mult_matrix(&a, &b, &result);
  ck_assert_int_eq(CALCULATION_ERROR, err);
  err = s21_mult_matrix(&a, NULL, &result);
  ck_assert_int_eq(ERROR, err);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
}
END_TEST

START_TEST(mtrx_det) {
  matrix_t a, b;
  double result, test = 1;
  s21_create_matrix(4, 4, &a);
  s21_create_matrix(3, 4, &b);
  int err = s21_determinant(&b, &result);
  s21_remove_matrix(&b);
  ck_assert_int_eq(CALCULATION_ERROR, err);
  s21_determinant(&a, &result);
  ck_assert_int_eq(s21_equal_double(0.0, result), TRUE);
  s21_fill_matrix(&a,
                  "5  7  6  5 "
                  "7 10  8  7 "
                  "6  8 10  9 "
                  "5  7  9 10 ");
  s21_determinant(&a, &result);
  ck_assert_int_eq(s21_equal_double(result, test), TRUE);
  s21_fill_matrix(&a,
                  " 0 4 1 2 "
                  " 4 0 2 0 "
                  " 1 1 2 1 "
                  " 0 2 0 1 ");
  s21_determinant(&a, &result);
  ck_assert_int_eq(s21_equal_double(-4.0, result), TRUE);
  s21_remove_matrix(&a);
  err = s21_determinant(&a, &result);
  ck_assert_int_eq(ERROR, err);
  /*s21_get_random_matrix(1000, 1000, &a);
  s21_determinant(&a, &result);
  s21_remove_matrix(&a);*/
}
END_TEST

START_TEST(mtrx_det_1x1) {
  matrix_t e;
  double result, test = 3.0;
  s21_create_matrix(1, 1, &e);
  s21_fill_matrix(&e, " 3 ");
  s21_determinant(&e, &result);
  ck_assert_int_eq(TRUE, s21_equal_double(result, test));
  s21_remove_matrix(&e);
}
END_TEST

START_TEST(transpose_mtrx) {
  matrix_t a = {0}, ta, test;
  int err = s21_transpose(&a, &ta);
  ck_assert_int_eq(err, ERROR);
  s21_create_matrix(5, 6, &a);
  s21_create_matrix(6, 5, &test);
  s21_fill_matrix(&a,
                  " 1 2 3 4 5 6 "
                  " 2 3 4 5 6 7 "
                  " 3 4 5 6 7 8 "
                  " 4 5 6 7 8 9 "
                  " 5 6 7 8 9 0 ");
  s21_fill_matrix(&test,
                  " 1 2 3 4 5 "
                  " 2 3 4 5 6 "
                  " 3 4 5 6 7 "
                  " 4 5 6 7 8 "
                  " 5 6 7 8 9 "
                  " 6 7 8 9 0 ");
  err = s21_transpose(&a, &ta);
  ck_assert_int_eq(s21_eq_matrix(&test, &ta), TRUE);
  ck_assert_int_eq(err, OK);
  s21_remove_matrix(&a);
  s21_remove_matrix(&ta);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(complem_mtrx) {
  matrix_t a, result, test;
  s21_create_matrix(4, 4, &a);
  s21_create_matrix(4, 4, &test);
  s21_fill_matrix(&a,
                  " 1 -1 -5  4 "
                  "-3  5  2  1 "
                  "10 -11 9 -7 "
                  "-2  0  4  1 ");
  s21_fill_matrix(&test,
                  "  163   45  31 202 "
                  "  268  249 103 124 "
                  "  107   39  44  38 "
                  " -171 -156  81 105 ");
  s21_calc_complements(&a, &result);
  ck_assert_int_eq(s21_eq_matrix(&test, &result), TRUE);
  s21_remove_matrix(&a);
  s21_remove_matrix(&test);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(complem_mtrx_1x1) {
  matrix_t a, result, test;
  s21_create_matrix(1, 1, &a);
  s21_create_matrix(1, 1, &test);
  s21_fill_matrix(&test, " 1.0 ");
  int err = s21_calc_complements(&a, &result);
  ck_assert_int_eq(err, OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &test), TRUE);
  s21_remove_matrix(&a);
  s21_remove_matrix(&result);
  s21_remove_matrix(&test);
}

START_TEST(complem_invalid_input) {
  matrix_t A = {0};
  int err = s21_calc_complements(&A, NULL);
  ck_assert_int_eq(err, ERROR);
}
END_TEST

START_TEST(s21_calc_complements_2) {
  // failure with vector matrix (rows or cols == 1)
  matrix_t A = {};
  matrix_t result = {};
  s21_create_matrix(1, 3, &A);
  ck_assert_int_eq(s21_calc_complements(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_calc_complements_3) {
  // success with task reference values
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 2;
  A.matrix[2][0] = 5, A.matrix[2][1] = 2, A.matrix[2][2] = 1;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  eq_matrix.matrix[0][0] = 0, eq_matrix.matrix[0][1] = 10,
  eq_matrix.matrix[0][2] = -20;
  eq_matrix.matrix[1][0] = 4, eq_matrix.matrix[1][1] = -14,
  eq_matrix.matrix[1][2] = 8;
  eq_matrix.matrix[2][0] = -8, eq_matrix.matrix[2][1] = -2,
  eq_matrix.matrix[2][2] = 4;

  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), TRUE);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(inverse_mtrx) {
  matrix_t a, test, inv, nonsquare;
  s21_create_matrix(2, 1, &nonsquare);
  s21_create_matrix(2, 2, &a);
  s21_create_matrix(2, 2, &test);
  s21_fill_matrix(&a,
                  " 3 4 "
                  " 2 3 ");
  s21_fill_matrix(&test,
                  "  3 -4 "
                  " -2  3 ");

  s21_inverse_matrix(&a, &inv);
  ck_assert_int_eq(TRUE, s21_eq_matrix(&test, &inv));
  s21_remove_matrix(&inv);
  int err = s21_inverse_matrix(&nonsquare, &inv);
  ck_assert_int_eq(err, CALCULATION_ERROR);
  s21_remove_matrix(&a);
  s21_remove_matrix(&test);
  s21_remove_matrix(&inv);
  s21_remove_matrix(&nonsquare);
  err = s21_inverse_matrix(&a, NULL);
  ck_assert_int_eq(err, ERROR);
}
END_TEST

START_TEST(inverse_mtrx_1x1) {
  matrix_t a, test, result;
  s21_create_matrix(1, 1, &a);
  s21_create_matrix(1, 1, &test);
  s21_fill_matrix(&a, " 21 ");
  s21_fill_matrix(&test, " 0.04761904761 ");
  int err = s21_inverse_matrix(&a, &result);
  ck_assert_int_eq(TRUE, s21_eq_matrix(&test, &result));
  ck_assert_int_eq(err, OK);

  s21_remove_matrix(&a);
  s21_remove_matrix(&result);
  s21_remove_matrix(&test);
}

Suite *test_s21_matrix_suite(void) {
  Suite *ret = suite_create("s21_matrix");

  TCase *tc_util = tcase_create("all");

  tcase_add_test(tc_util, s21_calc_complements_2);
  tcase_add_test(tc_util, s21_calc_complements_3);
  tcase_add_test(tc_util, complem_invalid_input);

  tcase_add_test(tc_util, create_mtrx);
  tcase_add_test(tc_util, remove_mtrx);
  tcase_add_test(tc_util, eq_mtrx);
  tcase_add_test(tc_util, not_eq_matr);
  tcase_add_test(tc_util, sum_mtrx);
  tcase_add_test(tc_util, sub_mtrx);
  tcase_add_test(tc_util, mult_num);
  tcase_add_test(tc_util, mult_mtrx);
  tcase_add_test(tc_util, mtrx_det);
  tcase_add_test(tc_util, mtrx_det_1x1);
  tcase_add_test(tc_util, transpose_mtrx);
  tcase_add_test(tc_util, inverse_mtrx);
  tcase_add_test(tc_util, complem_mtrx);
  tcase_add_test(tc_util, complem_mtrx_1x1);
  tcase_add_test(tc_util, mult_err_mtrx);
  tcase_add_test(tc_util, inverse_mtrx_1x1);

  suite_add_tcase(ret, tc_util);
  return ret;
}

int main(void) {
  int failed = 0, total = 0;
  SRunner *sr = srunner_create(test_s21_matrix_suite());
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);

  total = srunner_ntests_run(sr);
  failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  printf("\nTEST RESULT " GREEN "\nTOTAL:\t%d" RED
         "\nFAILED:\t"
         "%d" NOCOLOR "\n",
         total, failed);
  return failed;
}
