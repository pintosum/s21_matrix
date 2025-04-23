#include "s21_matrix.h"

#include <math.h>

int s21_is_valid_matrix_t(matrix_t *a) {
  return a && a->matrix && a->rows > 0 && a->columns > 0;
}

int s21_is_square_matrix(matrix_t *a) { return a->rows == a->columns; }

int s21_equal_double(double a, double b) {
  return a - b < 1e-7 && a - b > -1e-7;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int ret = OK;
  if (!result || rows <= 0 || columns <= 0) {
    ret = ERROR;
  } else {
    double **matrix = calloc(rows, sizeof(double *));
    if (matrix) {
      for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(columns, sizeof(double));
      }
      result->matrix = matrix;
      result->rows = rows;
      result->columns = columns;
    } else
      ret = ERROR;
  }
  return ret;
}

void s21_remove_matrix(matrix_t *a) {
  if (a && s21_is_valid_matrix_t(a)) {
    for (int i = 0; i < a->rows; i++) {
      free(a->matrix[i]);
    }
    free(a->matrix);
    a->matrix = NULL;
    a->rows = 0;
    a->columns = 0;
  }
}

int s21_equal_dims(matrix_t *a, matrix_t *b) {
  return s21_is_valid_matrix_t(a) && s21_is_valid_matrix_t(b) &&
         a->rows == b->rows && a->columns == b->columns;
}

int s21_eq_matrix(matrix_t *a, matrix_t *b) {
  int ret = FALSE;
  if (s21_equal_dims(a, b)) {
    ret = TRUE;
    for (int i = 0; i < a->rows && ret == TRUE; i++) {
      for (int j = 0; j < a->columns && ret == TRUE; j++) {
        if (!s21_equal_double(a->matrix[i][j], b->matrix[i][j])) {
          ret = FALSE;
        }
      }
    }
  }
  return ret;
}

int s21_sum_matrix(matrix_t *a, matrix_t *b, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a) || !s21_is_valid_matrix_t(b)) {
    ret = ERROR;
  } else if (a->rows == b->rows && a->columns == b->columns) {
    s21_create_matrix(a->rows, a->columns, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < a->columns; j++) {
        result->matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];
      }
    }
  } else {
    ret = CALCULATION_ERROR;
  }
  return ret;
}

int s21_sub_matrix(matrix_t *a, matrix_t *b, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a) || !s21_is_valid_matrix_t(b)) {
    ret = ERROR;
  } else if (a->rows == b->rows && a->columns == b->columns) {
    s21_create_matrix(a->rows, a->columns, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < a->columns; j++) {
        result->matrix[i][j] = a->matrix[i][j] - b->matrix[i][j];
      }
    }
  } else {
    ret = CALCULATION_ERROR;
  }
  return ret;
}

int s21_mult_number(matrix_t *a, double number, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a)) {
    ret = ERROR;
  } else {
    if (a != result) s21_create_matrix(a->rows, a->columns, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < a->columns; j++) {
        result->matrix[i][j] = a->matrix[i][j] * number;
      }
    }
  }
  return ret;
}

int s21_mult_matrix(matrix_t *a, matrix_t *b, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a) || !s21_is_valid_matrix_t(b)) {
    ret = ERROR;
  } else if (a->columns == b->rows) {
    s21_create_matrix(a->rows, b->columns, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < b->columns; j++) {
        for (int r = 0; r < b->rows; r++) {
          result->matrix[i][j] += a->matrix[i][r] * b->matrix[r][j];
        }
      }
    }
  } else {
    ret = CALCULATION_ERROR;
  }
  return ret;
}

int s21_transpose(matrix_t *a, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a))
    ret = ERROR;
  else {
    s21_create_matrix(a->columns, a->rows, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < a->columns; j++) {
        result->matrix[j][i] = a->matrix[i][j];
      }
    }
  }
  return ret;
}

void swap_ptr(void **p, void **r) {
  if (p && r) {
    void *temp = *p;
    *p = *r;
    *r = temp;
  }
}

int s21_new_triangle_matrix(matrix_t *a) {
  int ret = 0;
  int h = 0, k = 0;
  while (h < a->rows && k < a->columns) {
    double max = fabs(a->matrix[h][k]);
    int i_max = h;
    for (int i = i_max + 1; i < a->rows; i++) {
      double t = fabs(a->matrix[i][k]);
      if (t > max) {
        max = t;
        i_max = i;
      }
    }
    if (s21_equal_double(a->matrix[i_max][k], 0.0)) {
      k++;
    } else {
      if (h != i_max) {
        swap_ptr((void **)&a->matrix[h], (void **)&a->matrix[i_max]);
        ret++;
      }
      for (int i = h + 1; i < a->rows; i++) {
        double f = a->matrix[i][k] / a->matrix[h][k];
        a->matrix[i][k] = 0;
        for (int j = k + 1; j < a->columns; j++) {
          a->matrix[i][j] -= a->matrix[h][j] * f;
        }
      }
      h++;
      k++;
    }
  }
  return ret;
}

double s21_minor_matrix_det(matrix_t *a, int i, int j) {
  double ret = 0.0;
  matrix_t minor = {0};
  s21_create_matrix(a->rows - 1, a->columns - 1, &minor);
  int rows = 0;
  int columns = 0;
  for (int k = 0; k < a->rows; k++) {
    if (k == i) continue;
    for (int m = 0; m < a->columns; m++) {
      if (m == j) continue;
      minor.matrix[rows][columns] = a->matrix[k][m];
      columns++;
    }
    columns = 0;
    rows++;
  }
  if (rows == 0) {
    ret = 1.0;
  } else {
    s21_determinant(&minor, &ret);
    ret = (i + j) % 2 == 1 ? -ret : ret;
  }
  s21_remove_matrix(&minor);
  if (s21_equal_double(0.0, ret)) ret *= ret;
  return ret;
}

int s21_calc_complements(matrix_t *a, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a)) {
    ret = ERROR;
  } else if (s21_is_square_matrix(a)) {
    s21_create_matrix(a->rows, a->columns, result);
    for (int i = 0; i < a->rows; i++) {
      for (int j = 0; j < a->columns; j++) {
        result->matrix[i][j] = s21_minor_matrix_det(a, i, j);
      }
    }
  } else {
    ret = CALCULATION_ERROR;
  }
  return ret;
}

matrix_t s21_copy_matrix(matrix_t *a) {
  matrix_t ret = {0};
  s21_create_matrix(a->rows, a->columns, &ret);
  for (int i = 0; i < a->rows; i++) {
    for (int j = 0; j < a->columns; j++) {
      ret.matrix[i][j] = a->matrix[i][j];
    }
  }
  return ret;
}

int s21_determinant(matrix_t *a, double *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a))
    ret = ERROR;
  else if (s21_is_square_matrix(a)) {
    double det = 1.0;
    matrix_t copy = s21_copy_matrix(a);
    int swap = s21_new_triangle_matrix(&copy);
    for (int i = 0; i < copy.rows; i++) {
      det *= copy.matrix[i][i];
    }

    s21_remove_matrix(&copy);
    if (swap % 2 == 1) {
      det = -det;
    }
    if (s21_equal_double(0.0, det)) det *= det;
    *result = det;
  } else
    ret = CALCULATION_ERROR;
  return ret;
}

int s21_inverse_matrix(matrix_t *a, matrix_t *result) {
  int ret = OK;
  if (!result || !s21_is_valid_matrix_t(a)) {
    ret = ERROR;
  } else {
    double det = 0;
    ret = s21_determinant(a, &det);
    int det_is_zero = s21_equal_double(det, 0.0);
    if (!det_is_zero && ret == OK && a->rows == 1) {
      s21_create_matrix(1, 1, result);
      result->matrix[0][0] = 1. / det;
    } else if (!det_is_zero && ret == OK) {
      matrix_t temp = {0};
      s21_calc_complements(a, result);
      s21_transpose(result, &temp);
      s21_remove_matrix(result);
      s21_mult_number(&temp, 1. / det, result);
      s21_remove_matrix(&temp);
    } else {
      ret = CALCULATION_ERROR;
    }
  }
  return ret;
}
