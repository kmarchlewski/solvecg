#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

/* Test: Rprintf */
#define PRINT(args)
// #define PRINT(args) args

void swap_ptr(double **ptr1, double **ptr2) {
  double *cache_ptr;
  cache_ptr = *ptr1;
  *ptr1 = *ptr2;
  *ptr2 = cache_ptr;
}

void res(int N, double *A, double *x, double *b, double *r) {
  for (int i = 0; i < N; i++) {
    r[i] = b[i];
    for (int j = 0; j < N; j++)
      r[i] -= A[i + j*N] * x[j];
  }
}

double norm(int N, double *vec) {
  double cache = 0.;
  for (int i = 0; i < N; i++)
    cache += vec[i] * vec[i];
  return sqrt(cache);
}

void mult(int N, double *A, double *x, double *y) {
  for (int i = 0; i < N; i++) {
    y[i] = 0.;
    for (int j = 0; j < N; j++)
      y[i] += A[i + j*N] * x[j];
  }
}

double dprod(int N, double *x, double *y) {
  double cache = 0.;
  for (int i = 0; i < N; i++)
    cache += x[i] * y[i];
  return cache;
}

void subtract(int N, double *vec1, double *vec2, double *vec3) {
  for (int i = 0; i < N; i++)
    vec3[i] = vec1[i] - vec2[i];
}

void precond_i(int N, double *A, double *r, double *p) {
  for (int i = 0; i < N; i++)
    p[i] = r[i];
}

void precond_jacobi(int N, double *A, double *r, double *p) {
  for (int i = 0; i < N; i++)
    p[i] = 1. / A[i + i*N] * r[i];
}

void precond_gs(int N, double *A, double *r, double *p) {
  double cache = 0.;

  p[0] = 1. / A[0] * r[0];

  for (int i = 1; i < N; i++) {
    cache = 0.;
    for (int j = 0; j < i; j++)
      cache += A[i + j*N] * p[j];
    p[i] = 1. / A[i + i*N] * (r[i] - cache);
  }
}

SEXP solve_cg(SEXP mat_A, SEXP vec_x, SEXP vec_b, SEXP str_precond,
              SEXP num_iter_max, SEXP num_eps) {

  SEXP summary, iteration, residuum;
  int N, iter_max;
  double alpha, beta;
  double *A, *x, *b, eps;
  double *r, *p, *p_old, *Ap, *Ap_old, n;
  void (*precond)(int N, double *A, double *r, double *p);
  const char *names[] = {"iteration", "residuum", ""};


  if (!Rf_isMatrix(mat_A)) {
    Rf_error(" 'A' is not a matrix.\n");
    return R_NilValue;
  }

  if (!Rf_isNumeric(mat_A)) {
    Rf_error(" 'A' is not a numeric matrix.\n");
    return R_NilValue;
  }

  if (!(Rf_nrows(mat_A) == Rf_ncols(mat_A))) {
    Rf_error(" 'A' is not a square matrix.\n");
    return R_NilValue;
  }

  N = Rf_nrows(mat_A);

  if (!Rf_isReal(vec_x)) {
    Rf_error(" 'x' is not a numeric vector.\n");
    return R_NilValue;
  }

  if (Rf_length(vec_x) != N) {
    Rf_error(" The length of 'x' and the size of 'A' are not consistent.\n");
    return R_NilValue;
  }

  if (!Rf_isReal(vec_b)) {
    Rf_error(" 'b' is not a numeric vector.\n");
    return R_NilValue;
  }

  if (Rf_length(vec_b) != N) {
    Rf_error(" The length of 'b' and the size of 'A' are not consistent.\n");
    return R_NilValue;
  }

  if (!(Rf_isString(str_precond) && Rf_length(str_precond) == 1)) {
    Rf_error(" 'precond' is not a one string.\n");
    return R_NilValue;
  }

  if (!(strcmp(CHAR(STRING_ELT(str_precond, 0)), "I") == 0 ||
        strcmp(CHAR(STRING_ELT(str_precond, 0)), "Jacobi") == 0 ||
        strcmp(CHAR(STRING_ELT(str_precond, 0)), "GS") == 0)) {
    Rf_error(" Implemented preconditioners are: 'I', 'Jacobi' and 'GS'.\n");
    return R_NilValue;
  }

  if (!(Rf_isNumeric(num_iter_max) && Rf_length(num_iter_max) == 1)) {
    Rf_error(" 'iter_max' should be a number.\n");
    return R_NilValue;
  }

  if (!(Rf_isNumeric(num_eps) && Rf_length(num_eps) == 1)) {
    Rf_error(" 'eps' should be a number.\n");
    return R_NilValue;
  }


  summary = Rf_protect(Rf_mkNamed(VECSXP, names));
  iteration = Rf_protect(Rf_ScalarInteger(0));
  residuum = Rf_protect(Rf_ScalarReal(R_NaN));

  SET_VECTOR_ELT(summary, 0, iteration);
  SET_VECTOR_ELT(summary, 1, residuum);

  if (strcmp(CHAR(STRING_ELT(str_precond, 0)), "Jacobi") == 0) {
    precond = precond_jacobi;
  } else if (strcmp(CHAR(STRING_ELT(str_precond, 0)), "GS") == 0) {
    precond = precond_gs;
  } else {
    precond = precond_i;
  }

  r = (double *) Calloc(N, double);
  p = (double *) Calloc(N, double);
  p_old = (double *) Calloc(N, double);
  Ap = (double *) Calloc(N, double);
  Ap_old = (double *) Calloc(N, double);

  A = REAL(mat_A);
  x = REAL(vec_x);
  b = REAL(vec_b);
  eps = REAL(num_eps)[0];
  iter_max = Rf_asInteger(num_iter_max);

  res(N, A, x, b, r);
  n = norm(N, r);
  PRINT(Rprintf("CG solver: iter = %d, norm(res) = %6.3le\n", 0, n));

  if (n < eps) {
    Free(r);
    Free(p);
    Free(p_old);
    Free(Ap);
    Free(Ap_old);
    REAL(residuum)[0] = n;
    Rf_unprotect(3);
    return summary;
  }

  precond(N, A, r, p);
  mult(N, A, p, Ap);
  alpha = dprod(N, p, r) / dprod(N, p, Ap);
  PRINT(Rprintf("CG solver: iter = %d, alpha = %6.3lf\n", 1, alpha));

  for (int i = 0; i < N; i++) {
    x[i] += alpha * p[i];
  }

  res(N, A, x, b, r);
  n = norm(N, r);
  PRINT(Rprintf("CG solver: iter = %d, norm(res) = %6.3le\n", 1, n));

  if (n < eps) {
    Free(r);
    Free(p);
    Free(p_old);
    Free(Ap);
    Free(Ap_old);
    REAL(residuum)[0] = n;
    INTEGER(iteration)[0] = 1;
    Rf_unprotect(3);
    return summary;
  }

  swap_ptr(&p, &p_old);
  swap_ptr(&Ap, &Ap_old);

  for (int k = 2; k < iter_max + 1; k++) {

    precond(N, A, r, p);
    mult(N, A, p, Ap);
    beta = dprod(N, p_old, Ap) / dprod(N, p_old, Ap_old);

    for (int i = 0; i < N; i++)
      p[i] -= beta * p_old[i];

    alpha = dprod(N, p, r) / dprod(N, p, Ap);
    PRINT(Rprintf("CG solver: iter = %d, alpha = %6.3lf, beta = %6.3lf\n", k, alpha, beta));

    for (int i = 0; i < N; i++) {
      x[i] += alpha * p[i];
    }

    res(N, A, x, b, r);
    n = norm(N, r);
    PRINT(Rprintf("CG solver: iter = %d, norm(res) = %6.3le\n", k, n));

    if (n < eps || k == iter_max) {
      REAL(residuum)[0] = n;
      INTEGER(iteration)[0] = k;
      break;
    }

    swap_ptr(&p, &p_old);
    swap_ptr(&Ap, &Ap_old);

    if (k % 16)
      R_CheckUserInterrupt();
  }

  Free(r);
  Free(p);
  Free(p_old);
  Free(Ap);
  Free(Ap_old);

  Rf_unprotect(3);

  return summary;
}

SEXP invert_cg(SEXP mat_A, SEXP mat_A_inv, SEXP str_precond,
               SEXP num_iter_max, SEXP num_eps) {

  SEXP summary, iteration, residuum;
  int N, iter_max;
  double alpha, beta;
  double *A, *A_inv, eps;
  double *b, *x, *r, *p, *p_old, *Ap, *Ap_old, n;
  void (*precond)(int N, double *A, double *r, double *p);
  const char *names[] = {"iteration", "residuum", ""};


  if (!Rf_isMatrix(mat_A)) {
    Rf_error(" 'A' is not a matrix.\n");
    return R_NilValue;
  }

  if (!Rf_isNumeric(mat_A)) {
    Rf_error(" 'A' is not a numeric matrix.\n");
    return R_NilValue;
  }

  if (Rf_nrows(mat_A) != Rf_ncols(mat_A)) {
    Rf_error(" 'A' is not a square matrix.\n");
    return R_NilValue;
  }

  N = Rf_nrows(mat_A);

  if (!Rf_isMatrix(mat_A_inv)) {
    Rf_error(" 'A_inv' is not a matrix.\n");
    return R_NilValue;
  }

  if (!Rf_isNumeric(mat_A_inv)) {
    Rf_error(" 'A_inv' is not a numeric matrix.\n");
    return R_NilValue;
  }

  if (Rf_nrows(mat_A_inv) != Rf_ncols(mat_A_inv)) {
    Rf_error(" 'A_inv' is not a square matrix.\n");
    return R_NilValue;
  }

  if (Rf_nrows(mat_A_inv) != N) {
    Rf_error(" The sizes of 'A' and A_inv' are not consistent.\n");
    return R_NilValue;
  }

  if (!(Rf_isString(str_precond) && Rf_length(str_precond) == 1)) {
    Rf_error(" 'precond' is not a one string.\n");
    return R_NilValue;
  }

  if (!(strcmp(CHAR(STRING_ELT(str_precond, 0)), "I") == 0 ||
      strcmp(CHAR(STRING_ELT(str_precond, 0)), "Jacobi") == 0 ||
      strcmp(CHAR(STRING_ELT(str_precond, 0)), "GS") == 0)) {
    Rf_error(" Implemented preconditioners are: 'I', 'Jacobi' and 'GS'.\n");
    return R_NilValue;
  }

  if (!(Rf_isNumeric(num_iter_max) && Rf_length(num_iter_max) == 1)) {
    Rf_error(" 'iter_max' should be a number.\n");
    return R_NilValue;
  }

  if (!(Rf_isNumeric(num_eps) && Rf_length(num_eps) == 1)) {
    Rf_error(" 'eps' should be a number.\n");
    return R_NilValue;
  }


  summary = Rf_protect(Rf_mkNamed(VECSXP, names));
  iteration = Rf_protect(Rf_allocVector(INTSXP, N));
  residuum = Rf_protect(Rf_allocVector(REALSXP, N));

  SET_VECTOR_ELT(summary, 0, iteration);
  SET_VECTOR_ELT(summary, 1, residuum);

  if (strcmp(CHAR(STRING_ELT(str_precond, 0)), "Jacobi") == 0) {
    precond = precond_jacobi;
  } else if (strcmp(CHAR(STRING_ELT(str_precond, 0)), "GS") == 0) {
    precond = precond_gs;
  } else {
    precond = precond_i;
  }

  b = (double *) Calloc(N, double);
  x = (double *) Calloc(N, double);
  r = (double *) Calloc(N, double);
  p = (double *) Calloc(N, double);
  p_old = (double *) Calloc(N, double);
  Ap = (double *) Calloc(N, double);
  Ap_old = (double *) Calloc(N, double);

  A = REAL(mat_A);
  A_inv = REAL(mat_A_inv);
  eps = REAL(num_eps)[0];
  iter_max = Rf_asInteger(num_iter_max);


  for (int i = 0; i < N; i++)
    b[i] = 0.;

  for (int col = 0; col < N; col++) {

    b[col-1 < 0 ? 0 : col-1] = 0.;
    b[col] = 1.;

    for (int i = 0; i < N; i++)
      x[i] = A_inv[i + col*N];

    res(N, A, x, b, r);
    n = norm(N, r);
    PRINT(Rprintf("CG invert [%d]: iter = %d, norm(res) = %6.3le\n", col+1, 0, n));

    if (n < eps) {
      REAL(residuum)[col] = n;
      INTEGER(iteration)[col] = 0;
      continue;
    }

    precond(N, A, r, p);
    mult(N, A, p, Ap);
    alpha = dprod(N, p, r) / dprod(N, p, Ap);
    PRINT(Rprintf("CG invert [%d]: iter = %d, alpha = %6.3lf\n", col+1, 1, alpha));

    for (int i = 0; i < N; i++) {
      x[i] += alpha * p[i];
    }

    res(N, A, x, b, r);
    n = norm(N, r);
    PRINT(Rprintf("CG invert [%d]: iter = %d, norm(res) = %6.3le\n", col+1, 1, n));

    if (n < eps) {
      for (int i = 0; i < N; i++)
        A_inv[i + col*N] = x[i];

      REAL(residuum)[col] = n;
      INTEGER(iteration)[col] = 1;
      continue;
    }

    swap_ptr(&p, &p_old);
    swap_ptr(&Ap, &Ap_old);

    for (int k = 2; k < iter_max + 1; k++) {

      precond(N, A, r, p);
      mult(N, A, p, Ap);
      beta = dprod(N, p_old, Ap) / dprod(N, p_old, Ap_old);

      for (int i = 0; i < N; i++)
        p[i] -= beta * p_old[i];

      alpha = dprod(N, p, r) / dprod(N, p, Ap);
      PRINT(Rprintf("CG invert [%d]: iter = %d, alpha = %6.3lf, beta = %6.3lf\n", col+1, k, alpha, beta));

      for (int i = 0; i < N; i++) {
        x[i] += alpha * p[i];
      }

      res(N, A, x, b, r);
      n = norm(N, r);
      PRINT(Rprintf("CG invert [%d]: iter = %d, norm(res) = %6.3le\n", col+1, k, n));

      if (n < eps || k == iter_max) {
        for (int i = 0; i < N; i++)
          A_inv[i + col*N] = x[i];

        REAL(residuum)[col] = n;
        INTEGER(iteration)[col] = k;
        break;
      }

      swap_ptr(&p, &p_old);
      swap_ptr(&Ap, &Ap_old);

      if (k % 16)
        R_CheckUserInterrupt();
    }
  }

  Free(b);
  Free(x);
  Free(r);
  Free(p);
  Free(p_old);
  Free(Ap);
  Free(Ap_old);

  Rf_unprotect(3);

  return summary;
}

SEXP lehmer(SEXP size) {

  SEXP result, matrix, inverse;
  int N;
  double *ptr_matrix, *ptr_inverse;
  const char *names[] = {"matrix", "inverse", ""};

  if (!(Rf_isNumeric(size) && Rf_length(size) == 1)) {
    Rf_error(" 'size' should be a number.\n");
    return R_NilValue;
  }

  N = Rf_asInteger(size);

  result = Rf_protect(Rf_mkNamed(VECSXP, names));
  matrix = Rf_protect(Rf_allocMatrix(REALSXP, N, N));
  inverse = Rf_protect(Rf_allocMatrix(REALSXP, N, N));

  SET_VECTOR_ELT(result, 0, matrix);
  SET_VECTOR_ELT(result, 1, inverse);

  ptr_matrix = REAL(matrix);
  ptr_inverse = REAL(inverse);

  for (int i = 1; i <= N; i++)
    for (int j = 1; j <= N; j++) {
      if (j >= i) {
        ptr_matrix[i-1 + (j-1)*N] = (double)i / j;
      } else {
        ptr_matrix[i-1 + (j-1)*N] = (double)j / i;
      }
    }

  for (int i = 1; i <= N; i++)
    for (int j = 1; j <= N; j++) {

      if (i == j && i < N) {
        ptr_inverse[i-1 + (j-1)*N] = (4.*i*i*i) / (4.*i*i - 1.);
      } else if (i == j && i == N) {
        ptr_inverse[i-1 + (j-1)*N] = N*N / (2.*N - 1.);
      } else if (abs(i - j) == 1 && i < j) {
        ptr_inverse[i-1 + (j-1)*N] = -(i*(i + 1.)) / (2.*i + 1.);
      } else if (abs(i - j) == 1 && i > j) {
        ptr_inverse[i-1 + (j-1)*N] = ptr_inverse[j-1 + (i-1)*N];
      } else {
        ptr_inverse[i-1 + (j-1)*N] = 0.;
      }

    }

  Rf_unprotect(3);

  return result;
}

