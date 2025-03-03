# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True
# cython: language_level = 3

from scipy.linalg.cython_blas cimport *
from scipy.linalg.cython_lapack cimport *


cdef public double dnrm2_(int *n,double *x,int *incX) noexcept nogil:
    return dnrm2(n, x, incX)

cdef public float snrm2_(int *n,float *x,int *incX) noexcept nogil:
    return snrm2(n, x, incX)

cdef public void dcopy_(int *n,double *x,int *incX, double *y,int *incY) noexcept nogil:
    dcopy(n,x, incX, y, incY)

cdef public void scopy_(int *n,float *x,int *incX, float *y,int *incY) noexcept nogil:
    scopy(n, x, incX, y, incY)

cdef public void daxpy_(int *n,double* a, double *x,int *incX, double *y,int *incY) noexcept nogil:
    daxpy(n, a, x, incX, y, incY)

cdef public void saxpy_(int *n,float* a, float *x,int *incX, float *y,int *incY) noexcept nogil:
    saxpy(n, a,  x, incX, y, incY)

# cdef public void daxpby_(int *n,double* a, double *x,int *incX,double* b, double *y,int *incY) noexcept nogil:
#     daxpby(n, a, x, incX, b, y, incY)

# cdef public void saxpby_(int *n,float* a, float *x,int *incX, float* b, float *y,int *incY) noexcept nogil:
#     saxpby(n, a,  x, incX, b, y, incY)

cdef public void dscal_(int *n,double* a, double *x,int *incX) noexcept nogil:
    dscal(n, a, x, incX)

cdef public void sscal_(int *n,float* a, float *x,int *incX) noexcept nogil:
    sscal(n, a,  x, incX)

cdef public double dasum_(int *n,double *x,int *incX) noexcept nogil:
    return dasum(n,x, incX)

cdef public float sasum_(int *n,float *x,int *incX) noexcept nogil:
    return sasum(n, x, incX)

cdef public double ddot_(int *n,double *x,int *incX, double *y,int *incY) noexcept nogil:
    return ddot(n,x, incX, y, incY)

cdef public float sdot_(int *n,float *x,int *incX, float *y,int *incY) noexcept nogil:
    return sdot(n, x, incX, y, incY)

cdef public void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y,int *incy) noexcept nogil:
    dgemv(trans, m, n,  alpha,  a, lda, x, incx,  beta, y, incy)

cdef public void sgemv_(char *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y,int *incy) noexcept nogil:
    sgemv(trans, m, n,  alpha,  a, lda,  x, incx, beta, y, incy)

cdef public void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda) noexcept nogil:
    dger(m, n,  alpha, x, incx, y,  incy,  a, lda)

cdef public void sger_(int *m, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *a, int *lda) noexcept nogil:
    sger(m, n,  alpha,  x, incx, y,  incy,  a, lda)

cdef public void dtrmv_(char *uplo, char *trans, char *diag, int *n, double *a, int *lda, double *x, int *incx) noexcept nogil:
    dtrmv(uplo, trans, diag, n,  a, lda, x, incx)

cdef public void strmv_(char *uplo, char *trans, char *diag, int *n, float *a, int *lda, float *x, int *incx) noexcept nogil:
    strmv(uplo, trans, diag, n,  a, lda,  x, incx)

cdef public void dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda) noexcept nogil:
    dsyr(uplo, n,  alpha, x, incx,  a, lda)

cdef public void ssyr_(char *uplo, int *n, float *alpha, float *x, int *incx, float *a, int *lda) noexcept nogil:
    ssyr(uplo, n,  alpha,  x, incx,  a, lda)

cdef public void dsymv_(char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy) noexcept nogil:
    dsymv(uplo, n,  alpha,  a, lda, x, incx,  beta, y,  incy)

cdef public void ssymv_(char *uplo, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy) noexcept nogil:
    ssymv(uplo, n,  alpha,  a, lda,  x, incx, beta, y,  incy)

cdef public void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc) noexcept nogil:
    dgemm( transa, transb, m, n, k,  alpha,  a, lda,  b, ldb,  beta, c, ldc)

cdef public void sgemm_(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc) noexcept nogil:
    sgemm( transa, transb, m, n, k,  alpha,  a, lda, b, ldb, beta, c, ldc)

cdef public void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc) noexcept nogil:
    dsyrk(uplo, trans, n, k,  alpha,  a, lda,  beta, c, ldc)

cdef public void ssyrk_(char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *beta, float *c, int *ldc) noexcept nogil:
    ssyrk(uplo, trans, n, k,  alpha,  a, lda, beta, c, ldc)

cdef public void dtrmm_(char *side,char *uplo,char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb) noexcept nogil:
    dtrmm(side,uplo, transa, diag, m, n,  alpha,  a, lda,  b, ldb)

cdef public void strmm_(char *side,char *uplo,char *transa, char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb) noexcept nogil:
    strmm(side, uplo,  transa, diag, m, n,  alpha,  a, lda, b, ldb)

cdef public int idamax_(int *n, double *dx, int *incx) noexcept nogil:
    return idamax(n, dx, incx)

cdef public int isamax_(int *n, float *dx, int *incx) noexcept nogil:
    return isamax(n, dx, incx)

cdef public void dtrtri_(char* uplo, char* diag, int* n, double * a, int* lda, int* info) noexcept nogil:
    dtrtri(uplo, diag, n,  a, lda, info)

cdef public void strtri_(char* uplo, char* diag, int* n, float * a, int* lda, int* info) noexcept nogil:
    strtri(uplo, diag, n,  a, lda, info)

cdef public void dsytrf_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info) noexcept nogil:
    dsytrf(uplo, n,  a, lda, ipiv, work, lwork, info)

cdef public void ssytrf_(char* uplo, int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info) noexcept nogil:
    ssytrf(uplo, n,  a, lda, ipiv, work, lwork, info)

cdef public void dsytri_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work, int* info) noexcept nogil:
    dsytri(uplo, n,  a, lda, ipiv, work, info)

cdef public void ssytri_(char* uplo, int* n, float* a, int* lda, int* ipiv, float* work, int* info) noexcept nogil:
    ssytri(uplo, n,  a, lda, ipiv, work, info)

cdef public void dlasrt_(char* id, int* n, double *d, int* info) noexcept nogil:
    dlasrt(id, n, d, info)

cdef public void slasrt_(char* id, int* n, float*d, int* info) noexcept nogil:
    slasrt(id, n, d, info)

cdef public void dgesvd_(char*jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info) noexcept nogil:
    dgesvd(jobu, jobvt, m, n,  a, lda, s, u, ldu, vt, ldvt, work, lwork, info)

cdef public void sgesvd_(char*jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info) noexcept nogil:
    sgesvd(jobu, jobvt, m, n,  a, lda, s, u, ldu, vt, ldvt, work, lwork, info)

cdef public void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info) noexcept nogil:
    dsyev(jobz, uplo, n,  a, lda, w, work, lwork, info)

cdef public void ssyev_(char *jobz, char *uplo, int *n, float *a, int *lda, float *w, float *work, int *lwork, int *info) noexcept nogil:
    ssyev(jobz, uplo, n,  a, lda, w, work, lwork, info)
