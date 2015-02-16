#ifndef SPAMS_TST_H
#define SPAMS_TST_H
#include <cstdlib>

/* TESTS */
double v[6] = {1,3,5,2,4,8};
int indptr[4] = {0,1,3,5};
int r[5] = {0,1,2,0,1};


Matrix<double> *tst(Matrix<double> **path,bool return_path,AbstractMatrixB<double> *D) throw(const char *)
{
  *path = NULL;
  double *data = &v[0];
  int *pB = &indptr[0];
  int *pr = &r[0];
  Matrix<double> * A = new Matrix<double>(data,2,2);
  throw "ERROR";
  if (return_path) {
    Matrix<double> *p = new Matrix<double>(&v[1],1,2);
    *path = p;
  }
  //  sleep(5);
  return A;
}
SpMatrix<double> *xtst(Matrix<double> **path,bool return_path)  throw(const char *)
{
  *path = NULL;

  double *data = &v[0];
  int *pB = &indptr[0];
  int *pr = &r[0];
  SpMatrix<double> * A = new SpMatrix<double>(data,pr,pB,&indptr[1],3,3,5);
  //throw("ERROR\n");
  if (return_path) {
    Matrix<double> *p = new Matrix<double>(&v[1],2,2);
    *path = p;
  }
  cerr << "XTST\n";
  return A;
}
Matrix<double> *ztst(Data<double> *X,Matrix<double> **omA,Matrix<double> **omB,Vector<int> **omiter,bool ret_mat)
{
  double *data = &v[0];
  Matrix<double> * A = new Matrix<double>(data,2,2);
  // throw("ERROR\n");
  //  double *pr = X->rawX();
  //  printf("D= %f\n",*pr);
  if (ret_mat) {
    *omA = new Matrix<double>(&v[1],2,2);
    *omB = new Matrix<double>(&v[2],2,2);
    *omiter = new Vector<int>(1);
    int *p = (*omiter)->rawX();
    *p = 123;
  } else {
    *omA = NULL;
    *omB = NULL;
    *omiter = NULL;
  }
  return A;
}

#endif /* SPAMS_TST_H */
