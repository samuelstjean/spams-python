#include "linalg.h"

double data[] = {1.,2.,3.,4.};
int main(int argc, char**argv) {
  Matrix<double> A(data,2,2);
  A.invSym();
}
