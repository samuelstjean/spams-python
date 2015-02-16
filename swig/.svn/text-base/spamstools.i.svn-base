

/* ***************************
   macros to apply typemaps
************************** */


%define APPLY_INPLACE_VECTOR( ctype )
%apply Vector<ctype> *INPLACE_VECTOR {
       inplace_vectors(ctype)
};
%enddef
%define APPLY_ARGOUT_VECTOR( ctype )
%apply Vector<ctype> **ARGOUT_VECTOR {
       argout_vectors(ctype)
};
%enddef


%define APPLY_INPLACE_MATRIX( ctype )
%apply Matrix<ctype> *INPLACE_MATRIX {
       inplace_matrices(ctype)
};
%enddef
%define APPLY_ARGOUT_MATRIX( ctype )
%apply Matrix<ctype> **ARGOUT_MATRIX {
       argout_matrices(ctype)
};
%enddef

%define APPLY_INPLACE_SPMATRIX( ctype )
%apply SpMatrix<ctype> *INPLACE_SPMATRIX {
       inplace_spmatrices(ctype)
};
%enddef
%define APPLY_ARGOUT_SPMATRIX( ctype )
%apply SpMatrix<ctype> **ARGOUT_SPMATRIX {
       argout_spmatrices(ctype)
};
%enddef

%define APPLY_INPLACE_DSPMATRIX( ctype )
%apply AbstractMatrixB<ctype> *INPLACE_DSPMATRIX {
       inplace_dspmatrices(ctype)
};
%enddef
%define APPLY_INPLACE_DATAMATRIX( ctype )
%apply Data<ctype> *INPLACE_DATAMATRIX {
       inplace_datamatrices(ctype)
};
%enddef

%define APPLY_TREE( ctype )
%apply std::vector<StructNodeElem<ctype> *> *TREE {
       inplace_nodes(ctype)
};
%enddef

APPLY_INPLACE_VECTOR(float)
APPLY_INPLACE_VECTOR(double)
APPLY_INPLACE_VECTOR(int)
APPLY_ARGOUT_VECTOR(int)
APPLY_ARGOUT_VECTOR(double)

%apply Matrix<bool> *INPLACE_MATRIX {
       inplace_bool_matrices
};
APPLY_INPLACE_MATRIX(float)
APPLY_INPLACE_MATRIX(double)
APPLY_INPLACE_MATRIX(bool)
APPLY_ARGOUT_MATRIX(double)
APPLY_ARGOUT_MATRIX(float)

%apply SpMatrix<bool> *INPLACE_SPMATRIX {
       inplace_bool_spmatrices
};
APPLY_INPLACE_SPMATRIX(float)
APPLY_INPLACE_SPMATRIX(double)
APPLY_INPLACE_SPMATRIX(bool)
APPLY_ARGOUT_SPMATRIX(bool)
APPLY_ARGOUT_SPMATRIX(double)
APPLY_ARGOUT_SPMATRIX(float)

APPLY_INPLACE_DSPMATRIX(float)
APPLY_INPLACE_DSPMATRIX(double)

APPLY_INPLACE_DATAMATRIX(double)
APPLY_INPLACE_DATAMATRIX(float)

APPLY_TREE(double)
APPLY_TREE(float)

#ifdef SWIGPYTHON
%apply (int* INPLACE_ARRAY1,int DIM1) {
       (int *degr, int n)
};
#endif
%apply (int** ARGOUTVIEW_ARRAY1,int* DIM1) {
       (int **pperm, int *pnb_vars)
};

