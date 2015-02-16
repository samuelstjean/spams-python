%define DOCSTRING
"This module gives access to some functions of the spams C++ library.
The functions defined here should not be called directly.
Use of spams functions should only be done through module spams."
%enddef

#ifdef SWIGPYTHON
%module(docstring=DOCSTRING) spams_wrap
#endif
#ifdef SWIGR
%module spams
#endif

%{
#define SWIG_FILE_WITH_INIT

#include "spams.h"
#include "spams/prox/groups-graph.h"

#ifdef DEBUG
#include "spams-tst.h"
#endif
%}

%define argout_vectors(ctype)
	Vector<ctype> **omiter,
	Vector<ctype> **peta_g,
	Vector<ctype> **pown_variables,
	Vector<ctype> **pN_own_variables
%enddef

/*
   Following macros define which typemap must be applied to arguments of C++ functions are converted.
   args of type inplace_* may containt input data or are ready (allocated)
    to receive output data. Typemaps of type 'in' are applied to them
    before function call. 
   args of type argout_* are used when multiple return values are needed.
    they are of type **; their storage is allocated on return of the C++ function.
     Typemaps of type 'in' are applied to them  before function call,
     and typemaps of type 'argout' are applied to them afeter the call.
  Note : typemaps of type ('out') are applied to return values of C++ functions.
*****

*/
// list of arguments of type INPLACE_MATRIX
%define inplace_bool_matrices
	Matrix<bool> *B
%enddef
%define inplace_bool_spmatrices
    SpMatrix<bool> *groups,
    SpMatrix<bool> *groups_var
%enddef

%define inplace_matrices(ctype)
     Matrix<ctype> *A,
     Matrix<ctype> *B,
     Matrix<ctype> *D,
     Matrix<ctype> *D1,
     Matrix<ctype> *Q,
     Matrix<ctype> *q,
     Matrix<ctype> *U,
     Matrix<ctype> *V,
     Matrix<ctype> *X,
     Matrix<ctype> *m_A,
     Matrix<ctype> *m_B,
     Matrix<ctype> *XAt,
     Matrix<ctype> *Y,
     Matrix<ctype> *XY,
     Matrix<ctype> *AAt,
     Matrix<ctype> *alpha0,
     Matrix<ctype> *alpha,
     Matrix<ctype> *W,
     Matrix<ctype> *Z0,
     Matrix<ctype> *Z
%enddef
%define argout_matrices(ctype)
     Matrix<ctype> **path,
     Matrix<ctype> **omA,
     Matrix<ctype> **omB
%enddef
%define inplace_vectors(ctype)
    Vector<ctype> *v,
    Vector<ctype> *b,
    Vector<ctype> *x,
    Vector<ctype> *xCurr,
    Vector<ctype> *inner_weights,
    Vector<ctype> *L,
    Vector<ctype> *eps,
    Vector<ctype> *Lambda,
    Vector<ctype> *eta_g,
    Vector<ctype> *own_variables,
    Vector<ctype> *N_own_variables,
    Vector<ctype> *groups
%enddef
%define inplace_spmatrices(ctype)
    SpMatrix<ctype> *A,
    SpMatrix<ctype> *alpha,
    SpMatrix<ctype> *groups
%enddef
%define argout_spmatrices(ctype)
    SpMatrix<ctype> **pgroups,
    SpMatrix<ctype> **pgroups_var,
    SpMatrix<ctype> **spA,
    SpMatrix<ctype> **spB
%enddef

%define inplace_dspmatrices(ctype)
    AbstractMatrixB<ctype> *D
%enddef

%define inplace_datamatrices(ctype)
    Data<ctype> *X
%enddef

%define inplace_nodes(ctype)
	std::vector<StructNodeElem<ctype> *> *gstruct
%enddef

#ifdef SWIGPYTHON
%include "py_typemaps.i"
%init %{
    import_array();
%}
#endif

#ifdef SWIGR
%include "R_typemaps.i"
#endif
%include "spamstools.i"
%include "exception.i"


%include <spams.h>
%include <spams/prox/groups-graph.h>

//void im2col_sliding(Matrix<double>  *,Matrix<double>  *,int,int,bool);
//std::vector<NodeElem *> *simpleGroupTree(int *degr, int n)throw(const char *);
//std::vector<NodeElem *> *readGroupStruct(const char* file) throw(const char *);
//std::vector<NodeElem *> *groupStructOfString(const char* data) throw(const char *);
//Vector<double> * graphOfGroupStruct(std::vector<NodeElem *> *tree,SpMatrix<bool> **pgroups,SpMatrix<bool> **pgroups_var) throw(const char *);

//int  treeOfGroupStruct(std::vector<StructNodeElem<double> *> *tree,int **pperm,int *pnb_vars,Vector<double> **peta_g,SpMatrix<bool> **pgroups,Vector<int> **pown_variables,Vector<int> **pN_own_variables) throw(const char *);

#ifdef DEBUG
%include <spams-tst.h>
Matrix<double> *tst(Matrix<double> **,bool,AbstractMatrixB<double> *);
//int xtst(Matrix<double> **,bool );
SpMatrix<double> *xtst(Matrix<double> **,bool );
Matrix<double> *ztst(Data<double> *,Matrix<double> **omA,Matrix<double> **,Vector<int> **,bool) throw(const char *);

#endif

// linalg
INSTANTIATE_DATA(sort)
INSTANTIATE_DATA(mult)
INSTANTIATE_DATA(AAt)
INSTANTIATE_DATA(XAt)
INSTANTIATE_DATA(applyBayerPattern)
INSTANTIATE_DATA(conjugateGradient)
INSTANTIATE_DATA(invSym)
INSTANTIATE_DATA(normalize)

/**** decomp ****/
enum constraint_type { L1COEFFS, L2ERROR, PENALTY, SPARSITY, L2ERROR2, PENALTY2, FISTAMODE};

INSTANTIATE_DATA(sparseProject)
INSTANTIATE_DATA(lassoD)
INSTANTIATE_DATA(lassoQq)
INSTANTIATE_DATA(lassoMask)
INSTANTIATE_DATA(lassoWeighted)
INSTANTIATE_DATA(omp)
INSTANTIATE_DATA(ompMask)
INSTANTIATE_DATA(somp)
INSTANTIATE_DATA(cd)
INSTANTIATE_DATA(l1L2BCD)

/**** dictLearn ****/
enum constraint_type_D { L2,  L1L2, L1L2FL, L1L2MU};
INSTANTIATE_DATA(alltrainDL)
/* from arch */
INSTANTIATE_DATA(archetypalAnalysis)
INSTANTIATE_DATA(archetypalAnalysisInit)
INSTANTIATE_DATA(decompSimplex)


/**** prox ****/
INSTANTIATE_DATA(fistaFlat)
INSTANTIATE_DATA(fistaTree)
INSTANTIATE_DATA(fistaGraph)
INSTANTIATE_DATA(proximalFlat)
INSTANTIATE_DATA(proximalTree)
INSTANTIATE_DATA(proximalGraph)

/* groups-graph */
INSTANTIATE_DATA(simpleGroupTree)
INSTANTIATE_DATA(readGroupStruct)
INSTANTIATE_DATA(groupStructOfString)
INSTANTIATE_DATA(graphOfGroupStruct)
INSTANTIATE_DATA(treeOfGroupStruct)


/* misc */
INSTANTIATE_DATA(im2col_sliding)

