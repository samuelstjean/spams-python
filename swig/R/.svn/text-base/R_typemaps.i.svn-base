%{
extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rembedded.h>
#ifndef WIN32
#include <Rinterface.h>
#endif
#include <R_ext/RS.h>
#include <R_ext/Error.h>
}
#include "spams.h"

/* ********************
* this is a hack to handle multiple values outrput
* return value is a vector initailized by swig in typemap out
* R_result_pos is the index in the vector : it must be set to 0 in typemap(out)
* typemap(argout) cannot be used without typemap(ou) in a function
*/

static int R_result_pos = 0;

SEXP appendOutput(SEXP value,SEXP result) {
     R_result_pos++;
     if(LENGTH(result) > R_result_pos)
     	SET_VECTOR_ELT(result,R_result_pos,value);
     return result;
}
/*    end of hack */
%}

#define myerr(msg,n) {Rf_error(msg,n); return R_NilValue;}
#define MYADD_OUTPUT_ARG(result, value)  r_ans = appendOutput(value, R_OutputValues);

%typemap(throws) const char * %{
	Rf_error("Runtime Error %s",$1); 
	return R_NilValue;
%}


/* One dimensional input arrays */
%define %vector_typemaps(R_TYPE,R_CAST,DATA_TYPE)

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
    $1 = (TYPEOF($input) == R_TYPE && Rf_isVector($input) ) ? 1 : 0;
}

%typemap(in) (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
    SEXP rvec=$input;
    if (TYPEOF(rvec) != R_TYPE || ! Rf_isVector(rvec))
    {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected DATA_TYPE Vector as argument %d",$argnum);
    }

    $1 = new Vector<DATA_TYPE>((DATA_TYPE*) R_CAST(rvec), LENGTH(rvec));
}
%typemap(out) (Vector<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t n = result->n();
    DATA_TYPE *data = result->rawX();
    SEXP rmat;
    PROTECT(rmat = Rf_allocMatrix(REALSXP,1,n));
   if (!rmat) error_return("Cannot alloc R matrix for return value");
   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
   memcpy(rdata,data,n * sizeof(DATA_TYPE));
   delete result;
   $result = rmat;
   UNPROTECT(1);
 
}
%typemap(freearg) (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
	delete arg$argnum;
}
//  ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Vector<DATA_TYPE> **ARGOUT_VECTOR)
{
	$1 = 1;

}
%typemap(in,numinputs=0) (Vector<DATA_TYPE> **ARGOUT_VECTOR)
	     (Vector<DATA_TYPE>  *data_temp)
{
    $1 = &data_temp;
}
%typemap(argout) (Vector<DATA_TYPE> **ARGOUT_VECTOR ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    R_len_t n = data_temp$argnum->n();
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    SEXP rmat;
	    PROTECT(rmat = Rf_allocMatrix(R_TYPE,1,n));
	   if (!rmat) myerr("Cannot alloc R matrix for arg %d",$argnum);
	   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
	   memcpy(rdata,data,n * sizeof(DATA_TYPE));
	   delete data_temp$argnum;
	   MYADD_OUTPUT_ARG($result,rmat);
	   UNPROTECT(1);
        }
}

%enddef /* %vector_typemaps */

%define map_matrix(DATA_TYPE,R_TYPE,R_CAST)
     SEXP rmat=$input;
     SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
    if (TYPEOF(rmat) != R_TYPE || LENGTH(dims) != 2)	
    {	
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected DATA_TYPE dense matrix as argument %d",$argnum);
    }
    $1 = new Matrix<DATA_TYPE> ((DATA_TYPE *)R_CAST(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

%enddef /* map_matrix */

/* full matrix input/output */
%typemap(in) (Matrix<bool> *INPLACE_MATRIX)
{
     SEXP rmat=$input;
     SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
    if (TYPEOF(rmat) != LGLSXP || LENGTH(dims) != 2)	
    {	
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected bool dense matrix as argument %d",$argnum);
    }
      	 $1 = new Matrix<bool> (Rf_nrows(rmat),Rf_ncols(rmat));
	 int *pi = (int *)LOGICAL(rmat);
	 bool *po = $1->rawX();
	 for(int i =0;i < Rf_nrows(rmat) * Rf_ncols(rmat);i++)
	 	*po++ = (bool) *pi++;
         
}
%typemap(freearg)
  (Matrix<bool> *INPLACE_MATRIX)
{
	delete arg$argnum;
}

%define %matrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
           (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	$1 = (TYPEOF($input) == R_TYPE) ? 1 : 0;

}
%typemap(in) (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
     map_matrix(DATA_TYPE,R_TYPE,R_CAST)
}
%typemap(out) (Matrix<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t m = result->m();
    R_len_t n = result->n();
    DATA_TYPE *data = result->rawX();
    SEXP rmat;
    PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
   if (!rmat) error_return("Cannot alloc R matrix for return value");
   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
   memcpy(rdata,data,m * n * sizeof(DATA_TYPE));
   delete result;
   $result = rmat;
   UNPROTECT(1);
 
}

%typemap(freearg)
  (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	delete arg$argnum;
}
//  ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
{
	$1 = 1;

}
%typemap(in,numinputs=0) (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
	     (Matrix<DATA_TYPE>  *data_temp)
{
    $1 = &data_temp;
}
%typemap(argout) (Matrix<DATA_TYPE> **ARGOUT_MATRIX ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    R_len_t m = data_temp$argnum->m();
	    R_len_t n = data_temp$argnum->n();
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    SEXP rmat;
	    PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
	   if (!rmat) myerr("Cannot alloc R matrix for arg %d",$argnum);
	   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
	   memcpy(rdata,data,m * n * sizeof(DATA_TYPE));
	   delete data_temp$argnum;
	   MYADD_OUTPUT_ARG($result,rmat);
	   UNPROTECT(1);
        }
}

%enddef /* %matrix_typemaps */

%fragment("DSp_Check","header")
{
    const int check_sparse(SEXP arg) {
	SEXP ivec = Rf_getAttrib(arg,Rf_install("i"));
	SEXP xvec = Rf_getAttrib(arg,Rf_install("x"));
	SEXP dims = Rf_getAttrib(arg,Rf_install("Dim"));
//	return (Rf_isInteger(ivec) && (TYPEOF(xvec) == R_TYPE) && (LENGTH(dims) == 2) ? 1 : 0);
	return (Rf_isInteger(ivec) && (TYPEOF(xvec) == REALSXP) && (LENGTH(dims) == 2) ? 1 : 0);

    }
    const int check_matrix(SEXP arg) {
    	  SEXP dims = Rf_getAttrib(arg,Rf_install("dim"));
	  return ((TYPEOF(arg) == REALSXP) && (LENGTH(dims) == 2) ? 1 : 0);
    }
}
%define map_sparse(DATA_TYPE,R_TYPE,R_CAST)
     SEXP spmat=$input;
     SEXP ivec = Rf_getAttrib(spmat,Rf_install("i"));
     SEXP pvec = Rf_getAttrib(spmat,Rf_install("p"));
     SEXP xvec = Rf_getAttrib(spmat,Rf_install("x"));
     SEXP dims = Rf_getAttrib(spmat,Rf_install("Dim"));

     if (TYPEOF(xvec) != R_TYPE || LENGTH(dims) != 2)	
    {	
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected DATA_TYPE sparse matrix as argument %d",$argnum);
    }
     int *pB = INTEGER(pvec);
     int *pE = pB + 1;
     int *idims = INTEGER(dims);
     $1 = new SpMatrix<DATA_TYPE> ((DATA_TYPE *)R_CAST(xvec),INTEGER(ivec),pB,pE,idims[0],idims[1],LENGTH(xvec));
	
%enddef /* map_sparse */

/* special case for type bool (32bits in R , 8 bits in C++) */
%typemap(in) (SpMatrix<bool> *INPLACE_SPMATRIX)
{
     SEXP spmat=$input;
     SEXP ivec = Rf_getAttrib(spmat,Rf_install("i"));
     SEXP pvec = Rf_getAttrib(spmat,Rf_install("p"));
     SEXP xvec = Rf_getAttrib(spmat,Rf_install("x"));
     SEXP dims = Rf_getAttrib(spmat,Rf_install("Dim"));

     if (TYPEOF(xvec) != LGLSXP || LENGTH(dims) != 2)	
    {	
        myerr("Expected bool sparse matrix as argument %d",$argnum);
    }
    int nzmax = LENGTH(xvec);
     $1 = new SpMatrix<bool> (idims[0],idims[1],nzmax);
     int *pB = INTEGER(pvec);
     int *pE = pB + 1;
     int *idims = INTEGER(dims);
     int *pi = (int *)LOGICAL(xvec);
     memcpy($1->pB(),INTEGER(ivec),(idims[1] + 1) * sizeof(int));
     memcpy($1->r(),INTEGER(pvec),nzmax * sizeof(int));
     bool *po = $1->v();
     for(int i =0;i < nzmax;i++)
	*po++ = (bool) *pi++;
}
%typemap(freearg)
  (SpMatrix<bool> *INPLACE_SPMATRIX)
{
	delete arg$argnum;
}

%define %spmatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	$1 = check_sparse((SEXP) $input);

}
%typemap(in) (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
     map_sparse(DATA_TYPE,R_TYPE,R_CAST)
}
%typemap(freearg)
  (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	delete arg$argnum;
}

%typemap(out) (SpMatrix<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t m = result->m();
    R_len_t n = result->n();
    R_len_t nzmax = result->nzmax();
    SEXP indptr, indices, vdata, dims, output;
    PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
    PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
    PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
    PROTECT(dims = Rf_allocVector(VECSXP,2));
    SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
    SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

    DATA_TYPE *xdata = result->v();
    memcpy(R_CAST(vdata),xdata,nzmax * sizeof(DATA_TYPE));
    int *pB = result->pB();
    int *r = result->r();
    memcpy(INTEGER(indices),r,nzmax * sizeof(int));
    memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));
    
    PROTECT(output = Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(output,0,indptr);
    SET_VECTOR_ELT(output,1,indices);
    SET_VECTOR_ELT(output,2,vdata);
    SET_VECTOR_ELT(output,3,dims);
    delete result;
    $result = output;
    UNPROTECT(5);
}

%typemap(freearg)
  (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	delete arg$argnum;
}
// ARGOUT
%typemap(in,numinputs=0) (SpMatrix<DATA_TYPE> **ARGOUT_SPMATRIX ) 
(SpMatrix<DATA_TYPE> *data_temp)
{
	$1 = &data_temp;
}
%typemap(argout) (SpMatrix<DATA_TYPE> **ARGOUT_SPMATRIX ) 
{
  if(data_temp$argnum != NULL) {
    R_len_t m = data_temp$argnum->m();
    R_len_t n = data_temp$argnum->n();
    R_len_t nzmax = data_temp$argnum->nzmax();
    SEXP indptr, indices, vdata, dims, output;
    PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
    PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
    PROTECT(vdata = Rf_allocVector(R_TYPE,nzmax));
    PROTECT(dims = Rf_allocVector(VECSXP,2));
    SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
    SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

    DATA_TYPE *xdata = data_temp$argnum->v();
    memcpy(R_CAST(vdata),xdata,nzmax * sizeof(DATA_TYPE));
    int *pB = data_temp$argnum->pB();
    int *r = data_temp$argnum->r();
    memcpy(INTEGER(indices),r,nzmax * sizeof(int));
    memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));
    delete data_temp$argnum;
    PROTECT(output = Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(output,0,indptr);
    SET_VECTOR_ELT(output,1,indices);
    SET_VECTOR_ELT(output,2,vdata);
    SET_VECTOR_ELT(output,3,dims);
    MYADD_OUTPUT_ARG($result,output);
    UNPROTECT(5);
  }
}
%enddef /* %spmatrix_typemaps */

%define %dspmatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
	if(pvec == R_NilValue) {
	    $1 = check_matrix((SEXP) $input);
	 } else {
	    $1 = check_sparse((SEXP) $input);
	 }

}
%typemap(in) (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
    SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
    if(pvec ==  R_NilValue) {
       map_matrix(DATA_TYPE,R_TYPE,R_CAST)
    } else {
       map_sparse(DATA_TYPE,R_TYPE,R_CAST)
    }
}
%typemap(freearg)
  (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	delete arg$argnum;
}

%enddef /* %dspmatrix_typemaps */

%define %datamatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
	if(pvec == R_NilValue) {
	    $1 = check_matrix((SEXP) $input);
	 } else {
	    $1 = check_sparse((SEXP) $input);
	 }

}
%typemap(in) (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
    SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
    if(pvec ==  R_NilValue) {
       map_matrix(DATA_TYPE,R_TYPE,R_CAST)
    } else {
       map_sparse(DATA_TYPE,R_TYPE,R_CAST)
    }
}
%typemap(freearg)
  (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	delete arg$argnum;
}

%enddef /* %datamatrix_typemaps */

%define %node_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typemap(out) (std::vector<StructNodeElem<DATA_TYPE> *> *)
{	      
  int n = result->size();
  SEXP node_list;
  PROTECT(node_list = Rf_allocVector(VECSXP,n));
  int i = 0;
  for(std::vector<StructNodeElem<DATA_TYPE> *>::iterator it = result->begin();it != result->end();it++, i++) {
    SEXP vars, children, elem;
    PROTECT(elem = Rf_allocVector(VECSXP,4));
    StructNodeElem<DATA_TYPE> *node = *it;
    int inode = node->node_num;
    SET_VECTOR_ELT(elem,0,Rf_ScalarInteger(inode));
    SET_VECTOR_ELT(elem,1,Rf_ScalarReal(node->weight));
    int k = node->vars->size();
    PROTECT(vars = Rf_allocVector(INTSXP,k));
    int * rvars = (int *)INTEGER(vars);
    memcpy(rvars,node->vars->data(),k * sizeof(int));
    SET_VECTOR_ELT(elem,2,vars);

    k = node->children->size();
    PROTECT(children = Rf_allocVector(INTSXP,k));
    int * rchildren = (int *)INTEGER(children);
    memcpy(rchildren,node->children->data(),k * sizeof(int));
    SET_VECTOR_ELT(elem,3,children);
    
    SET_VECTOR_ELT(node_list,i,elem);
    UNPROTECT(3);
  }
  del_gstruct(result);
  $result = node_list;
  UNPROTECT(1);
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    (std::vector<StructNodeElem<DATA_TYPE> *> *TREE)
{
    $1 = (TYPEOF($input) == VECSXP && Rf_isVector($input) ) ? 1 : 0;
}
%typemap(in)
  (std::vector<StructNodeElem<DATA_TYPE> *> *TREE) 
{
  SEXP rtree = $input;
  if(TYPEOF(rtree) != VECSXP || ! Rf_isVector(rtree)) {
    myerr("Expected VECSXP Vector as argument %d",$argnum);
  }
  $1 = new std::vector<StructNodeElem<DATA_TYPE> *>;
  int n = Rf_length(rtree);
  for (int i = 0;i < n;i++) {
    SEXP rnode = VECTOR_ELT(rtree,i);
    if(TYPEOF(rnode) != VECSXP || ! Rf_isVector(rnode)) {
      myerr("Bad type for node elem in arg %d",$argnum);
    }
    if(Rf_length(rnode) != 4) {
      myerr("A node must have 4 elements in arg %d",$argnum);
    }
    int inode = *((int *) INTEGER(VECTOR_ELT(rnode,0)));
    double w = *((double *)REAL(VECTOR_ELT(rnode,1)));
    std::vector<int> *vars = new std::vector<int>;
    std::vector<int> *children = new std::vector<int>;
    SEXP rvars = VECTOR_ELT(rnode,2);
    SEXP rchildren = VECTOR_ELT(rnode,3);
    int *pvars = (int *) INTEGER(rvars);
    int *pchildren = (int *) INTEGER(rchildren);
    for(int i = 0;i < Rf_length(rvars);i++)
      vars->push_back(*pvars++);
    for(int i = 0;i < Rf_length(rchildren);i++)
      children->push_back(*pchildren++);
    StructNodeElem<DATA_TYPE> *node = new StructNodeElem<DATA_TYPE>(inode,w,vars,children);
    $1->push_back(node);
  }
}

%typemap(freearg)
  (std::vector<StructNodeElem<DATA_TYPE> *> *TREE)
{
  del_gstruct(arg$argnum);
}
%enddef /* %node_typemaps */

/* ********************** */

%vector_typemaps(REALSXP,REAL,float)
%vector_typemaps(REALSXP,REAL,double)
%vector_typemaps(INTSXP,INTEGER,int)

%matrix_typemaps(REALSXP,REAL,float)
%matrix_typemaps(REALSXP,REAL,double)
#%matrix_typemaps(LGLSXP,LOGICAL,bool)

%spmatrix_typemaps(REALSXP,REAL,float)
%spmatrix_typemaps(REALSXP,REAL,double)
%spmatrix_typemaps(LGLSXP,LOGICAL,bool)

%dspmatrix_typemaps(REALSXP,REAL,float)
%dspmatrix_typemaps(REALSXP,REAL,double)

%datamatrix_typemaps(REALSXP,REAL,float)
%datamatrix_typemaps(REALSXP,REAL,double)

%node_typemaps(REALSXP,REAL,float)
%node_typemaps(REALSXP,REAL,double)

// must instatiate templates with different names, else generated R file will do bad dispatching
%define INSTANTIATE_DATA( f_name )
%template(## f_name) _ ## f_name<double>;
//%template(f_ ## f_name) _ ## f_name<float>;
%enddef



// ARGOUTVIEW_ARRAY1
%typemap(in,numinputs=0) (int** ARGOUTVIEW_ARRAY1,int* DIM1)
  (int *data_temp, int dim_temp)
{
	$1 = &data_temp;
	$2 = &dim_temp;
}
%typemap(argout) (int** ARGOUTVIEW_ARRAY1,int* DIM1) 
{
  SEXP perm;
  if(data_temp$argnum != NULL) {
    int n = dim_temp$argnum;
    PROTECT(perm = Rf_allocVector(INTSXP,n));
    int *pperm = INTEGER(perm);
    memcpy(pperm,data_temp$argnum,n * sizeof(int));
    delete []data_temp$argnum;
    UNPROTECT(1);
  } else perm = R_NilValue;
  MYADD_OUTPUT_ARG($result,perm);
}
