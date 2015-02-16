//#ifdef HAVE_PYTHON
%{
#ifndef SWIG_FILE_WITH_INIT
#  define NO_IMPORT_ARRAY
#endif
#include <stdio.h>
//#include "spams.h"
#undef _POSIX_C_SOURCE
extern "C" {
#include <Python.h>
#include <numpy/arrayobject.h>
}
#define check_array(a,npy_type) (!is_array(a) || !require_contiguous(a) || !require_dimensions(a,1) || !require_native(a) || array_type(a)!=npy_type)
%}

//#endif
%include "numpy.i"

%typemap(throws) const char * %{
  PyErr_SetString(PyExc_RuntimeError, $1);
  SWIG_fail;
%}

%define %vector_typemaps(DATA_TYPE,DATA_TYPECODE)

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros") (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
	$1 = is_array($input) && (array_numdims($input) == 1) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);

}
%typemap(in,fragment="NumPy_Fragments") (Vector<DATA_TYPE> *INPLACE_VECTOR)
	(PyArrayObject* array=NULL)
{
    array = obj_to_array_no_conversion($input, DATA_TYPECODE);
    if (!array || !require_dimensions(array,1) || !require_contiguous(array) || !require_native(array)) SWIG_fail;
	 $1 = new Vector<DATA_TYPE> ((DATA_TYPE *)array_data(array),(int)array_size(array,0));
}
%typemap(out) (Vector<DATA_TYPE> *) 
{
    npy_intp n = result->n();
    npy_intp dims[1] = {n};
    PyArrayObject * array = (PyArrayObject * )PyArray_SimpleNew(1, dims, DATA_TYPECODE);
    DATA_TYPE *data = (DATA_TYPE *)array->data;
    DATA_TYPE *idata = result->rawX();
    memcpy(data,idata,n * sizeof(DATA_TYPE));
    delete result;
    $result = SWIG_Python_AppendOutput($result,(PyObject*)array);
    	
}
%typemap(freearg)
  (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
	delete arg$argnum;
}

// ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Vector<DATA_TYPE> **ARGOUT_VECTOR)
{
	$1 = 1;

}

%typemap(in,numinputs=0,fragment="NumPy_Fragments") (Vector<DATA_TYPE> **ARGOUT_VECTOR)
(Vector<DATA_TYPE>  *data_temp)
{
	# argout in
	$1 = &data_temp;	
}

%typemap(argout) (Vector<DATA_TYPE> **ARGOUT_VECTOR ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    npy_intp n = data_temp$argnum->n();
	    npy_intp dims[1] = {n};
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    PyArrayObject * array = (PyArrayObject * )PyArray_SimpleNewFromData(1, dims, DATA_TYPECODE,(void*)data);
	   if (!array) SWIG_fail;
           $result = SWIG_Python_AppendOutput($result,(PyObject*)array);
	  }
}

%enddef /* %vector_typemaps */

%define map_matrix(DATA_TYPE,DATA_TYPECODE)
	array = obj_to_array_no_conversion($input, DATA_TYPECODE);
	/* !!!!!
	WARNING! bug (?) : the variable name choosen above must not appear
	in the string, otherwise swig will not correctly generate
	final variable names (above name + number)
	*/
	/* we cannot use require_fortran, because it convert a numpy C array to a numpy 
	fortran array by just modifying the strides */
	if (!array || !require_dimensions(array,2) || !array_is_fortran(array) || !require_native(array)) {
	SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"matrix arg $argnum must be a 2d DATA_TYPE Fortran Array"); SWIG_fail;
	}
	$1 = new Matrix<DATA_TYPE> ((DATA_TYPE *)array_data(array),(int)array_size(array,0),(int)array_size(array,1));

%enddef /* map_matrix */

%define %matrix_typemaps(DATA_TYPE,DATA_TYPECODE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros") (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	$1 = is_array($input) && (array_numdims($input) == 2) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);

}
%typemap(in,fragment="NumPy_Fragments") (Matrix<DATA_TYPE> *INPLACE_MATRIX)
	(PyArrayObject* array=NULL)
{
	map_matrix(DATA_TYPE,DATA_TYPECODE)
}
%typemap(out) (Matrix<DATA_TYPE> *) 
{
    npy_intp m = result->m();
    npy_intp n = result->n();
    npy_intp dims[2] = {m,n};
    
    PyArrayObject * array = (PyArrayObject * )PyArray_SimpleNew(2, dims, DATA_TYPECODE);
    DATA_TYPE *data = (DATA_TYPE *)array->data;
    DATA_TYPE *idata = result->rawX();
    memcpy(data,idata,m * n * sizeof(DATA_TYPE));
    delete result;
    if (! require_fortran(array)) {
       SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"Cannot make a fortran out matrix"); SWIG_fail;}
    $result = SWIG_Python_AppendOutput($result,(PyObject*)array);
 	
}
%typemap(freearg)
  (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	delete arg$argnum;
}
// ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
{
	$1 = 1;

}

%typemap(in,numinputs=0,fragment="NumPy_Fragments") (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
(Matrix<DATA_TYPE>  *data_temp)
{
	# argout in
	$1 = &data_temp;	
}

%typemap(argout) (Matrix<DATA_TYPE> **ARGOUT_MATRIX ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    npy_intp m = data_temp$argnum->m();
	    npy_intp n = data_temp$argnum->n();
	    npy_intp dims[2] = {m,n};
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    PyArrayObject * array = (PyArrayObject * )PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE,(void*)data);
	   if (!array) SWIG_fail;
	   if (! require_fortran(array)) {
       	      SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"Cannot make a fortran argout matrix"); SWIG_fail;}

           $result = SWIG_Python_AppendOutput($result,(PyObject*)array);
        }
}


%enddef /* %matrix_typemaps */

%fragment("DSp_Check","header")
{
   const int check_sparse(PyObject* input) {
    return (PyObject_HasAttrString(input, "indptr") &&
        PyObject_HasAttrString(input, "indices") &&
        PyObject_HasAttrString(input, "data") &&
        PyObject_HasAttrString(input, "shape")
        ) ? 1 : 0;
  }
  const int check_matrix(PyObject* input,int data_typecode) {
   return (is_array(input) && (array_numdims(input) == 2) && PyArray_EquivTypenums(array_type(input),data_typecode));
  }
}

%define map_sparse(DATA_TYPE,DATA_TYPECODE)
        sparray = $input;
	if ( !( PyObject_HasAttrString(sparray, "indptr") &&
            PyObject_HasAttrString(sparray, "indices") &&
            PyObject_HasAttrString(sparray, "data") &&
		PyObject_HasAttrString(sparray, "shape"))) {
	  SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"arg $argnum : not a column compressed sparse matrix");
	  return NULL;
	}
	
        /* fetch sparse attributes */
        PyArrayObject* indptr = (PyArrayObject *) PyObject_GetAttrString(sparray, "indptr");
        PyArrayObject* indices = (PyArrayObject *) PyObject_GetAttrString(sparray, "indices");
        PyArrayObject* data = (PyArrayObject *) PyObject_GetAttrString(sparray, "data");
        PyObject* shape = PyObject_GetAttrString(sparray, "shape");

        /* check that types are OK */
	if (check_array(indptr,NPY_INT))
        {
            PyErr_SetString(PyExc_TypeError,"spmatrix arg$argnum: indptr array should be 1d int's");
            return NULL;
        }

	  if check_array(indices,NPY_INT)
        {
            PyErr_SetString(PyExc_TypeError,"spmatrix arg$argnum: indices array should be 1d int's");
            return NULL;
        }

	    if check_array(data, DATA_TYPECODE)
        {
            PyErr_SetString(PyExc_TypeError,"spmatrix arg$argnum: data array should be 1d and match datatype");
            return NULL;
        }

        if (!PyTuple_Check(shape))
        {
            PyErr_SetString(PyExc_TypeError,"shape should be a tuple");
            return NULL;
        }

        /* get array dimensions */
        int32_t m =PyInt_AsLong(PyTuple_GetItem(shape, 0));
        int32_t n =PyInt_AsLong(PyTuple_GetItem(shape, 1));


	int *pB = (int *)array_data(indptr);
	int *pE = pB + 1;
	int nzmax = (int)array_size(data,0);
	Py_DECREF(indptr);
        Py_DECREF(indices);
        Py_DECREF(data);
        Py_DECREF(shape);

 
	$1 = new SpMatrix<DATA_TYPE> ((DATA_TYPE *)array_data(data),(int *)array_data(indices),pB,pE,m,n,nzmax);
%enddef /* map_sparse */

%define %spmatrix_typemaps(DATA_TYPE,DATA_TYPECODE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros",fragment="DSp_Check") (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	$1 = check_sparse($input);
/*
	$1 = ( PyObject_HasAttrString($input, "indptr") &&
        PyObject_HasAttrString($input, "indices") &&
        PyObject_HasAttrString($input, "data") &&
        PyObject_HasAttrString($input, "shape")
        ) ? 1 : 0;
*/
}
%typemap(in,fragment="NumPy_Fragments") (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
	(PyObject* sparray=NULL)
{
	    /* a column compressed storage sparse matrix in python scipy
       looks like this

       A = csc_matrix( ... )
       A.indptr # pointer array
       A.indices # indices array
       A.data # nonzero values array
       A.shape # size of matrix

       >>> type(A.indptr)
       <type 'numpy.ndarray'> #int32
       >>> type(A.indices)
       <type 'numpy.ndarray'> #int32
       >>> type(A.data)
       <type 'numpy.ndarray'>
       >>> type(A.shape)
       <type 'tuple'>
     */
     map_sparse(DATA_TYPE,DATA_TYPECODE)

}
%typemap(out) (SpMatrix<DATA_TYPE> *) 
{
    npy_intp m = result->m();
    npy_intp n = result->n();
    npy_intp nzmax = result->nzmax();
    npy_intp dims[2] = {m,n};
    dims[0] = n + 1;
    PyArrayObject *indptr = (PyArrayObject * )PyArray_SimpleNew(1,dims, NPY_INT);
    dims[0] = nzmax;
    PyArrayObject *indices = (PyArrayObject * )PyArray_SimpleNew(1,dims, NPY_INT);
    PyArrayObject *vdata = (PyArrayObject * )PyArray_SimpleNew(1,dims, DATA_TYPECODE);
    int i;
    DATA_TYPE *xdata = result->v();
    DATA_TYPE *data = (DATA_TYPE *)array_data(vdata);
    memcpy(data,xdata,nzmax * sizeof(DATA_TYPE));
    npy_int *pi = (npy_int *)array_data(indices);
    int *r = result->r();
    int *pB = result->pB();
    if(sizeof(npy_int) == sizeof(int)) {
       memcpy(pi,r,nzmax * sizeof(int));
       pi = (npy_int *)array_data(indptr);
       memcpy(pi,pB,(n + 1) * sizeof(int));
    } else {
      for(i = 0;i< nzmax;i++) 
    	  *(pi+i) = (npy_int) *(r+i);
      pi = (npy_int *)array_data(indptr);
      for(i = 0;i< n + 1;i++) 
    	  *(pi+i) = (npy_int) *(pB+i);
    }
    PyObject* tuple = PyTuple_New(4);
    PyObject* shape = PyTuple_New(2);
    PyTuple_SetItem(shape, 0,  PyInt_FromLong((long)m));
    PyTuple_SetItem(shape, 1,  PyInt_FromLong((long)n));
    PyTuple_SetItem(tuple,0, (PyObject* )indptr);
    PyTuple_SetItem(tuple,1,(PyObject* )indices);
    PyTuple_SetItem(tuple,2,(PyObject* )vdata);
    PyTuple_SetItem(tuple,3,shape);
    delete result;
    $result = SWIG_Python_AppendOutput($result,tuple);
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
# test argout
  if(data_temp$argnum != NULL) {
    npy_intp m = data_temp$argnum->m();
    npy_intp n = data_temp$argnum->n();
    npy_intp nzmax = data_temp$argnum->nzmax();
    npy_intp dims[2] = {m,n};
    dims[0] = n + 1;
    PyArrayObject *indptr = (PyArrayObject * )PyArray_SimpleNew(1,dims, NPY_INT);
    dims[0] = nzmax;
    PyArrayObject *indices = (PyArrayObject * )PyArray_SimpleNew(1,dims, NPY_INT);
    PyArrayObject *vdata = (PyArrayObject * )PyArray_SimpleNew(1,dims, DATA_TYPECODE);
    if (! indptr || !indices || !vdata) SWIG_fail;
    int i;
    DATA_TYPE *xdata = data_temp$argnum->v();
    DATA_TYPE *data = (DATA_TYPE *)array_data(vdata);
    memcpy(data,xdata,nzmax * sizeof(DATA_TYPE));
    npy_int *pi = (npy_int *)array_data(indices);
    int *r = data_temp$argnum->r();
    int *pB = data_temp$argnum->pB();
    if(sizeof(npy_int) == sizeof(int)) {
       memcpy(pi,r,nzmax * sizeof(int));
       pi = (npy_int *)array_data(indptr);
       memcpy(pi,pB,(n + 1) * sizeof(int));
    } else {
      for(i = 0;i< nzmax;i++) 
    	  *(pi+i) = (npy_int) *(r+i);
      pi = (npy_int *)array_data(indptr);
      for(i = 0;i< n + 1;i++) 
    	  *(pi+i) = (npy_int) *(pB+i);
    }
    PyObject* tuple = PyTuple_New(4);
    PyObject* shape = PyTuple_New(2);
    PyTuple_SetItem(shape, 0,  PyInt_FromLong((long)m));
    PyTuple_SetItem(shape, 1,  PyInt_FromLong((long)n));
    PyTuple_SetItem(tuple,0, (PyObject* )indptr);
    PyTuple_SetItem(tuple,1,(PyObject* )indices);
    PyTuple_SetItem(tuple,2,(PyObject* )vdata);
    PyTuple_SetItem(tuple,3,shape);
    $result = SWIG_Python_AppendOutput($result,tuple);
  }
}
%enddef /* %spmatrix_typemaps */

%define %dspmatrix_typemaps(DATA_TYPE,DATA_TYPECODE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros",fragment="DSp_Check") (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	if( PyObject_HasAttrString($input, "indptr")) 
	    $1 = check_sparse($input);
	else
	    $1 = check_matrix($input,DATA_TYPECODE);
#if 0
	$1 = ( PyObject_HasAttrString($input, "indptr") &&
        PyObject_HasAttrString($input, "indices") &&
        PyObject_HasAttrString($input, "data") &&
        PyObject_HasAttrString($input, "shape")
        ) ? 1 : 0;
#endif
}
%typemap(in,fragment="NumPy_Fragments") (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)

{
	    /* a column compressed storage sparse matrix in python scipy
       looks like this

       A = csc_matrix( ... )
       A.indptr # pointer array
       A.indices # indices array
       A.data # nonzero values array
       A.shape # size of matrix

       >>> type(A.indptr)
       <type 'numpy.ndarray'> #int32
       >>> type(A.indices)
       <type 'numpy.ndarray'> #int32
       >>> type(A.data)
       <type 'numpy.ndarray'>
       >>> type(A.shape)
       <type 'tuple'>
     */
     if ( PyObject_HasAttrString($input, "indptr")) {
       	PyObject* sparray =$input;
	map_sparse(DATA_TYPE,DATA_TYPECODE)
     } else {
       	PyArrayObject* array = NULL;
        map_matrix(DATA_TYPE,DATA_TYPECODE)
     }
}
%typemap(freearg)
  (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	delete arg$argnum;
}

%enddef /* %dspmatrix_typemaps */

%define %datamatrix_typemaps(DATA_TYPE,DATA_TYPECODE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros",fragment="DSp_Check") (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	if( PyObject_HasAttrString($input, "indptr")) 
	    $1 = check_sparse($input);
	else
	    $1 = check_matrix($input,DATA_TYPECODE);
}
%typemap(in,fragment="NumPy_Fragments") (Data<DATA_TYPE> *INPLACE_DATAMATRIX)

{
     if ( PyObject_HasAttrString($input, "indptr")) {
       	PyObject* sparray =$input;
	map_sparse(DATA_TYPE,DATA_TYPECODE)
     } else {
       	PyArrayObject* array = NULL;
        map_matrix(DATA_TYPE,DATA_TYPECODE)
     }
}
%typemap(freearg)
  (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	delete arg$argnum;
}

%enddef /* %datamatrix_typemaps */

%define %node_typemaps(DATA_TYPE,DATA_TYPECODE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros") (std::vector<StructNodeElem<DATA_TYPE> *> *TREE)
{
	$1 = PyList_Check($input);

}

%typemap(in) (std::vector<StructNodeElem<DATA_TYPE> *> *TREE) 
{
  PyObject* pytree = $input;
  if(!PyList_Check(pytree)) {
    SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"arg $argnum must be a list");SWIG_fail;
  }
  $1 = new std::vector<StructNodeElem<DATA_TYPE> *>;
  for(Py_ssize_t i = 0;i < PyList_Size(pytree);i++) {
    PyObject* pynode = PyList_GetItem(pytree,i);
    if(! PyTuple_Check(pynode) || (PyTuple_Size(pynode) != 4)) {
      SWIG_Python_SetErrorMsg(PyExc_RuntimeError,"List elements of arg $argnum must be tuples of size 4");SWIG_fail;
    }
    long inode = PyInt_AsLong(PyTuple_GetItem(pynode,(Py_ssize_t)0));
    DATA_TYPE w = PyFloat_AsDouble(PyTuple_GetItem(pynode,(Py_ssize_t)1));
    std::vector<int> *vars = new std::vector<int>;
    std::vector<int> *children = new std::vector<int>;
    PyObject* pyvars = PyTuple_GetItem(pynode,(Py_ssize_t)2);
    PyObject* pychildren = PyTuple_GetItem(pynode,(Py_ssize_t)3);
    for(Py_ssize_t j = 0;j < PyList_Size(pyvars);j++)
      vars->push_back(static_cast<int>(PyInt_AsLong(PyList_GetItem(pyvars,j))));
    for(Py_ssize_t j = 0;j < PyList_Size(pychildren);j++)
      children->push_back(static_cast<int>(PyInt_AsLong(PyList_GetItem(pychildren,j))));
    StructNodeElem<DATA_TYPE> *node = new StructNodeElem<DATA_TYPE>(inode,w,vars,children);
    $1->push_back(node);
  }
    
}

%typemap(out) (std::vector<StructNodeElem<DATA_TYPE> *> *)
{	      
  //int n = result->size();
  PyObject* node_list = PyList_New(0);
  for(std::vector<StructNodeElem<DATA_TYPE> *>::iterator it = result->begin();it != result->end();it++) {
    PyObject* tuple = PyTuple_New(4);
    StructNodeElem<DATA_TYPE> *node = *it;
    int inode = node->node_num;
    PyTuple_SetItem(tuple,0, PyInt_FromLong((long)inode));
    PyTuple_SetItem(tuple,1, PyFloat_FromDouble(node->weight));
    int k = node->vars->size();
    PyObject *vars = PyList_New(0);
    std::vector<int> *pvars = node->vars;
    for(int i = 0;i < k;i++)
      PyList_Append(vars,PyInt_FromLong((long)(*pvars)[i]));
    PyTuple_SetItem(tuple,2, (PyObject* )vars);
    k = node->children->size();
    pvars = node->children;
    PyObject *children = PyList_New(0);
    for(int i = 0;i < k;i++)
      PyList_Append(children,PyInt_FromLong((long)(*pvars)[i]));
    
    PyTuple_SetItem(tuple,3,(PyObject* )children );
    PyList_Append(node_list,tuple);
  }
  del_gstruct(result);
  $result = SWIG_Python_AppendOutput($result,node_list);
}

%typemap(freearg)
  (std::vector<StructNodeElem<DATA_TYPE> *> *TREE)
{
  del_gstruct(arg$argnum);
}
%enddef /* %node_typemaps */

%vector_typemaps(float, NPY_FLOAT)
%vector_typemaps(double, NPY_DOUBLE)
%vector_typemaps(int, NPY_INT)

%matrix_typemaps(float, NPY_FLOAT)
%matrix_typemaps(double, NPY_DOUBLE)
%matrix_typemaps(bool, NPY_BOOL)

%spmatrix_typemaps(float, NPY_FLOAT)
%spmatrix_typemaps(double, NPY_DOUBLE)
%spmatrix_typemaps(bool, NPY_BOOL)

%dspmatrix_typemaps(float, NPY_FLOAT)
%dspmatrix_typemaps(double, NPY_DOUBLE)

%datamatrix_typemaps(float, NPY_FLOAT)
%datamatrix_typemaps(double, NPY_DOUBLE)

%node_typemaps(float, NPY_FLOAT)
%node_typemaps(double, NPY_DOUBLE)

// In case of multiple instantiation :
// template with same name : OK in spite of warning
// But the args are checked to choose between implementations
// special type need a typecheck typemap
#ifndef WINDOWS
       %define INSTANTIATE_DATA( f_name )
       %feature("autodoc","1") _ ## f_name;
       %template(f_name) _ ## f_name<double>;
       %template(f_name) _ ## f_name<float>;
       %enddef
#else
	%define INSTANTIATE_DATA( f_name )
	%feature("autodoc","1") _ ## f_name;
	%template(f_name) _ ## f_name<double>;
	%enddef
#endif

