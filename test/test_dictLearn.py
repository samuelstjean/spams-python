import os.path
import numpy as np
import scipy.sparse as ssp

try:
    from PIL import Image
except Exception as e:
    print("No module PIL.\nYou need to install it if you want to run TrainDL tests\n")
    raise e

import spams
import time

ssprand = ssp.rand
imgpath = os.path.dirname(os.path.realpath(__file__))
import pytest

from scipy.optimize import nnls

def _extract_lasso_param(f_param):
    lst = [ 'L','lambda1','lambda2','mode','pos','ols','numThreads','length_path','verbose','cholesky']
    l_param = {'return_reg_path' : False}
    for x in lst:
        if x in f_param:
            l_param[x] = f_param[x]
    return l_param

def _objective(X,D,param,imgname = None):
    print('Evaluating cost function...')
    lparam = _extract_lasso_param(param)
    alpha = spams.lasso(X,D = D,**lparam)
    # NB : as alpha is sparse, D*alpha is the dot product
    xd = X - D * alpha
    R = np.mean(0.5 * (xd * xd).sum(axis=0) + param['lambda1'] * np.abs(alpha).sum(axis=0))
    print("objective function: %f" %R)
    #* display ?
    if imgname is not None:
        img = spams.displayPatches(D)
        print("IMG %s" %str(img.shape))
        x = np.uint8(img[:,:,0] * 255.)
        image = Image.fromarray(x,mode = 'L')
        image.save("%s.png" %imgname)

@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_trainDL(myfloat):
    img_file = os.path.join(imgpath, 'boat.png')
    try:
        img = Image.open(img_file)
    except:
        print("Cannot load image %s : skipping test" %img_file)
        return None
    I = np.array(img) / 255
    if I.ndim == 3:
        A = np.asfortranarray(I.reshape((I.shape[0],I.shape[1] * I.shape[2])))
        rgb = True
    else:
        A = np.asfortranarray(I)
        rgb = False

    m = 8;n = 8;
    X = spams.im2col_sliding(A,m,n,rgb)

    X = X - np.tile(np.mean(X,0),(X.shape[0],1))
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)),dtype = myfloat)
    param = { 'K' : 100, # learns a dictionary with 100 elements
              'lambda1' : 0.15, 'numThreads' : 4, 'batchsize' : 400,
              'iter' : 1000}

    ########## FIRST EXPERIMENT ###########
    tic = time.time()
    D = spams.trainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    ##param['approx'] = 0
    # save dictionnary as dict.png
    _objective(X,D,param,'dict')

    #### SECOND EXPERIMENT ####
    print("*********** SECOND EXPERIMENT ***********")

    X1 = X[:,0:X.shape[1]//2]
    X2 = X[:,X.shape[1]//2 -1:]
    param['iter'] = 500
    tic = time.time()
    (D,model) = spams.trainDL(X1,return_model = True,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f\n' %t)

    _objective(X,D,param,'dict1')

    # Then reuse the learned model to retrain a few iterations more.
    param2 = param.copy()
    param2['D'] = D
    tic = time.time()
    (D,model) = spams.trainDL(X2,return_model = True,model = model,**param2)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)
    _objective(X,D,param,'dict2')

    #################### THIRD & FOURTH EXPERIMENT ######################
    # let us add sparsity to the dictionary itself

    print('*********** THIRD EXPERIMENT ***********')
    param['modeParam'] = 0
    param['iter'] = 1000
    param['gamma1'] = 0.3
    param['modeD'] = 1

    tic = time.time()
    D = spams.trainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)
    _objective(X,D,param)

    #* DISPLAY
    print('*********** FOURTH EXPERIMENT ***********')
    param['modeParam'] = 0
    param['iter'] = 1000
    param['gamma1'] = 0.3
    param['modeD'] = 3

    tic = time.time()
    D = spams.trainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)
    _objective(X,D,param)

    return None

@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_trainDL_Memory(myfloat):
    img_file = os.path.join(imgpath, 'lena.png')
    try:
        img = Image.open(img_file)
    except:
        print("Cannot load image %s : skipping test" %img_file)
        return None
    I = np.array(img) / 255
    if I.ndim == 3:
        A = np.asfortranarray(I.reshape((I.shape[0],I.shape[1] * I.shape[2])))
        rgb = True
    else:
        A = np.asfortranarray(I)
        rgb = False

    m = 8;n = 8;
    X = spams.im2col_sliding(A,m,n,rgb)

    X = X - np.tile(np.mean(X,0),(X.shape[0],1))
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)))
    X = np.asfortranarray(X[:,np.arange(0,X.shape[1],10)],dtype = myfloat)

    param = { 'K' : 200, # learns a dictionary with 100 elements
          'lambda1' : 0.15, 'numThreads' : 4,
          'iter' : 100}

    ############# FIRST EXPERIMENT  ##################
    tic = time.time()
    D = spams.trainDL_Memory(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    print('Evaluating cost function...')
    lparam = _extract_lasso_param(param)
    alpha = spams.lasso(X,D = D,**lparam)
    xd = X - D * alpha
    R = np.mean(0.5 * (xd * xd).sum(axis=0) + param['lambda1'] * np.abs(alpha).sum(axis=0))
    print("objective function: %f" %R)
    #* ? DISPLAY

    ############# SECOND EXPERIMENT  ##################
    tic = time.time()
    D = spams.trainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)
    print('Evaluating cost function...')
    alpha = spams.lasso(X,D = D,**lparam)
    xd = X - D * alpha
    R = np.mean(0.5 * (xd * xd).sum(axis=0) + param['lambda1'] * np.abs(alpha).sum(axis=0))
    print("objective function: %f" %R)

    #* ? DISPLAY

    return None

@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_structTrainDL(myfloat):
    img_file = os.path.join(imgpath, 'lena.png')
    try:
        img = Image.open(img_file)
    except Exception as e:
        print("Cannot load image %s (%s) : skipping test" %(img_file,e))
        return None
    I = np.array(img) / 255
    if I.ndim == 3:
        A = np.asfortranarray(I.reshape((I.shape[0],I.shape[1] * I.shape[2])),dtype = myfloat)
        rgb = True
    else:
        A = np.asfortranarray(I,dtype = myfloat)
        rgb = False

    m = 8;n = 8;
    X = spams.im2col_sliding(A,m,n,rgb)

    X = X - np.tile(np.mean(X,0),(X.shape[0],1))
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)),dtype = myfloat)
    param = { 'K' : 64, # learns a dictionary with 100 elements
              'lambda1' : 0.05, 'tol' : 1e-3,
              'numThreads' : 4, 'batchsize' : 400,
              'iter' : 20}
    paramL = {'lambda1' : 0.05, 'numThreads' : 4}

    param['regul'] = 'l1'
    print("with Fista Regression %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)
    _objective(X,D,param)

#
    param['regul'] = 'l2'
    print("with Fista Regression %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)

#
    param['regul'] = 'elastic-net'
    print("with Fista %s" %param['regul'])
    param['lambda2'] = 0.1
    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)
## if we want a pause :
##    s = raw_input("graph> ")
########### GRAPH
    param['lambda1'] = 0.1
    param['tol'] = 1e-5
    param['K'] = 10
    eta_g = np.array([1, 1, 1, 1, 1],dtype=myfloat)

    groups = ssp.csc_matrix(np.array([[0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 1, 0, 0]],dtype=bool),dtype=bool)

    groups_var = ssp.csc_matrix(np.array([[1, 0, 0, 0, 0],
                           [1, 0, 0, 0, 0],
                           [1, 0, 0, 0, 0],
                           [1, 1, 0, 0, 0],
                           [0, 1, 0, 1, 0],
                           [0, 1, 0, 1, 0],
                           [0, 1, 0, 0, 1],
                           [0, 0, 0, 0, 1],
                           [0, 0, 0, 0, 1],
                           [0, 0, 1, 0, 0]],dtype=bool),dtype=bool)

    graph = {'eta_g': eta_g,'groups' : groups,'groups_var' : groups_var}
    param['graph'] = graph
    param['tree'] = None

    param['regul'] = 'graph'
    print("with Fista %s" %param['regul'])

    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)

    param['regul'] = 'graph-ridge'
    print("with Fista %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)
## if we want a pause :
##    s = raw_input("tree> ")
##### TREE
    tree_data = """0 1. [] -> 1 4
1 1. [0 1 2] -> 2 3
4 2. [] -> 5 6
2 1. [3 4]
3 2. [5]
5 2. [6 7]
6 2.5 [8] -> 7
7 2.5 [9]
"""
    param['lambda1'] = 0.001
    param['tol'] = 1e-5
    own_variables =  np.array([0,0,3,5,6,6,8,9],dtype=np.int32)
    N_own_variables =  np.array([0,3,2,1,0,2,1,1],dtype=np.int32)
    eta_g = np.array([1,1,1,2,2,2,2.5,2.5],dtype=myfloat)
    groups = np.asfortranarray([[0, 0, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 1, 0]],dtype = bool)
    groups = ssp.csc_matrix(groups,dtype=bool)
    tree = {'eta_g': eta_g,'groups' : groups,'own_variables' : own_variables,
            'N_own_variables' : N_own_variables}
    param['tree'] = tree
    param['graph'] = None
    param['regul'] = 'tree-l0'
    print("with Fista %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X,**param)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)

    gstruct = spams.groupStructOfString(tree_data)
    (perm,tree,nbvars) = spams.treeOfGroupStruct(gstruct)
    param['tree'] = tree

    param['regul'] = 'tree-l2'
    print("with Fista %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X.astype(np.float64),**param).astype(myfloat)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)

    param['regul'] = 'tree-linf'
    print("with Fista %s" %param['regul'])
    tic = time.time()
    D = spams.structTrainDL(X.astype(np.float64),**param).astype(myfloat)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    _objective(X,D,param)


@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_nmf(myfloat):
    img_file = os.path.join(imgpath, 'boat.png')
    try:
        img = Image.open(img_file)
    except:
        print("Cannot load image %s : skipping test" %img_file)
        return None
    I = np.array(img) / 255
    if I.ndim == 3:
        A = np.asfortranarray(I.reshape((I.shape[0],I.shape[1] * I.shape[2])),dtype = myfloat)
        rgb = True
    else:
        A = np.asfortranarray(I,dtype = myfloat)
        rgb = False

    m = 16;n = 16;
    X = spams.im2col_sliding(A,m,n,rgb)
    X = X[:,::10]
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)),dtype = myfloat)
    ########## FIRST EXPERIMENT ###########
    tic = time.time()
    (U,V) = spams.nmf(X,return_lasso= True,K = 49,numThreads=4,iter = -5)
    tac = time.time()
    t = tac - tic
    print('time of computation for Dictionary Learning: %f' %t)

    print('Evaluating cost function...')
    Y = X - U * V
    R = np.mean(0.5 * (Y * Y).sum(axis=0))
    print('objective function: %f' %R)
    return None


# Archetypal Analysis, run first steps with FISTA and run last steps with activeSet,
@pytest.mark.parametrize("myfloat", [np.float64])
def test_archetypalAnalysis(myfloat):
    img_file = os.path.join(imgpath, 'lena.png')
    try:
        img = Image.open(img_file)
    except Exception as e:
        print("Cannot load image %s (%s) : skipping test" %(img_file,e))
        return None
    I = np.array(img) / 255
    if I.ndim == 3:
        A = np.asfortranarray(I.reshape((I.shape[0],I.shape[1] * I.shape[2])),dtype = myfloat)
        rgb = True
    else:
        A = np.asfortranarray(I,dtype = myfloat)
        rgb = False

    m = 8;n = 8;
    X = spams.im2col_sliding(A,m,n,rgb)

    X = X - np.tile(np.mean(X,0),(X.shape[0],1))
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)),dtype = myfloat)
    K = 64 # learns a dictionary with 64 elements
    robust = False # use robust archetypal analysis or not, default parameter(True)
    epsilon = 1e-3 # width in Huber loss, default parameter(1e-3)
    computeXtX = True # memorize the product XtX or not default parameter(True)
    stepsFISTA = 0 # 3 alternations by FISTA, default parameter(3)
    # a for loop in FISTA is used, we stop at 50 iterations
    # remember that we are not guarantee to descent in FISTA step if 50 is too small
    stepsAS = 10 # 7 alternations by activeSet, default parameter(50)
    randominit = True # random initilazation, default parameter(True)

    ############# FIRST EXPERIMENT  ##################
    tic = time.time()
    # learn archetypes using activeSet method for each convex sub-problem
    (Z,A,B) = spams.archetypalAnalysis(np.asfortranarray(X[:, :10000]), returnAB= True, p = K, robust = robust, epsilon = epsilon, computeXtX = computeXtX,  stepsFISTA = stepsFISTA , stepsAS = stepsAS, numThreads = -1)
    tac = time.time()
    t = tac - tic
    print('time of computation for Archetypal Dictionary Learning: %f' %t)

    print('Evaluating cost function...')
    alpha = spams.decompSimplex(np.asfortranarray(X[:, :10000]),Z = Z, computeXtX = True, numThreads = -1)
    xd = X[:,:10000] - Z * alpha
    R = np.sum(xd*xd)
    print("objective function: %f" %R)

    ############# FIRST EXPERIMENT  ##################
    tic = time.time()
    # learn archetypes using activeSet method for each convex sub-problem
    Z2 = spams.archetypalAnalysis(np.asfortranarray(X[:, :10000]), Z0 = Z, robust = robust, epsilon = epsilon, computeXtX = computeXtX , stepsFISTA = stepsFISTA,stepsAS = stepsAS, numThreads = -1)
    tac = time.time()
    t = tac - tic
    print('time of computation for Archetypal Dictionary Learning (Continue): %f' %t)

    print('Evaluating cost function...')
    alpha = spams.decompSimplex(np.asfortranarray(X[:, :10000]),Z = np.asfortranarray(Z2), computeXtX = True, numThreads = -1)
    xd = X[:,:10000] - Z2 * alpha
    R = np.sum(xd*xd)
    print("objective function: %f" %R)

    # learn archetypes using activeSet method for each convex sub-problem
    (Z3,A3,B3) = spams.archetypalAnalysis(np.asfortranarray(X[:, :10000]), returnAB= True, p = K, robust = True, epsilon = epsilon, computeXtX = computeXtX,  stepsFISTA = stepsFISTA , stepsAS = stepsAS, numThreads = -1)
    tac = time.time()
    t = tac - tic
    print('time of computation for Robust Archetypal Dictionary Learning: %f' %t)

@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_lasso_weighted_pos(myfloat):
    rng = np.random.default_rng(123456)
    m, n, p = 50, 100, 200

    alpha = rng.standard_normal([n, p])
    D = rng.standard_normal([m, n])
    D = np.abs(D) * 10
    X = D @ alpha

    X = np.asfortranarray(X).astype(myfloat)
    D = np.asfortranarray(D).astype(myfloat)
    W = np.ones(alpha.shape, dtype=myfloat, order='F')

    alpha1 = spams.lasso(X, D, lambda1=0.0, pos=True).toarray()
    alpha2 = spams.lassoWeighted(X, D, W=W, lambda1=0.0, pos=True).toarray()
    alpha3 = [nnls(D, X[:, i])[0] for i in range(X.shape[1])]
    alpha3 = np.asarray(alpha3).T

    assert np.all(alpha1 >= 0)
    assert np.all(alpha2 >= 0)
    assert np.all(alpha3 >= 0)

    np.testing.assert_allclose(alpha1, alpha2)

    # only tested on double precision because nnls casts internally to double
    # while spams runs at the given precision
    if myfloat is np.float64:
        np.testing.assert_allclose(alpha2, alpha3)
