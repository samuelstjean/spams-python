import numpy as np
import scipy.sparse as ssp
import spams
import time

from test_utils import *

ssprand = ssp.rand


def test_sort():
    n = 2000000
    X = np.asfortranarray(np.random.normal(size=(n,)), dtype=myfloat)
    locs = locals()
    return Xtest("np.sort(X)", "spams.sort(X,True)", locs)


def test_calcAAt():
    """
    test A * A'
    """
    m = 200
    n = 200000
    d = 0.05
    A = ssprand(m, n, density=d, format="csc", dtype=myfloat)
    locs = locals()
    return Xtest("(A * A.T).todense()", "spams.calcAAt(A)", locs)


def test_calcXAt():
    m = 200
    n = 200000
    d = 0.05
    A = ssprand(m, n, density=d, format="csc", dtype=myfloat)
    X = np.asfortranarray(np.random.normal(size=(64, n)), dtype=myfloat)

    # * dot is very very slow betewwen a full and a sparse matrix
    locs = locals()
    return Xtest("np.dot(X,A.T.todense())", "spams.calcXAt(X,A)", locs)


def test_calcXY():
    X = np.asfortranarray(np.random.normal(size=(64, 200)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(200, 20000)), dtype=myfloat)
    locs = locals()
    return Xtest("np.dot(X,Y)", "spams.calcXY(X,Y)", locs)


def test_calcXYt():
    X = np.asfortranarray(np.random.normal(size=(64, 200)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(20000, 200)), dtype=myfloat)
    locs = locals()
    return Xtest("np.dot(X,Y.T)", "spams.calcXYt(X,Y)", locs)


def test_calcXtY():
    X = np.asfortranarray(np.random.normal(size=(200, 64)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(200, 20000)), dtype=myfloat)
    locs = locals()
    return Xtest("np.dot(X.T,Y)", "spams.calcXtY(X,Y)", locs)


def test_bayer():
    n = 2000000
    X = np.asfortranarray(np.random.normal(size=(n,)), dtype=myfloat)

    locs = locals()
    Z = Xtest1("spams", "spams.bayer(X,0)", locs)
    return None


def test_conjGrad():
    A = np.asfortranarray(np.random.normal(size=(5000, 500)))
    # *    np.random.seed(0)
    # *    A = np.asfortranarray(np.random.normal(size = (10,5)))
    A = np.asfortranarray(np.dot(A.T, A), dtype=myfloat)
    b = np.ones((A.shape[1],), dtype=myfloat, order="F")
    x0 = b
    tol = 1e-4
    itermax = int(0.5 * len(b))

    tic = time.time()
    for i in range(0, 20):
        y1 = np.linalg.solve(A, b)
    tac = time.time()
    print("  Time (numpy): ", tac - tic)
    x1 = np.abs(b - np.dot(A, y1))
    print("Mean error on b : %f" % (x1.sum() / b.shape[0]))

    tic = time.time()
    for i in range(0, 20):
        y2 = spams.conjGrad(A, b, x0, tol, itermax)
    # *        y2 = spams.conjGrad(A,b)
    tac = time.time()
    print("  Time (spams): ", tac - tic)
    x1 = np.dot(A, y2)
    x2 = np.abs(b - x1)
    print("Mean error on b : %f" % (x2.sum() / b.shape[0]))

    err = np.abs(y1 - y2)
    return err.max()


def test_invSym():
    A = np.asfortranarray(np.random.random(size=(1000, 1000)))
    A = np.asfortranarray(np.dot(A.T, A), dtype=myfloat)
    locs = locals()
    return Xtest("np.linalg.inv(A)", "spams.invSym(A)", locs)


def test_normalize():
    A = np.asfortranarray(np.random.random(size=(100, 1000)), dtype=myfloat)
    locs = locals()
    res2 = Xtest1("spams", "spams.normalize(A)", locs)
    return None

tests = [
    'sort' , test_sort,
    'calcAAt' , test_calcAAt,
    'calcXAt' , test_calcXAt,
    'calcXY' , test_calcXY,
    'calcXYt' , test_calcXYt,
    'calcXtY' , test_calcXtY,
    'bayer' , test_bayer,
    'conjGrad' , test_conjGrad,
    'invSym' , test_invSym,
    'normalize' , test_normalize,
    ]
