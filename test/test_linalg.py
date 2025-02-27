import numpy as np
import scipy.sparse as ssp
import scipy.linalg
import spams
import time

ssprand = ssp.rand

import pytest
from numpy.testing import assert_allclose

@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_sort(myfloat):
    n = 2000000
    X = np.asfortranarray(np.random.normal(size=(n,)), dtype=myfloat)
    return assert_allclose(np.sort(X), spams.sort(X,True))


@pytest.mark.parametrize("myfloat", [np.float64])
def test_calcAAt(myfloat):
    """
    test A * A'
    """
    m = 200
    n = 200000
    d = 0.05
    A = ssprand(m, n, density=d, format="csc", dtype=myfloat)
    return assert_allclose((A @ A.T).todense(), spams.calcAAt(A))


@pytest.mark.parametrize("myfloat", [np.float64])
def test_calcXAt(myfloat):
    m = 200
    n = 200000
    d = 0.05
    A = ssprand(m, n, density=d, format="csc", dtype=myfloat)
    X = np.asfortranarray(np.random.normal(size=(64, n)), dtype=myfloat)
    return assert_allclose((X @ A.T) , spams.calcXAt(X,A))


@pytest.mark.parametrize("myfloat", [np.float64])
def test_calcXY(myfloat):
    X = np.asfortranarray(np.random.normal(size=(64, 200)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(200, 20000)), dtype=myfloat)
    return assert_allclose(np.dot(X,Y), spams.calcXY(X,Y))


@pytest.mark.parametrize("myfloat", [np.float64])
def test_calcXYt(myfloat):
    X = np.asfortranarray(np.random.normal(size=(64, 200)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(20000, 200)), dtype=myfloat)
    return assert_allclose(np.dot(X,Y.T), spams.calcXYt(X,Y))


@pytest.mark.parametrize("myfloat", [np.float64])
def test_calcXtY(myfloat):
    X = np.asfortranarray(np.random.normal(size=(200, 64)), dtype=myfloat)
    Y = np.asfortranarray(np.random.normal(size=(200, 20000)), dtype=myfloat)
    return assert_allclose(np.dot(X.T,Y), spams.calcXtY(X,Y))


@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_bayer(myfloat):
    n = 2000000
    X = np.asfortranarray(np.random.normal(size=(n,)), dtype=myfloat)

    Z = spams.bayer(X,0)
    return None


@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_conjGrad(myfloat):
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
    print("Max error", err.max())


@pytest.mark.parametrize("myfloat", [np.float64])
def test_invSym(myfloat):
    A = np.asfortranarray(np.random.random(size=(1000, 1000)))
    A = np.asfortranarray(np.dot(A.T, A), dtype=myfloat)
    assert_allclose(spams.invSym(A), scipy.linalg.pinvh(A), atol=1e-5)


@pytest.mark.parametrize("myfloat", [np.float32, np.float64])
def test_normalize(myfloat):
    A = np.asfortranarray(np.random.random(size=(100, 1000)), dtype=myfloat)
    res2 = spams.normalize(A)
    return None
