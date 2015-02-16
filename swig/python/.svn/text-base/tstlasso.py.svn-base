import spams
import numpy as np
import sys
import scipy
import scipy.sparse as ssp
import time

myfloat = np.float32

m = 100;n = 10000;nD = 200
#m = 5;n = 10;nD = 5

np.random.seed(0)
X = np.asfortranarray(np.random.normal(size=(m,n)))
# X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)),dtype= myfloat)
D = np.asfortranarray(np.random.normal(size=(m,nD)))
D = np.asfortranarray(D / np.tile(np.sqrt((D*D).sum(axis=0)),(D.shape[0],1)),dtype= myfloat)
# parameter of the optimization procedure are chosen
#param.L=20; # not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
param = {
    'lambda1' : 0.15, # not more than 20 non-zeros coefficients
    'numThreads' : -1, # number of processors/cores to use; the default choice is -1
    # and uses all the cores of the machine
    'mode' : 2}        # penalized formulation
##
X1 = X.reshape(m * n)
f = open('datax','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
##
if False:
    X1 = D.reshape(m * nD)
    f = open('datay','w')
    for x in X1:
        print >> f,"%f" %x
    f.close()
##

tic = time.time()
(A,path) = spams.lasso(X,D = D,return_reg_path = True,**param)
tac = time.time()
print "XX A %s, path %s" %(str(A.dtype),str(path.dtype))
t = tac - tic
print "T= %f,%f signals processed per second\n" %(t,float(X.shape[1]) / t)
##
exit
X1 = np.asfortranarray(A.todense()).reshape(n * nD)
f = open('dataw','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
