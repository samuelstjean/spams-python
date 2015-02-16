import spams
import numpy as np
import sys
import scipy
import scipy.sparse as ssp
import myscipy
ssprand = myscipy.rand

np.random.seed(0)
X = np.asfortranarray(np.random.normal(size = (2000,100)))
# matlab : X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)))
param = {'numThreads' : -1, # number of processors/cores to use (-1 => all cores)
         'pos' : False,
         'mode': 1, # projection on the l1 ball
         'thrs' : 2}
param['mode'] = 2  # projection on the Elastic-Net
param['lambda1'] = 0.15
param['mode'] = 6       # projection on the FLSA
param['lambda1'] = 0.7
param['lambda2'] = 0.7
param['lambda3'] = 1.0
X1 = spams.SparseProject(X,param)
print "PASS1"
X2 = spams.SparseProject(X,param)
print "ERR %f" %abs(X2 - X1).max()
