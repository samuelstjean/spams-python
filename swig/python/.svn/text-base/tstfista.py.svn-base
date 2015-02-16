import spams
import numpy as np
import sys
import scipy
import scipy.sparse as ssp

def tsteq(x,y):
    a = (x == y)
    if a.all():
        return True
    else:
        return False

param = {'numThreads' : 1,'verbose' : False,
         'lambda1' : 0.05, 'it0' : 10, 'max_it' : 200,
         'L0' : 0.1, 'tol' : 1e-3, 'intercept' : False,
         'pos' : False}
np.random.seed(0)
m = 100;n = 200
X = np.asfortranarray(np.random.normal(size = (m,n)))

X = np.asfortranarray(X - np.tile(np.mean(X,0),(X.shape[0],1)))
##
X1 = X.reshape(m * n)
f = open('datax','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
X = spams.normalize(X)
Y = np.asfortranarray(np.random.normal(size = (m,1)))
##
X1 = Y.reshape(m)
f = open('datay','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
Y = np.asfortranarray(Y - np.tile(np.mean(Y,0),(Y.shape[0],1)))
Y = spams.normalize(Y)
W0 = np.zeros((X.shape[1],Y.shape[1]),dtype=np.float64,order="FORTRAN")
param['compute_gram'] = True
param['verbose'] = True
param['loss'] = 'square'
param['regul'] = 'l1'
if False:
    (W, optim_info) = spams.fistaFlat(Y,X,W0,True,**param)
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
param['regul'] = 'group-lasso-l2'
param2=param
param2['groups'] = np.array(np.random.random_integers(1,5,X.shape[1]),dtype = np.int32)
param2['lambda1'] *= 10
(W, optim_info) = spams.fistaFlat(Y,X,W0,True,**param)
exit()

param['ista'] = False
param['subgrad'] = True
param['a'] = 0.1
param['b'] = 1000 # arbitrary parameters
max_it = param['max_it']
it0 = param['it0']
param['max_it'] = 500
param['it0'] = 50
X0 = np.copy(X)
Y0 = np.copy(Y)
(W, optim_info) = spams.fistaFlat(Y,X,W0,param,True)
print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
print "XX %s" %str(W.shape)
##
X1 = W.reshape(n)
f = open('dataw','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##

if not tsteq(X0,X):
    print "X MODIFIE"
if not tsteq(Y0,Y):
    print "Y MODIFIE"
#W00 = np.zeros((X.shape[1],Y.shape[1]),dtype=np.float64,order="FORTRAN")
(W, optim_info) = spams.fistaFlat(Y,X,W0,param,True)
print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
print "ERR %f" %abs(W).max()
Y = np.asfortranarray(np.ceil(5 * np.random.random(size = (100,1000))) - 1)
param['loss'] = 'multi-logistic'
print '\nFISTA + Multi-Class Logistic l1'
nclasses = np.max(Y[:])+1
W0 = np.zeros((X.shape[1],nclasses * Y.shape[1]),dtype=np.float64,order="FORTRAN")
print "X %s, Y %s, W0 %s" %(str(X.shape),str(Y.shape),str(W0.shape))
(W, optim_info) = spams.fistaFlat(Y,X,W0,param,True)

print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
 
