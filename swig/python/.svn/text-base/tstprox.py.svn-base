import spams
import numpy as np
import sys

param = {'numThreads' : -1,'verbose' : True,
         'lambda1' : 0.1 }
m = 5;n = 10
U = np.asfortranarray(np.random.normal(size = (m,n)))

# test L0
print "\nprox lo"
param['regul'] = 'l0'
param['pos'] = False       # false by default
param['intercept'] = False # false by default
(alpha,X) = spams.proximalFlat(U,True,**param);
print "ALPHA %s" %str(alpha.shape)
print "XX %s" %str(X.shape)
##
X1 = U.reshape(m * n)
f = open('datax','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
##
X1 = X.reshape(n)
f = open('datay','w')
for x in X1:
    print >> f,"%f" %x
f.close()
##
##
X1 = alpha.reshape(m * n)
f = open('dataw','w')
for x in X1:
    print >> f,"%f" %x
f.close()
