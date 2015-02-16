import numpy as np
import sys

def main(argv):
    N = 96  # index size
    m = 100
    X = np.random.random_sample(m * N)
    X = X.reshape((N,m))
    y = np.sum(X,0)
    dX = np.divide(X,np.tile(y,(N,1)))
    z = np.sum(dX,0)
    print "%s y %s, z %s\n%s" %(str(dX.shape),str(y.shape),str(z.shape),str(z))
    f = open("textes.txt","w")
    for j in xrange(0,dX.shape[1]):
        for i in xrange(0,dX.shape[0]):
            print >>f, "%f" %dX[i,j]
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
