library("spams",lib.loc = "./lib")
library('Matrix')

.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}
#m = 5;n = 10;nD = 5
m = 100; n = 100;k = 20
set.seed(0)
#  X = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
#  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
l = scan(file= "../../xxX",what =  double(0))
X = matrix(l,nrow = m,ncol = n,byrow = TRUE)

#  D = matrix(rnorm(100 * 20),nrow = 100,ncol = 20,byrow = FALSE)
#  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
l = scan(file= "../../xxD",what =  double(0))
D = matrix(l,nrow = m,ncol = k,byrow = TRUE)

  mask = (X > 0) # generating a binary mask
  tic = proc.time()
  alpha = spams.lassoMask(X,D,mask,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
#alpha
