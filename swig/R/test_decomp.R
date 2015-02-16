
test_sparseProject <- function () {
  set.seed(0)
  X = matrix(rnorm(20000 * 100),nrow = 20000,ncol = 100,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  .printf( "\n  Projection on the l1 ball\n")
  tic = proc.time()
  X1 = spams.sparseProject(X,numThreads = -1, pos = FALSE,mode= 1,thrs = 2)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("  Time : %f\n", t)
  if (t != 0)
    .printf("%f signals of size %d projected per second\n",(ncol(X) / t),nrow(X))
  s = colSums(abs(X1))
  .printf("Checking constraint: %f, %f\n",min(s),max(s))

  .printf("\n  Projection on the Elastic-Net\n")
  lambda1 = 0.15
  tic = proc.time()
  X1 = spams.sparseProject(X,numThreads = -1, pos = FALSE,mode= 2,thrs = 2,lambda1= lambda1)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("  Time : %f\n", t)
  if (t != 0)
    .printf("%f signals of size %d projected per second\n",(ncol(X) / t),nrow(X))
  constraints = colSums(X1*X1) + lambda1 * colSums(abs(X1))
  .printf("Checking constraint: %f, %f (Projection is approximate : stops at a kink)\n",min(constraints),max(constraints))

  .printf("\n  Projection on the FLSA\n")
  lambda1 = 0.7
  lambda2 = 0.7
  lambda3 = 0.7
  X = matrix(runif(2000 * 100,0,1),nrow = 2000,ncol = 100,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
 #* matlab : X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
  tic = proc.time()
  X1 = spams.sparseProject(X,numThreads = -1, pos = FALSE,mode= 6,thrs = 2,lambda1= lambda1,lambda2= lambda2,lambda3= lambda3)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("  Time : %f\n", t)
  if (t != 0)
    .printf("%f  signals of size %d projected per second\n",ncol(X) / t,nrow(X))
#  constraints = 0.5 * param[['lambda3']] * colSums(X1*X1) + param[['lambda1']] * colSums(abs(X1)) + param[['lambda2']] * abs(X1[3:nrow(X1),]) - colSums(X1[2,])
  constraints = 0.5 * lambda3 * colSums(X1*X1) + lambda1 * colSums(abs(X1)) + lambda2 * abs(X1[3:nrow(X1),])  - colSums(matrix(X1[2,],nrow=1))
  .printf("Checking constraint: %f, %f (Projection is approximate : stops at a kink)\n",min(constraints),max(constraints))
  return(NULL)
}

test_cd <- function() {
  set.seed(0)
  X = matrix(rnorm(64 * 100),nrow = 64,ncol = 100,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  D = matrix(rnorm(64 * 100),nrow = 64,ncol = 100,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  # parameter of the optimization procedure are chosen
  lambda1 = 0.015
  mode = 'PENALTY'
  tic = proc.time()
  alpha = spams.lasso(X,D,lambda1 = lambda1,mode = mode,numThreads = 4)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second for LARS\n", (ncol(X) / t))
  E = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("Objective function for LARS: %g\n",E)
  tol = 0.001
  itermax = 1000
  A0 = as(matrix(c(0),nrow = nrow(alpha),ncol = ncol(alpha)),'CsparseMatrix')
  tic = proc.time()
  alpha2 = spams.cd(X,D,A0,lambda1 = lambda1,mode = mode,tol = tol, itermax = itermax,numThreads = 4)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second for CD\n", (ncol(X) / t))
  E = mean(0.5 * colSums((X - D %*% alpha2) ^ 2) + lambda1 * colSums(abs(alpha2)))
  .printf("Objective function for CD: %g\n",E)
  .printf("With Random Design, CD can be much faster than LARS\n")
  return(NULL)
}

test_l1L2BCD <- function() {
  set.seed(0)
  X = matrix(rnorm(64 * 100),nrow = 64,ncol = 100,byrow = FALSE)
  D = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  ind_groups = as.vector(seq(from = 0,to = ncol(X) - 1,by = 10),mode= 'integer')#indices of the first signals in each group
  itermax = 100
  tol = 1e-3
  mode = 'PENALTY'
  lambda1 = 0.15 # squared norm of the residual should be less than 0.1
  numThreads = -1 # number of processors/cores to use the default choice is -1
                    # and uses all the cores of the machine

  alpha0 = matrix(c(0),nrow = ncol(D), ncol = ncol(X),byrow = FALSE)
  tic = proc.time()
  alpha = spams.l1L2BCD(X,D,alpha0,ind_groups,lambda1 = lambda1,mode = mode,itermax = itermax,tol = tol,numThreads = numThreads)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
  
  return(NULL)
}

test_lasso <- function() {
  set.seed(0)
  .printf("test lasso\n")
##############################################
# Decomposition of a large number of signals
##############################################
# data generation
  X = matrix(rnorm(100 * 100000),nrow = 100,ncol = 100000,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  D = matrix(rnorm(100 * 200),nrow = 100,ncol = 200,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  tic = proc.time()
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY' )
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
########################################
# Regularization path of a single signal 
########################################
  X = matrix(rnorm(64),nrow = 64,ncol = 1,byrow = FALSE)
  D = matrix(rnorm(64 * 10),nrow = 64,ncol = 10,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  res = spams.lasso(X,D,return_reg_path = TRUE,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY' )
  alpha = res[[1]]
  path = res[[2]]
  return(NULL)
}
test_lassoMask <- function() {
  set.seed(0)
  .printf("test lassoMask\n")
##############################################
# Decomposition of a large number of signals
##############################################
# data generation
  X = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  D = matrix(rnorm(100 * 20),nrow = 100,ncol = 20,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  mask = (X > 0) # generating a binary mask
  tic = proc.time()
  alpha = spams.lassoMask(X,D,mask,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
  return(NULL)
}
test_lassoWeighted <- function() {
  set.seed(0)
  .printf("test lasso weighted\n")
##############################################
# Decomposition of a large number of signals
##############################################
# data generation
  X = matrix(rnorm(64 * 10000),nrow = 64,ncol = 10000,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  D = matrix(rnorm(64 * 256),nrow = 64,ncol = 256,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  W = matrix(runif(ncol(D) * ncol(X),0,1),nrow = ncol(D),ncol = ncol(X),byrow = FALSE)
  tic = proc.time()
  alpha = spams.lassoWeighted(X,D,W,lambda1 = 0.15,numThreads = -1,mode = 'PENALTY')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
  
  return(NULL)
}

test_omp <- function() {
  set.seed(0)
  .printf("test omp\n")
  X = matrix(rnorm(64 * 100000),nrow = 64,ncol = 100000,byrow = FALSE)
  D = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  L = 10
  eps =0.1
#*  L = as.vector(c(10),mode='integer')
#*  eps = as.vector(c(0.1),mode='double')
  numThreads = -1

  tic = proc.time()
  alpha = spams.omp(X,D,L=L,eps=eps,return_reg_path = FALSE,numThreads = numThreads)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)

########################################
# Regularization path of a single signal 
########################################
  X = matrix(rnorm(64 * 1),nrow = 64,ncol = 1,byrow = FALSE)
  D = matrix(rnorm(64 * 10),nrow = 64,ncol = 10,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  L = 5
  res = spams.omp(X,D,L = L,eps = eps,return_reg_path = TRUE,numThreads = numThreads)
  alpha = res[[1]]
  path = res[[2]]
  return(NULL)
}
test_ompMask <- function() {
  set.seed(0)
  .printf("test ompMask\n")
  X = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  D = matrix(rnorm(100 * 20),nrow = 100,ncol = 20,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  mask = (X > 0) # generating a binary mask
  L = 20
  eps =0.01
  numThreads = -1

  tic = proc.time()
  alpha = spams.ompMask(X,D,mask,L = L,eps = eps,return_reg_path = FALSE,numThreads = numThreads)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)

}

test_somp <- function() {
  set.seed(0)
  X = matrix(rnorm(64 * 10000),nrow = 64,ncol = 10000,byrow = FALSE)
  D = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
  D = D / matrix(rep(sqrt(colSums(D*D)),nrow(D)),nrow(D),ncol(D),byrow=T)
  ind_groups = as.vector(seq(from = 0,to = 9999,by = 10),mode= 'integer')
  tic = proc.time()
  alpha = spams.somp(X,D,ind_groups,L = 10,eps = 0.1,numThreads=-1)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("%f signals processed per second\n",as.double(ncol(X)) / t)
  return(NULL)
}

#
test_decomp.tests =list( 
  'sparseProject' = test_sparseProject,
  'cd' = test_cd,
  'l1L2BCD' = test_l1L2BCD,
  'lasso' = test_lasso,
  'lassoMask' = test_lassoMask,
  'lassoWeighted' = test_lassoWeighted,
  'omp' = test_omp,
  'ompMask' = test_ompMask,
  'somp' = test_somp
  )
