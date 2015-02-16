if (! (library(png,logical.return= TRUE))) {
  cat("\nNo module png!\nIf you want to run TrainDL tests, you need to install R package png.\n\n")
  LOADED <- FALSE
}
.imagesDir = "../extdata"
.objective <- function(X,D,lambda1,imgname = NULL) {
  .printf("Evaluating cost function...\n")
  alpha = spams.lasso(X,D,return_reg_path = FALSE,lambda1 = lambda1, numThreads = 4)

  R = mean(0.5 * colSums((X - D %*% alpha) ^ 2) + lambda1 * colSums(abs(alpha)))
  .printf("objective function: %f\n",R)
  if(! is.null(imgname)) {
    img = spams.displayPatches(D)
    cat("IMG ",dim(img),"\n")
    .printf("XX %f - %f\n",min(img),max(img))
    writePNG(imgname + '.png')
  }
 
}

test_archetypalAnalysis <- function() {
  I = readPNG(paste(.imagesDir,'boat.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)
  #X = X[,seq(from = 1,to = ncol(X),by = 10)]
  X = X[,1:10000]

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)

########## FIRST EXPERIMENT ###########
  tic = proc.time()
  outputAA <- spams.archetypalAnalysis(X,p = 64,returnAB = TRUE, numThreads = -1, stepsFISTA = 0, stepsAS = 10)
  tac = proc.time()
  Z <- outputAA[[1]]
  A <- outputAA[[2]]
  B <- outputAA[[3]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Archetypal Analysis: %f\n",t)
  alpha <- spams.decompSimplex(X, Z, computeXtX = TRUE, numThreads = -1)
  R = sum(colSums((X - Z %*% alpha) ^ 2))
  .printf("objective function: %f\n",R)
  
  tic = proc.time()
  outputAA2 <- spams.archetypalAnalysis(X,Z0 = Z,returnAB = TRUE, numThreads = -1, stepsFISTA = 0, stepsAS = 10)
  tac = proc.time()
  Z2 <- outputAA[[1]]
  A2 <- outputAA[[2]]
  B2 <- outputAA[[3]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Archetypal Analysis: %f\n",t)
  alpha <- spams.decompSimplex(X, Z2, computeXtX = TRUE, numThreads = -1)
  R = sum(colSums((X - Z %*% alpha) ^ 2))
  .printf("objective function: %f\n",R)

  tic = proc.time()
  outputAA <- spams.archetypalAnalysis(X,p = 64,returnAB = TRUE, robust=TRUE, numThreads = -1, stepsFISTA = 0, stepsAS = 10)
  tac = proc.time()
  Z3 <- outputAA[[1]]
  A3 <- outputAA[[2]]
  B3 <- outputAA[[3]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Robust Archetypal Analysis: %f\n",t)

  return(NULL)
}



test_trainDL <- function() {
  I = readPNG(paste(.imagesDir,'boat.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)

  lambda1 = 0.15

########## FIRST EXPERIMENT ###########
  tic = proc.time()
  D <- spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)

#### SECOND EXPERIMENT ####
  .printf("*********** SECOND EXPERIMENT ***********\n")

  X1 = X[,1:as.integer(ncol(X) / 2 )]
  X2 = X[,as.integer(ncol(X) / 2):ncol(X)]

  tic = proc.time()
  res = spams.trainDL(X1,return_model = TRUE,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 500)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)

                                        # Then reuse the learned model to retrain a few iterations more.

  tic = proc.time()
  res = spams.trainDL(X2,model = model,return_model = TRUE,D = D,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 500)
  tac = proc.time()
  D = res[[1]]
  model = res[[2]]
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)


#################### THIRD & FOURTH EXPERIMENT ######################
                                        # let us add sparsity to the dictionary itself

  .printf('*********** THIRD EXPERIMENT ***********\n')

  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000,modeParam = 0,gamma1 = 0.3,modeD = 'L1L2')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)

  .printf('*********** FOURTH EXPERIMENT ***********\n')

  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4, batchsize = 400,iter = 1000,modeParam = 0,gamma1 = 0.3,modeD = 'L1L2MU')
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)

  return(NULL)
}

#


test_trainDL_Memory <- function() {
  I = readPNG(paste(.imagesDir,'lena.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
 #*#!!   X = X[:,np.arange(0,X.shape[1],10)]
  X = X[,seq(from = 1,to = ncol(X),by = 10)]
  lambda1 = 0.15
  
  ############# FIRST EXPERIMENT  ##################
  tic = proc.time()
  D = spams.trainDL_Memory(X,K = 100,lambda1 = lambda1, numThreads = 4,iter = 100)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)

#### SECOND EXPERIMENT ####
  tic = proc.time()
  D = spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 4,iter = 100)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)

  .objective(X,D,lambda1)

  return(NULL)
}

#####
test_structTrainDL <- function() {
  I = readPNG(paste(.imagesDir,'boat.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }

  m = 8;n = 8;
  X = spams.im2col_sliding(A,m,n)

  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)

  lambda1 = 0.05
  K = 64
  tol = 1e-3
  iter = 20
########## FIRST EXPERIMENT ###########
  regul = 'l1'
  .printf("with Fista Regression %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
#
  regul = 'l2'
  .printf("with Fista Regression %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
#
  regul = 'elastic-net'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)

   # pause :
   # readline("pause ")
   ########### GRAPH
  lambda1 = 0.1
  tol = 1e-5
  K = 10
  eta_g = as.vector(c(1, 1, 1, 1, 1),mode='double')
  groups = as(matrix(as.vector(c(0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 1, 0, 0),mode='logical'),ncol = 5,byrow = T),'CsparseMatrix')

  groups_var = as(matrix(as.vector(c(1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 1, 0, 0, 0,
    0, 1, 0, 1, 0,
    0, 1, 0, 1, 0,
    0, 1, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 1, 0, 0),mode='logical'),ncol = 5,byrow = T),'CsparseMatrix')
  graph = list('eta_g'= eta_g,'groups' = groups,'groups_var' = groups_var)
  
  regul = 'graph'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul, graph = graph)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
#
  regul = 'graph-ridge'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul, graph = graph)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
   # pause :
   # readline("pause ")
   ###########   TREE

  lambda1 = 0.001
  tol = 1e-5
  treedata = "0 1. [] -> 1 4
1 1. [0 1 2] -> 2 3
4 2. [] -> 5 6
2 1. [3 4]
3 2. [5]
5 2. [6 7]
6 2.5 [8] -> 7
7 2.5 [9]
"
  own_variables = as.vector(c(0, 0, 3, 5, 6, 6, 8, 9),mode= 'integer')
  N_own_variables = as.vector(c(0,3,2,1,0,2,1,1),mode= 'integer')
  eta_g = as.vector(c(1,1,1,2,2,2,2.5,2.5),mode = 'double')
  groups = matrix(as.vector(c(0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,
    0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,1,0),mode='logical'),ncol = 8,byrow = T)
  groups = as(groups,'CsparseMatrix')
  tree = list('eta_g'= eta_g,'groups' = groups,'own_variables' = own_variables,
    'N_own_variables' = N_own_variables)

  #
  regul = 'tree-l0'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul, tree = tree)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
#
  gstruct = spams.groupStructOfString(treedata)
  x = spams.treeOfGroupStruct(gstruct)
  tree = x[[2]]
  regul = 'tree-l2'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul, tree = tree)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
#
  regul = 'tree-linf'
  .printf("with Fista  %s\n",regul)
  tic = proc.time()
  D <- spams.structTrainDL(X,K = K,lambda1 = lambda1, tol = tol,numThreads = 4, batchsize = 400,iter = iter,regul = regul, tree = tree)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  .objective(X,D,lambda1)
  
  return(NULL)
}

#
test_nmf <- function() {
  I = readPNG(paste(.imagesDir,'boat.png',sep= '/'))
  if (length(dim(I)) == 3) {
    A = matrix(I,nrow = nrow(I),ncol = 3 * ncol(I))
  } else {
    A = I
  }
  m = 16;n = 16;
  X = spams.im2col_sliding(A,m,n)
  X = X[,seq(from = 1,to = ncol(X),by = 10)]
  X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
  
########## FIRST EXPERIMENT ###########
  tic = proc.time()
  res <- spams.nmf(X,return_lasso= TRUE,K = 49,numThreads=4,iter = -5)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  .printf("time of computation for Dictionary Learning: %f\n",t)
  U = res[[1]]
  V = res[[2]]
  .printf("Evaluating cost function...\n")
  R = mean(0.5 * colSums((X - U %*% V) ^ 2))
  .printf("objective function: %f\n",R)
  return(NULL)
}


test_dictLearn.tests = list('archetypalAnalysis' = test_archetypalAnalysis,
  'trainDL' = test_trainDL,
  'trainDL_Memory' = test_trainDL_Memory,
  'structTrainDL' = test_structTrainDL,
  'nmf' = test_nmf
  )
