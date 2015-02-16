test_fistaFlat <- function() {
  set.seed(0)
  m = 100;n = 200
  X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = spams.normalize(X)
  Y = matrix(rnorm(m),nrow = m,ncol = 1,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  # Regression experiments 
  # 100 regression problems with the same design matrix X.
  .printf("\nVarious regression experiments\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  #*# optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
 ###

  .printf("\nFISTA + Regression l1\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = TRUE)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nSubgradient Descent + Regression l1\n")
  max_it = 200
  it0 = 10
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 50, max_it = 500,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = TRUE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nFISTA + Regression l2\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = 50, max_it = 500,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l2',ista = FALSE,subgrad = TRUE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
  .printf("\nFISTA + Regression l2 + sparse feature matrix\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,as(X,'CsparseMatrix'),W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
#######
  .printf("\nFISTA + Regression Elastic-Net\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'elastic-net',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])

  .printf("\nFISTA + Group Lasso L2\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'group-lasso-l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,size_group = 2)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])

  .printf("\nFISTA + Group Lasso L2 with variable size of groups\n")
  groups = matrix(sample(1:5,ncol(X),replace=TRUE),nrow = 1)
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.5, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'group-lasso-l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,size_group = 2,
    groups = groups)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  .printf("\nFISTA + Trace Norm\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'trace-norm-vec',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  
###
  .printf("\nFISTA + Regression Fused-Lasso\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'fused-lasso',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression no regularization\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'none',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept\n")
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.fistaFlat(Y,x1,W01,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = TRUE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
  .printf("\nFISTA + Regression l1 with intercept+ non-negative\n")
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.fistaFlat(Y,x1,W01,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = TRUE,pos = TRUE,compute_gram = TRUE, loss = 'square',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
  
  .printf("\nISTA + Regression l0\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.05, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'l0',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
# Classification
    
  .printf("\nOne classification experiment\n")
#*    Y = 2 * double(randn(100,1) > 0)-1
  Y = matrix(2. * as.double(rnorm(100) > 0.) - 1.,nrow = 100,ncol = 1,byrow = FALSE)
  .printf("\nFISTA + Logistic l1\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'weighted-logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
#!    pause
    
  .printf("\nFISTA + Logistic l1 + sparse matrix\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,as(X,'CsparseMatrix'),W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...

# Multi-Class classification
#*    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  .printf("\nFISTA + Multi-Class Logistic l1\n")
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'multi-logistic',regul = 'l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
    
# Multi-Task regression
  Y = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  .printf("\nFISTA + Regression l1l2\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  #*# optim_info = vecteur colonne de 4 lignes
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  
  .printf("\nFISTA + Regression l1linf\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
    
    
  .printf("\nFISTA + Regression l1l2 + l1\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1l2+l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + l1\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf+l1',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
    
  .printf("\nFISTA + Regression l1linf + row + columns\n")
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'l1linf-row-column',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])

 # Multi-Task Classification
    
  .printf("\nFISTA + Logistic + l1l2\n")
#*    Y = 2*double(randn(100,100) > 0)-1
  Y = matrix(2. * as.double(rnorm(100 * 100) > 1.) - 1.,nrow = 100,ncol = 100,byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'logistic',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# Multi-Class + Multi-Task Regularization
  
  .printf("\nFISTA + Multi-Class Logistic l1l2\n")
#*    Y = double(ceil(5*rand(100,1000))-1)
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  Y = spams.normalize(Y)
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaFlat(Y,X,W0,TRUE,numThreads = 1,verbose = TRUE,lambda1 = 0.01, it0 = it0, max_it = max_it,L0 = 0.1, tol = 1e-3, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'multi-logistic',regul = 'l1l2',ista = FALSE,subgrad = FALSE,a = 0.1, b = 1000,lambda2 = 0.1,lambda3 = 0.1,size_group = 5)),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...
 
  return(NULL)
}

#
test_fistaGraph <- function() {
  set.seed(0)
  num_threads = -1 # all cores (-1 by default)
  verbose = FALSE   # verbosity, false by default
  lambda1 = 0.1 # regularization ter
  it0 = 1      # frequency for duality gap computations
  max_it = 100 # maximum number of iterations
  L0 = 0.1
  tol = 1e-5
  intercept = FALSE
  pos = FALSE
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
  verbose = TRUE
  X = matrix(rnorm(1000),nrow = 100,ncol = 10,byrow = FALSE)
  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = spams.normalize(X)
  Y = matrix(rnorm(100),nrow = 100,ncol = 1,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
 # Regression experiments 
  # 100 regression problems with the same design matrix X.
  .printf('\nVarious regression experiments\n')
  compute_gram = TRUE
#
  .printf('\nFISTA + Regression graph\n')
  loss = 'square'
  regul = 'graph'
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul)
  
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
  .printf('\nADMM + Regression graph\n')
  admm = TRUE
  lin_admm = TRUE
  c = 1
  delta = 1
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,delta = delta)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
  admm = FALSE
  max_it = 5
  it0 = 1
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,delta = delta)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
#  works also with non graph-structured regularization. graph is ignored
  .printf('\nFISTA + Regression Fused-Lasso\n')
  regul = 'fused-lasso'
  lambda2 = 0.01
  lambda3 = 0.01
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)
  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f,  time: %f, number of iterations: %f\n",optim_info[1],t,optim_info[4])
#
  .printf("\nFISTA + Regression graph with intercept\n")
  regul = 'graph'
  intercept = TRUE
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
  intercept = FALSE
# Classification
  .printf('\nOne classification experiment\n')
  Y = matrix(2. * as.double(rnorm(100 * ncol(Y)) > 0.) - 1.,nrow = 100,ncol = ncol(Y),byrow = FALSE)
  .printf('\nFISTA + Logistic + graph-linf\n')
  loss = 'logistic'
  regul = 'graph'
  lambda1 = 0.01
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
# can be used of course with other regularization functions, intercept,...

# Multi-Class classification
  Y = ceiling(5 * matrix(runif(100 * ncol(Y),0,1),nrow = 100,ncol = ncol(Y),byrow = FALSE)) - 1
  loss = 'multi-logistic'
  regul = 'graph'
  .printf('\nFISTA + Multi-Class Logistic + graph\n')
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
# can be used of course with other regularization functions, intercept,...
  
# Multi-Task regression
  Y = matrix(rnorm(100 * ncol(Y)),nrow = 100,ncol = ncol(Y),byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  compute_gram = FALSE
  verbose = TRUE
  loss = 'square'
  .printf('\nFISTA + Regression  multi-task-graph\n')
  regul = 'multi-task-graph'
  lambda2 = 0.01
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
# Multi-Task Classification
  .printf('\nFISTA + Logistic + multi-task-graph\n')
  regul = 'multi-task-graph'
  lambda2 = 0.01
  loss = 'logistic'
  Y = matrix(rnorm(100 * ncol(Y)),nrow = 100,ncol = ncol(Y),byrow = FALSE)
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
# Multi-Class + Multi-Task Regularization
  verbose = FALSE
  .printf('\nFISTA + Multi-Class Logistic +multi-task-graph\n')
  Y = ceiling(5 * matrix(runif(100 * ncol(Y),0,1),nrow = 100,ncol = ncol(Y),byrow = FALSE)) - 1
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  loss = 'multi-logistic'
  regul = 'multi-task-graph'
  tic = proc.time()
  res = spams.fistaGraph(
    Y,X,W0,graph,TRUE,numThreads = num_threads,verbose = verbose,
    lambda1 = lambda1,it0 = it0,max_it = max_it,L0 = L0,tol = tol,
    intercept = intercept,pos = pos,compute_gram = compute_gram,
    loss = loss,regul = regul,admm = admm,lin_admm = lin_admm,
    c = c,lambda2 = lambda2,lambda3 = lambda3,delta = delta)

  tac = proc.time()
  t = (tac - tic)[['elapsed']]
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n",optim_info[1],optim_info[3],t,optim_info[4])
#
# can be used of course with other regularization functions, intercept,...

  return(NULL)
}

test_fistaTree <- function() {
  set.seed(0)
  m = 100;n = 10
  X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = spams.normalize(X)
  Y = matrix(rnorm(m),nrow = m,ncol = 1,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  own_variables = as.vector(c(0,0,3,5,6,6,8,9),mode= 'integer')
  N_own_variables = as.vector(c(0,3,2,1,0,2,1,1),mode= 'integer')
  eta_g = as.vector(c(1,1,1,2,2,2,2.5,2.5),mode = 'double')
  groups = matrix(as.vector(c(0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0),mode='logical'),ncol = 8,byrow = T)
  groups = as(groups,'CsparseMatrix')
  tree = list('eta_g'= eta_g,'groups' = groups,'own_variables' = own_variables,
            'N_own_variables' = N_own_variables)
  .printf('\nVarious regression experiments\n')
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l2')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
 ###
  .printf('\nFISTA + Regression tree-linf\n')
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-linf')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
# works also with non tree-structured regularization. tree is ignored
  .printf('\nFISTA + Regression Fused-Lasso\n')
   res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'fused-lasso')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
 ###
  .printf('\nISTA + Regression tree-l0\n')
   res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l0')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
 ###
  .printf('\nFISTA + Regression tree-l2 with intercept\n')
  x1 = cbind(X,matrix(1,nrow = nrow(X),ncol = 1))
  W01 = rbind(W0,matrix(0,nrow = 1, ncol = ncol(W0)))
  res = Xtest1('spams',quote(spams.fistaTree(Y,x1,W01,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = TRUE,pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l2')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
 
#    Classification

  .printf('\nOne classification experiment')
  Y = matrix(2. * as.double(rnorm(100 * ncol(Y)) > 0.) - 1.,nrow = 100,ncol = ncol(Y),byrow = FALSE)
  .printf('\nFISTA + Logistic + tree-linf\n')
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'logistic',regul = 'tree-linf')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
###
# can be used of course with other regularization functions, intercept,...

#  Multi-Class classification
  Y = ceiling(5 * matrix(runif(100 * ncol(Y),0,1),nrow = 100,ncol = ncol(Y),byrow = FALSE)) - 1
  .printf('\nFISTA + Multi-Class Logistic + tree-l2\n')
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = TRUE, loss = 'multi-logistic',regul = 'tree-l2')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, number of iterations: %f\n",optim_info[1],optim_info[4])
# can be used of course with other regularization functions, intercept,...

# Multi-Task regression
  Y = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  .printf('\nFISTA + Regression  multi-task-tree\n')
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = TRUE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'square',regul = 'multi-task-tree')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  
# Multi-Task Classification
  .printf('\nFISTA + Logistic + multi-task-tree\n')
  Y = matrix(rnorm(100 * ncol(Y)),nrow = 100,ncol = ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = TRUE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'logistic',regul = 'multi-task-tree')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])

#  Multi-Class + Multi-Task Regularization
  .printf('\nFISTA + Multi-Class Logistic +multi-task-tree\n')
  Y = ceiling(5 * matrix(runif(100 * ncol(Y),0,1),nrow = 100,ncol = ncol(Y),byrow = FALSE)) - 1
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaTree(Y,X,W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'multi-logistic',regul = 'multi-task-tree')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
# can be used of course with other regularization functions, intercept,...

  .printf('\nFISTA + Multi-Class Logistic +multi-task-tree + sparse matrix\n')
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = Xtest1('spams',quote(spams.fistaTree(Y,as(X,'CsparseMatrix'),W0,tree,TRUE,numThreads = -1,verbose = FALSE,lambda1 = 0.001,lambda2 = 0.001,lambda3 = 0.001, it0 = 10, max_it = 200,L0 = 0.1, tol = 1e-5, intercept = FALSE,pos = FALSE,compute_gram = FALSE, loss = 'multi-logistic',regul = 'multi-task-tree')),n = 1)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  

  return(NULL)
}

test_proximalFlat <- function() {
  m = 100;n = 1000
  U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)

  # test l0
  .printf("\nprox l0\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l0', pos = FALSE,intercept = FALSE)),n = 1)
  
  # test l1
  .printf("\nprox l1, intercept, positivity constraint\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1', pos = TRUE,intercept = TRUE)),n = 1)
  
  # test l2
  .printf("\nprox squared-l2\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l2', pos = FALSE,intercept = FALSE)),n = 1)

  # test elastic-net
  .printf("\nprox elastic-net\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'elastic-net', pos = FALSE,intercept = FALSE,lambda2 = 0.1)),n = 1)
  
  # test fused-lasso
  .printf("\nprox fused-lasso\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'fused-lasso', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1l2
  .printf("\nprox mixed norm l1/l2\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1l2', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf
  .printf("\nprox mixed norm l1/linf\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1l2+l1
  .printf("\nprox mixed norm l1/l2 + l1\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1l2+l1', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf+l1
  .printf("\nprox mixed norm l1/linf + l1\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf+l1', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test l1linf-row-column
  .printf("\nprox mixed norm l1/linf on rows and columns\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1linf-row-column', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)
  
  # test none
  .printf("\nprox no regularization\n")
  alpha = Xtest1('spams',quote(spams.proximalFlat(U,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'none', pos = FALSE,intercept = FALSE,lambda2 = 0.1,lambda3 = 0.1)),n = 1)

  return(NULL)
}

test_proximalGraph <- function() {
  set.seed(0)
  lambda1 = 0.1 # regularization parameter
  num_threads = -1 # all cores (-1 by default)
  verbose = TRUE   # verbosity, false by default
  pos = FALSE       # can be used with all the other regularizations
  intercept = FALSE # can be used with all the other regularizations     
  
  U = matrix(rnorm(1000),nrow = 10,ncol = 100,byrow = FALSE)
  .printf('First graph example\n')
# Example 1 of graph structure
# groups:
# g1= {0 1 2 3}
# g2= {3 4 5 6}
# g3= {6 7 8 9}
  eta_g = as.vector(c(1, 1, 1),mode='double')
  groups = as(matrix(as.vector(rep(0,9),mode='logical'),nrow = 3,ncol = 3),'CsparseMatrix')
  groups_var = as(matrix(as.vector(c(1, 0, 0,
    1, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
    0, 1, 0,
    0, 1, 1,
    0, 0, 1,
    0, 0, 1,
    0, 0, 1),mode='logical'),ncol = 3,byrow = T),'CsparseMatrix')
  graph = list('eta_g'= eta_g,'groups' = groups,'groups_var' = groups_var)
  .printf('\ntest prox graph\n')
  regul='graph'
  alpha = Xtest1('spams',quote(spams.proximalGraph(U,graph,FALSE,lambda1 = lambda1,numThreads  = num_threads ,verbose = verbose,pos = pos,intercept = intercept,regul = regul)),n = 1)
  
# Example 2 of graph structure
# groups:
# g1= {0 1 2 3}
# g2= {3 4 5 6}
# g3= {6 7 8 9}
# g4= {0 1 2 3 4 5}
# g5= {6 7 8}
  
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
  
  .printf('\ntest prox graph\n')
  alpha = Xtest1('spams',quote(spams.proximalGraph(U,graph,FALSE,lambda1 = lambda1,numThreads  = num_threads ,verbose = verbose,pos = pos,intercept = intercept,regul = regul)),n = 1)
 #
  .printf('\ntest prox multi-task-graph\n')
  regul = 'multi-task-graph'
  lambda2 = 0.1
  alpha = Xtest1('spams',quote(spams.proximalGraph(U,graph,FALSE,lambda1 = lambda1,lambda2 = lambda2,numThreads  = num_threads ,verbose = verbose,pos = pos,intercept = intercept,regul = regul)),n = 1)
 #
  .printf('\ntest no regularization\n')
  regul = 'none'
  alpha = Xtest1('spams',quote(spams.proximalGraph(U,graph,FALSE,lambda1 = lambda1,lambda2 = lambda2,numThreads  = num_threads ,verbose = verbose,pos = pos,intercept = intercept,regul = regul)),n = 1)
  
  return(NULL)
}

test_proximalTree <- function() {
  m = 10;n = 1000
  U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)

  .printf("First tree example\n")
# Example 1 of tree structure
# tree structured groups:
# g1= {0 1 2 3 4 5 6 7 8 9}
# g2= {2 3 4}
# g3= {5 6 7 8 9}
  own_variables = as.vector(c(0,2,5),mode= 'integer')
  N_own_variables = as.vector(c(2,3,5),mode= 'integer')
  eta_g = as.vector(c(1,1,1),mode = 'double')
  groups = matrix(as.vector(c(0,0,0,
    1,0,0,
    1,0,0),mode='logical'),ncol = 3,byrow = T)
  groups = as(groups,'CsparseMatrix')
  tree = list('eta_g'= eta_g,'groups' = groups,'own_variables' = own_variables,
    'N_own_variables' = N_own_variables)
  .printf('\ntest prox tree-linf\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-l2', pos = FALSE,intercept = FALSE)),n = 1)

  .printf('\ntest prox tree-linf\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-linf', pos = FALSE,intercept = FALSE)),n = 1)

  .printf('Second tree example\n')
# Example 2 of tree structure
# tree structured groups:
# g1= {0 1 2 3 4 5 6 7 8 9}    root(g1) = { };
# g2= {0 1 2 3 4 5}            root(g2) = {0 1 2};
# g3= {3 4}                    root(g3) = {3 4};
# g4= {5}                      root(g4) = {5};
# g5= {6 7 8 9}                root(g5) = { };
# g6= {6 7}                    root(g6) = {6 7};
# g7= {8 9}                    root(g7) = {8};
# g8 = {9}                     root(g8) = {9};
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
  .printf('\ntest prox tree-l0\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-l0', pos = FALSE,intercept = FALSE)),n = 1)
 
  .printf('\ntest prox tree-l2\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-l2', pos = FALSE,intercept = FALSE)),n = 1)
 
  .printf('\ntest prox tree-linf\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-linf', pos = FALSE,intercept = FALSE)),n = 1)
 
# mexProximalTree also works with non-tree-structured regularization functions
  .printf('\nprox l1, intercept, positivity constraint\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'l1', pos = TRUE,intercept = TRUE)),n = 1)

  .printf('\nprox multi-task tree\n')
  alpha = Xtest1('spams',quote(spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,lambda2 = 0.1,regul = 'multi-task-tree', pos = FALSE,intercept = FALSE)),n = 1)
  
  return(NULL)
}

test_prox.tests = list ('fistaFlat' = test_fistaFlat,
    'fistaGraph' = test_fistaGraph,
    'fistaTree' = test_fistaTree,
    'proximalFlat' = test_proximalFlat,
    'proximalGraph' = test_proximalGraph,
    'proximalTree' = test_proximalTree
  )
