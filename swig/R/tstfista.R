library("spams",lib.loc = "./lib")
.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}

  param = list('numThreads' = 1,'verbose' = FALSE,
    'lambda' = 0.05, 'it0' = 10, 'max_it' = 200,
    'L0' = 0.1, 'tol' = 1e-3, 'intercept' = FALSE,
    'pos' = FALSE)
  set.seed(0)
  m = 100;n = 200
  X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
#  l = scan(file= "../python/datax",what =  double(0))
#  X = matrix(l,nrow = m,ncol = n,byrow = TRUE)
  
  X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
  X = spams.Normalize(X)
  Y = matrix(rnorm(m),nrow = m,ncol = 1,byrow = FALSE)
#  l = scan(file= "../python/datay",what =  double(0))
#  Y = matrix(l,nrow = m,ncol = 1,byrow = TRUE)
  Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
  Y = spams.Normalize(Y)
  W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))
  param['compute_gram'] = TRUE
  param['verbose'] = TRUE
  param['loss'] = 'square'
#  param['regul'] = 'l1'
  param['ista'] = FALSE
  param['subgrad'] = TRUE
  param['a'] = 0.1
  param['b'] = 1000 # arbitrary parameters
  max_it = param['max_it']
  it0 = param['it0']
  param['max_it'] = 500
  param['it0'] = 50
  res = spams.FistaFlat(Y,X,W0,param,TRUE)
  W = res[[1]]
  optim_info = res[[2]]
  .printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4])
  l = scan(file= "../python/dataw",what =  double(0))
  W1 = matrix(l,nrow = n,ncol = 1,byrow = TRUE)
  err = max(abs(W1 -W))
  cat(sprintf("ERR %f\n",err))
####
  Y = ceiling(5 * matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000,byrow = FALSE)) - 1
  param['loss'] = 'multi-logistic'
  .printf("\nFISTA + Multi-Class Logistic l1\n")
  nclasses = max(Y) + 1
  W0 = matrix(0,nrow = ncol(X),nclasses * ncol(Y),byrow = FALSE)
  res = spams.FistaFlat(Y,X,W0,param,TRUE)
  W = res[[1]]
  optim_info = res[[2]]
  cat(sprintf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info[1],optim_info[3],optim_info[4]))
# can be used of course with other regularization functions, intercept,...
