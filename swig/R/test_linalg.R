test_sort <- function() {
   x = rnorm(2000000,0,1)
   return(Xtest(quote(sort(x)),quote(spams.sort(x,TRUE))))
}

test_calcAAt <- function () {
  m = 200;n = 200000; d= 0.05
  A = rSpMatrix(m,n,floor(m * n * d))
  return(Xtest(quote(A %*% t(A)),quote(spams.calcAAt(A))))
}

test_calcXAt <- function () {
  m = 200;n = 200000; d= 0.05
  A = rSpMatrix(m,n,floor(m * n * d))
  X = matrix(rnorm(64 * n),nrow = 64,ncol = n)
  return(Xtest(quote(X %*% t(A)),quote(spams.calcXAt(X,A))))
}
matprod <- function(x,y) "%*%"(x,y)
test_calcXY <- function() {
    X = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 200,ncol = 20000,byrow = FALSE)
    return(Xtest(quote(X %*% Y),quote(spams.calcXY(X,Y))))
}

test_calcXYt <- function () {
    X = matrix(rnorm(64 * 200),nrow = 64,ncol = 200,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 20000,ncol = 200,byrow = FALSE)
    return(Xtest(quote(X %*% t(Y)),quote(spams.calcXYt(X,Y))))
  
}

test_calcXtY <- function () {
    X = matrix(rnorm(64 * 200),nrow = 200,ncol = 64,byrow = FALSE)
    Y = matrix(rnorm(200 * 20000),nrow = 200,ncol = 20000,byrow = FALSE)
    return(Xtest(quote(t(X) %*% Y),quote(spams.calcXtY(X,Y))))
}

test_bayer <- function () {
  X = rnorm(2000000,0,1)
  y2 = Xtest1('spams',quote(spams.bayer(X,0)),n=1)
  return(NULL)
}

test_conjGrad <- function () {
  A = matrix(rnorm(500 * 5000),nrow = 5000,ncol = 500)
  A = t(A) %*% A
  b = as.vector(rep(1,ncol(A)))
  x0 = b
  tol = 1e-4
  itermax = floor(0.5 * length(b))
  CGtest <- function(txt,expr) {
    tic = proc.time()
    for (i in 1:20)
      y = eval(expr)
    tac = proc.time()
    cat(sprintf("  Time (%s): %f\n",txt,(tac - tic)[['user.self']]))
    x1 = abs(b - (A %*% y))
    cat(sprintf("Mean error on b : %f\n\n",sum(x1) / length(b)))
    return(y)
  }
  y2 = CGtest("R",quote(solve(A,b)))
  y1 = CGtest("spams",quote(spams.conjGrad(A,b,x0,tol,itermax)))
  return(max(abs(y1 - y2)))
}

test_invSym <- function() {
    A = matrix(runif(1000 * 1000,0,1),nrow = 1000,ncol = 1000,byrow = FALSE)
    A = t(A) %*% A
    return(Xtest(quote(solve(A)),quote(spams.invSym(A))))
}
test_normalize <- function() {
    A = matrix(runif(100 * 1000,0,1),nrow = 100,ncol = 1000)
    y2 = Xtest1("spams",quote(spams.normalize(A)),n=1)
    return(NULL)
}

test_linalg.tests =list( 'sort' = test_sort,
  'calcAAt' = test_calcAAt,
  'calcXAt' = test_calcXAt,
  'calcXY' = test_calcXY,
  'calcXYt' = test_calcXYt,
  'calcXtY' = test_calcXtY,
  'bayer' = test_bayer,
  'conjGrad' = test_conjGrad,
  'invSym' = test_invSym,
  'normalize' = test_normalize
  )
