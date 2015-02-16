library("spams",lib.loc = "./lib")
.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}

  m = 10;n = 1000
#  U = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
l = scan(file= "../../xxD",what =  double(0))
U = matrix(l,nrow = m,ncol = n,byrow = TRUE)

# Example 1 of tree structure
# tree structured groups:
# g1= {0 1 2 3 4 5 6 7 8 9}
# g2= {2 3 4}
# g3= {5 6 7 8 9}
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
  alpha = spams.proximalTree(U,tree,FALSE,numThreads = -1,verbose = TRUE,lambda1 = 0.1,regul = 'tree-l0', pos = FALSE,intercept = FALSE)
 
