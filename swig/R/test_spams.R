#library("spams",lib.loc = "./lib")
# R_LIBS must be set
library("spams")

# for tests
simul <- FALSE

.printf <- function(...) {
  argv <- list(...)
  cat(sprintf(...))
}

usage <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  name <- substring(argv[grep("--file=", argv)], 8)
  sname <- basename(name)
  .printf("\nUsage : %s [test-or-group-name]+\n",name)
  .printf("  Run specified test or group of tests (all by default)\n")
  .printf("  Available groups and tests are:\n")
  for (m in modules) {
    .printf("%s :\n    ",m)
    s = sprintf("lst <- test_%s.tests",m)
    eval(parse(text= s))
    for (s in names(lst)) 
      .printf("%s ",s)
    .printf("\n")
  }
  .printf("\nExamples:\n")
  .printf("  Rscript %s linalg\n",name)
  .printf("  Rscript %s sort calcAAt\n",name)
  quit()
}

Xtest1 <- function(txt,expr,n = 2) {
  tic = proc.time()
  res = eval.parent(expr,n=n)
  tac = proc.time()
  t = tac-tic
  .printf("  Time (%s) : %.3fs\n",txt,t[['elapsed']])
  return(res)
}

Xtest <- function(expr1,expr2) {
  y1 = Xtest1('R',expr1)
  y2 = Xtest1('spams',expr2)
  return(max(abs(y2 - y1)))
}

allmodules <- list('dictLearn','linalg','decomp','prox')
modules = list()

for (s in allmodules) {
  LOADED <- TRUE
  source(paste("test_",s,".R",sep=''))
  if(! LOADED) {
    .printf("REMOVING %s\n",s)
  } else {
    modules = c(modules,s)
  }
}


rSpMatrix <- function(nrow, ncol, nnz,
                       rand.x = function(n) round(rnorm(nnz), 2))
  {
    ## Purpose: random sparse matrix
    ## --------------------------------------------------------------
    ## Arguments: (nrow,ncol): dimension
    ##          nnz  :  number of non-zero entries
    ##         rand.x:  random number generator for 'x' slot
    ## --------------------------------------------------------------
    ## Author: Martin Maechler, Date: 14.-16. May 2007
    # adapated by Jean-Paul Chieze to build a csc sparse matrix
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0,
              nnz <= nrow * ncol)
    A = spMatrix(nrow, ncol,
             i = sample(nrow, nnz, replace = TRUE),
             j = sample(ncol, nnz, replace = TRUE),
             x = rand.x(nnz))
    return(as(A,"dgCMatrix"))
  }

                    
run_test  <- function(testname,prog) {
  cat("** Running ",testname,"\n")
  if(simul) return()
  err = prog()
  if(! is.null(err))
    cat("  ERR = ",err,"\n")
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  lst = c(args)
  if(length(lst) == 0)
    lst = modules
  tic = proc.time()
  for (name in lst) {
    if (substr(name,1,1) == '-')
      usage()
    if(name %in% modules) {
      print(sprintf("**** %s ****",name))
      s = sprintf("lst <- test_%s.tests",name)
      eval(parse(text= s))
      for (m in names(lst)) {
        run_test(m,lst[[m]])
      }
    } else {
      found = FALSE
      for (m in modules) {
        s = sprintf("lst <- test_%s.tests",m)
        eval(parse(text= s))
        for (m in names(lst)) {
          if (m == name) {
            found = TRUE
            run_test(m,lst[[m]])
            break
          }
        }
        if (found) break
      }
      if(! found)
        print(sprintf("Test %s not found!",name))
    }

  }
  tac = proc.time()
  cat("\nTotal time:\n")
  print(tac - tic)
}

main()
