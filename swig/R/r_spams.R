#################
.verif_enum <- function(arg,ename,msg) {
  defName = paste(".__E___", ename, sep = "")
  l = eval.parent(parse(text = sprintf("l <- %s",defName)))
  if (! (arg %in% names(l)))
    stop("ERROR : bad enum value ",msg,"\n")
}
###########  linalg ##############

spams.sort <- function(X,mode=T) {
     Y = c(X)
     sort(Y,mode)
     return(Y)
}

spams.calcAAt <- function(A) {
  m = nrow(A)
# the matrix creation is about 10 times faster wtih c(0)
#  AAt = matrix(rep(0,m * m),nrow = m,ncol = m)
  AAt = matrix(c(0),nrow = m,ncol = m)
  AAt(A,AAt)
  return(AAt)
}
spams.calcXAt <- function(X,A) {
  m = nrow(X)
  n = nrow(A)
#  XAt = matrix(rep(0,m * n),nrow = m,ncol = n)
  XAt = matrix(c(0),nrow = m,ncol = n)
  XAt(A,X,XAt)
  return(XAt)
}

spams.mult <- function(X,Y,transX = FALSE, transY = FALSE) {
    if(! is.matrix(X) || ! is.matrix(Y)) {
        stop("Arguments  must be matrices")
    }
    if(transX)
        m = ncol(X)
    else
	m = nrow(X)
    if(transY)
        n = nrow(Y)
    else
	n = ncol(Y)
#    XY = matrix(rep(0,m * n),nrow = m,ncol = n)
    XY = matrix(c(0),nrow = m,ncol = n)
    mult(X,Y,XY,transX,transY,1,0)
    return(XY)
}

spams.calcXY <- function(X,Y) {
    return(spams.mult(X,Y,F,F))
}
spams.calcXYt <- function(X,Y) {
    return(spams.mult(X,Y,F,T))
}
spams.calcXtY <- function(X,Y) {
    return(spams.mult(X,Y,T,F))
}
spams.bayer <- function(X,offset) {
  y =c(X)
  applyBayerPattern(y,offset)
  return(y)
}


spams.conjGrad <- function(A,b,x0 = NULL,tol = 1e-10,itermax = NULL) {
  n = ncol(A)
  if(is.null(x0)) {
    x = as.vector(matrix(c(0),ncol = n))
  } else {
    x = c(x0)
  }
  if(is.null(itermax)) {
    itermax = n
  }
  conjugateGradient(A,b,x,tol,itermax)
  return(x)
}

spams.invSym <- function(A) {
      B = matrix(A,nrow = nrow(A),ncol = ncol(A))
      invSym(B)
      return(B)
}

spams.normalize <- function(X) {
      B = matrix(X,nrow = nrow(X),ncol = ncol(X))
      normalize(B)
      return(B)
}

##################################################

###########  decomp ##################

spams.sparseProject <- function(U,thrs = 1.0,mode = 1,lambda1 = 0.0,lambda2 = 0.0,lambda3 = 0.0,pos = FALSE,numThreads = -1) {
  m = nrow(U)
  n = ncol(U)
##  paramlist = list('thrs' = 1.0,'mode' = 1,'lambda1' = 0.0,'lambda2' = 0.0,'lambda3' = 0.0,'pos' = FALSE,'numThreads' = -1)
##  params = .param_struct(paramlist,param)
#  V = matrix(rep(0,m * n),nrow = m,ncol = n)
  V = matrix(c(0),nrow = m,ncol = n)

##  .mycall('sparseProject',c('U','V',params))
  sparseProject(U,V,thrs,mode,lambda1,lambda2,lambda3,pos,numThreads)
  return(V)
}

spams.lasso <- function(X,D= NULL,Q = NULL,q = NULL,return_reg_path = FALSE,L= -1,lambda1= NULL,lambda2= 0.,
                        mode= 'PENALTY',pos= FALSE,ols= FALSE,numThreads= -1,
                 max_length_path= -1,verbose=FALSE,cholesky= FALSE) {
#  require('Matrix')
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h

  if (! is.null(Q)) {
    if (is.null(q)) {
      stop("ERROR lasso : q is needed when Q is given\n")
    }
  } else {
    if(is.null(D)) {
      stop("ERROR lasso : you must give D or Q and q\n")
    }
  }
  if(is.null(lambda1)) {
    stop("ERROR lasso : lambda1 must be defined\n")
  }
  .verif_enum(mode,'constraint_type','mode in Lasso')
  path = NULL
  x = NULL
  if(! is.null(q)) {
##    x = do.call(spams_wrap.lassoQq,c(list(X,D,q,0,return_reg_path),params))
    x = lassoQq(X,Q,q,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
##    x = .mycall('lassoQq',c('X','D','q',0,'return_reg_path',params))
  } else {
##    x = do.call(spams_wrap.lassoD,c(list(X,D,0,return_reg_path),params))
    x = lassoD(X,D,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
##    x = .mycall('lassoD',c('X','D',0,'return_reg_path',params))
  }
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)

}

spams.lassoMask <- function(X,D,B,L= -1,lambda1= NULL,lambda2= 0.,
                        mode= 'PENALTY',pos= FALSE,numThreads= -1,
                 verbose=FALSE) {
#  require('Matrix')
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h

  if(is.null(lambda1)) {
    stop("ERROR lassoMask : lambda1 must be defined\n")
  }
  .verif_enum(mode,'constraint_type','mode in Lasso')
  x = lassoMask(X,D,B,L,lambda1,lambda2,mode,pos,numThreads,verbose)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  cat("LASSO : ", length(shape),"\n")
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)
}

spams.lassoWeighted <- function(X,D,W,L= -1,lambda1= NULL,
                        mode= 'PENALTY',pos= FALSE,numThreads= -1,
                 verbose=FALSE) {
#  require('Matrix')
    # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
    # will be set in spams.h

  if(is.null(lambda1)) {
    stop("ERROR lassoWeighted : lambda1 must be defined\n")
  }
  .verif_enum(mode,'constraint_type','mode in Lasso')
  x = lassoWeighted(X,D,W,L,lambda1,mode,pos,numThreads,verbose)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)
}

spams.omp <- function(X,D,L = NULL,eps = NULL,lambda1 = NULL,return_reg_path = FALSE, numThreads = -1) {
  path = NULL
  given_L = FALSE
  given_eps = FALSE
  given_lambda1 = FALSE
  if (is.null(L)) {
    L = as.vector(c(0),mode='integer')
  } else {
    given_L = TRUE
    if(length(L) == 1 && ! is.integer(L)) {
      L = as.vector(c(L),mode='integer')
    }
  }
  if (is.null(eps)) {
    eps = as.vector(c(0.),mode='double')
  } else {
    given_eps = TRUE
  }
  if (is.null(lambda1)) {
    lambda1 = as.vector(c(0.),mode='double')
  } else {
    given_lambda1 = TRUE
  }
  
#  if(! is.vector(eps)) {
#    eps = as.vector(c(eps),mode='double')
#  }
  x = omp(X,D,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1, numThreads)
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)
}

spams.ompMask <- function(X,D,B,L = NULL,eps = NULL,lambda1 = NULL,return_reg_path = FALSE, numThreads = -1) {
  path = NULL
  given_L = FALSE
  given_eps = FALSE
  given_lambda1 = FALSE
  if (is.null(L)) {
    L = as.vector(c(0),mode='integer')
  } else {
    given_L = TRUE
    if(length(L) == 1 && ! is.integer(L)) {
      L = as.vector(c(L),mode='integer')
    }
  }
  if (is.null(eps)) {
    eps = as.vector(c(0.),mode='double')
  } else {
    given_eps = TRUE
  }
  if (is.null(lambda1)) {
    lambda1 = as.vector(c(0.),mode='double')
  } else {
    given_lambda1 = TRUE
  }

  x = ompMask(X,D,B,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1, numThreads)
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)
}

spams.cd  <- function(X,D,A0,lambda1 = NULL,mode= 'PENALTY',itermax=100,tol = 0.001,numThreads =-1) {
  if(is.null(lambda1)) {
    stop("ERROR cd : lambda1 must be defined\n")
  }
  x = cd(X,D,A0,lambda1,mode,itermax,tol,numThreads)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)

}

spams.somp  <- function(X,D,list_groups,L = NULL,eps = 0.,numThreads = -1){
  if(is.null(L)) {
    stop("somp : L must be defined\n")
  }
  x = somp(X,D,list_groups,L,eps,numThreads)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)
  
}

spams.l1L2BCD  <- function(X,D,alpha0,list_groups,lambda1 = NULL, mode= 'PENALTY', itermax = 100, tol = 1e-3,numThreads = -1) {
  if(is.null(lambda1)) {
    stop("l1L2BCD : lambda1 must be defined\n")
  }
  alpha = matrix(alpha0,nrow = nrow(alpha0),ncol = ncol(alpha0))
  l1L2BCD(X,D,alpha,list_groups,lambda1,mode,itermax,tol,numThreads)
  return(alpha)
  
}

  

###########  END decomp ##############
##################################################

###########  prox ##################
# W = FistaFlat(Y,X,W0,param,return_optim_info = False)
# (W,optim_info) = FistaFlat(Y,X,W0,param,return_optim_info = True)
spams.fistaFlat <- function(Y,X,W0,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),size_group=1,groups = NULL,sqrt_step=TRUE,transpose=FALSE,linesearch_mode=0) {

  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  W = matrix(c(0),nrow = m,ncol = n)
#  optim_info = do.call(solver,c(list(Y,X,W0,W),params))
##  optim_info = .mycall('fistaFlat',c('Y','X','W0','W',params))
  if (is.null(groups)) {
    groups = vector(mode = 'integer')
  }
  optim_info = fistaFlat(Y,X,W0,W,groups,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose,linesearch_mode)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
}

spams.fistaTree <- function(Y,X,W0,tree,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),size_group=1,sqrt_step=TRUE,transpose=FALSE,linesearch_mode=0) {
  if (length(tree) != 4) {
    stop("fistaTree : tree should be a list of 4 elements")
  }
  eta_g = tree[['eta_g']]
  groups = tree[['groups']]
  own_variables = tree[['own_variables']]
  N_own_variables = tree[['N_own_variables']]
  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  W = matrix(c(0),nrow = m,ncol = n)
  optim_info = fistaTree(Y,X,W0,W,eta_g,groups,own_variables,N_own_variables,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose,linesearch_mode)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
}

spams.fistaGraph <- function(Y,X,W0,graph,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),size_group=1,sqrt_step=TRUE,transpose=FALSE,linesearch_mode=0) {
  if (length(graph) != 3) {
    stop("fistaGraph : graph should be a list of 3 elements")
  }
  eta_g = graph[['eta_g']]
  groups = graph[['groups']]
  groups_var = graph[['groups_var']]
  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  W = matrix(c(0),nrow = m,ncol = n)
  optim_info = fistaGraph(Y,X,W0,W,eta_g,groups,groups_var,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose,linesearch_mode)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
  
}
  
spams.proximalFlat <- function(U,return_val_loss = FALSE,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=FALSE,resetflow=FALSE,regul="",verbose=FALSE,
                 pos=FALSE,clever=TRUE,size_group=1,groups = NULL,transpose=FALSE) {

  if (is.null(groups)) {
    groups = vector(mode = 'integer')
  }
  eval = return_val_loss
  m = nrow(U)
  n = ncol(U)
#  alpha = matrix(rep(0,m * n),nrow = m,ncol = n)
  alpha = matrix(c(0),nrow = m,ncol = n)
##  val_loss = .mycall('proximalFlat',c('U','alpha',params))
  val_loss = proximalFlat(U,alpha,groups,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
  if(return_val_loss == TRUE)
    return(list(alpha,val_loss))
  else
    return (alpha)
  
}

spams.proximalTree <- function(U,tree,return_val_loss = FALSE,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=FALSE,resetflow=FALSE,regul="",verbose=FALSE,
                 pos=FALSE,clever=TRUE,size_group=1,transpose=FALSE) {
  eval = return_val_loss
  if (length(tree) != 4) {
    stop("proximalTree : tree should be a list of 4 elements")
  }
  eta_g = tree[['eta_g']]
  groups = tree[['groups']]
  own_variables = tree[['own_variables']]
  N_own_variables = tree[['N_own_variables']]
  m = nrow(U)
  n = ncol(U)
#  alpha = matrix(rep(0,m * n),nrow = m,ncol = n)
  alpha = matrix(c(0),nrow = m,ncol = n)
##  val_loss = .mycall('proximalFlat',c('U','alpha',params))
  val_loss = proximalTree(U,alpha,eta_g,groups,own_variables,N_own_variables,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
  if(return_val_loss == TRUE)
    return(list(alpha,val_loss))
  else
    return (alpha)
  
}

spams.proximalGraph <- function(U,graph,return_val_loss = FALSE,numThreads =-1,lambda1=1.0,lambda2=0.,
                 lambda3=0.,intercept=FALSE,resetflow=FALSE,regul="",verbose=FALSE,
                 pos=FALSE,clever=TRUE,size_group=1,transpose=FALSE) {
  eval = return_val_loss
  if (length(graph) != 3) {
    stop("proximalGraph : graph should be a list of 4 elements")
  }
  eta_g = graph[['eta_g']]
  groups = graph[['groups']]
  groups_var = graph[['groups_var']]
  m = nrow(U)
  n = ncol(U)
  alpha = matrix(c(0),nrow = m,ncol = n)
  val_loss = proximalGraph(U,alpha,eta_g,groups,groups_var,numThreads ,lambda1,lambda2,lambda3,intercept,resetflow,regul,verbose,pos,clever,eval,size_group,transpose)
  if(return_val_loss == TRUE)
    return(list(alpha,val_loss))
  else
    return (alpha)
  
}

###########  END prox ##############
##################################################

###########  dictLearn ##################
.TrainDL <- function(X,return_model= NULL,model = NULL,in_memory= FALSE,D = NULL,
                     graph= NULL, tree = NULL,numThreads = -1,
                     tol = 0.000001,fixed_step = True,ista = False,
                     batchsize = -1,K= -1,lambda1= NULL,lambda2= 10e-10,lambda3 = 0.,iter=-1,t0=1e-5,
                     mode='PENALTY',regul= "none",posAlpha=FALSE,posD=FALSE,
                     expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,verbose=TRUE,
                     gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1,stochastic_deprecated=FALSE,
                     modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  # We can only have simple objects in the param list of .mycall

  if (is.null(D)) {
    D = matrix(c(0.),nrow = 0,ncol=0)
  }
  if(is.null(lambda1)) {
    stop("ERROR trainDL : lambda1 must be defined\n")
  }
  
  .verif_enum(modeD,'constraint_type_D','modeD in TrainDL')
  .verif_enum(mode,'constraint_type','mode in TrainDL')
  if(is.null(graph) && is.null(tree)) {
    eta_g = vector(mode = 'double',length = 0)
    groups = as(matrix(c(FALSE),nrow = 1,ncol = 1),'CsparseMatrix')
  }
  if( is.null(tree)) {
    own_variables = vector(mode= 'integer',length = 0)
    N_own_variables = vector(mode= 'integer',length = 0)
  } else {
    eta_g = tree[['eta_g']]
    groups = tree[['groups']]
    own_variables = tree[['own_variables']]
    N_own_variables = tree[['N_own_variables']]
    if (is.null(eta_g) || is.null(groups) || is.null(own_variables) ||
        is.null(N_own_variables)) {
      stop("structTrainDL : incorrect tree structure\n")
    }
    if(! is.null(graph)) {
      stop("structTrainDL : only one of tree or graph can be given\n")
    }
  }
  if( is.null(graph)) {
    groups_var = as(matrix(c(FALSE),nrow = 1,ncol = 1),'CsparseMatrix')
  } else {
    eta_g = graph[['eta_g']]
    groups = graph[['groups']]
    groups_var = graph[['groups_var']]
    if(is.null(eta_g) || is.null(groups) || is.null(groups_var)) {
      stop("structTrainDL : incorrect graph structure\n")
    }
  }
  if (is.null(model)) {
    m_A = matrix(c(0.),nrow = 0,ncol=0)
    m_B = matrix(c(0.),nrow = 0,ncol=0)
    m_iter = 0
  } else {
    m_A = model[['A']]
    m_B = model[['B']]
    m_iter = model[['iter']]
  }
#  x =  .mycall('alltrainDL',c('X',in_memory,0,0,0,'return_model','m_A','m_B','m_iter','D',params))
  x = alltrainDL(
    X,in_memory,return_model,m_A,m_B,m_iter,D,
    eta_g, groups, groups_var, own_variables, N_own_variables,
    numThreads,tol,fixed_step,ista,batchsize,K,lambda1,lambda2,lambda3,
    iter,t0,mode,regul,posAlpha,posD,
    expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,
    stochastic_deprecated,modeParam,batch,log_deprecated,logName)
  
  D = x[[1]]
  if(return_model) {
    A = x[[2]]
    B = x[[3]]
    iter = x[[4]]
    model = list( 'A' = A, 'B' = B, 'iter' = iter[1])
    return(list(D,model))
  } else {
    return(D)
  }
  
}

spams.trainDL <- function(X,return_model= FALSE,model= NULL,D = NULL,numThreads = -1,batchsize = -1,
            K= -1,lambda1= NULL,lambda2= 10e-10,iter=-1,t0=1e-5,mode='PENALTY',
                 posAlpha=FALSE,posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,verbose=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=NULL,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  
  if (is.null(iter_updateD)) {
    if(batch) {
      iter_updateD = 5
    } else {
      iter_updateD = 1
    }
  }
  lambda3 = 0.
  regul = "undef"

  return (.TrainDL(X,return_model,model,FALSE,D,NULL,NULL,numThreads,0.000001,TRUE,FALSE,batchsize,K,lambda1,lambda2,lambda3,iter,t0,mode,regul,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName))
}

spams.trainDL_Memory <- function(X,D = NULL,numThreads = -1,batchsize = -1,
            K= -1,lambda1= NULL,iter=-1,t0=1e-5,mode='PENALTY',
                 posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=1,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  lambda2 = 10e-10
  verbose = FALSE
  posAlpha = FALSE
  lambda3 = 0.
  regul = "undef"

  return (.TrainDL(X,FALSE,NULL,TRUE,D,NULL,NULL,numThreads,0.000001,TRUE,FALSE,batchsize,K,lambda1,lambda2,lambda3,iter,t0,mode,regul,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName))
}

spams.structTrainDL <- function(X,return_model= FALSE,model= NULL,D = NULL,
               graph = NULL, tree = NULL,numThreads = -1,
               tol = 0.000001,fixed_step = TRUE,ista = FALSE,batchsize = -1,
               K= -1,lambda1= NULL,lambda2= 10e-10,lambda3 = 0.,iter=-1,t0=1e-5,regul='none',
               posAlpha=FALSE,posD=FALSE,expand=FALSE,modeD='L2',whiten=FALSE,clean=TRUE,
               verbose=TRUE,gamma1=0.,gamma2=0.,rho=1.0,iter_updateD=NULL,stochastic_deprecated=FALSE,modeParam=0,batch=FALSE,log_deprecated=FALSE,logName='') {
  
  if (is.null(iter_updateD)) {
    if(batch) {
      iter_updateD = 5
    } else {
      iter_updateD = 1
    }
  }

  return (.TrainDL(X,return_model,model,FALSE,D,graph,tree,numThreads,0.000001,TRUE,FALSE,batchsize,K,lambda1,lambda2,lambda3,iter,t0,'FISTAMODE',regul,posAlpha,posD,expand,modeD,whiten,clean,verbose,gamma1,gamma2,rho,iter_updateD,stochastic_deprecated,modeParam,batch,log_deprecated,logName))
}

spams.archetypalAnalysis <- function(X,p = 10, Z0 = NULL, returnAB = FALSE, robust=FALSE, epsilon=1e-3, computeXtX=TRUE, stepsFISTA=3, stepsAS=50, randominit=TRUE,numThreads=-1)  {
   if (is.null(Z0)) {
      x = archetypalAnalysis(X=X, p=p, robust=robust, epsilon=epsilon, computeXtX = computeXtX, stepsFISTA = stepsFISTA, stepsAS = stepsAS, randominit = randominit, numThreads = numThreads)
   } else {
      x = archetypalAnalysisInit(X=X, Z0=Z0, robust=robust, epsilon=epsilon, computeXtX = computeXtX, stepsFISTA = stepsFISTA, stepsAS = stepsAS, numThreads = numThreads)
   }
   Z = x[[1]]
   if (returnAB) {
      indptr = x[[2]][[1]]
      indices = x[[2]][[2]]
      data = x[[2]][[3]]
      shape = x[[2]][[4]]
      A = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
      indptr2 = x[[3]][[1]]
      indices2 = x[[3]][[2]]
      data2 = x[[3]][[3]]
      shape2 = x[[3]][[4]]
      B = sparseMatrix(i = indices2, p = indptr2, x = data2,dims = shape2, index1 = FALSE)
      return (list(Z,A,B))
   } else {
      return(Z)
   }
}

spams.decompSimplex <- function(X, Z, computeXtX = FALSE, numThreads = -1) {
   x = decompSimplex(X,Z,computeXtX,numThreads)
   indptr = x[[1]]
   indices = x[[2]]
   data = x[[3]]
   shape = x[[4]]
   A = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
   return(A)
}

spams.nmf <- function(X,return_lasso= FALSE,model= NULL,numThreads = -1,batchsize = -1,K= -1,
    iter=-1,t0=1e-5,clean=TRUE,rho=1.0,modeParam=0,batch=FALSE) {
  lambda1 = 0

  U = spams.trainDL(X,model = model,numThreads = numThreads,batchsize = batchsize,
    K = K,iter = iter, t0 = t0, clean = clean, rho = rho,verbose=FALSE, 
    modeParam = modeParam,batch = batch, lambda1 = lambda1,
    mode = 'PENALTY', posAlpha=TRUE,posD=TRUE,whiten=FALSE)
  if (! return_lasso) {
    return(U)
  }
  if (! is.matrix(X)) {
    stop("sparse matrix for lasso not yet implemented")
  } else {
    V = spams.lasso(X,D = U,return_reg_path = FALSE, numThreads = numThreads,
                  lambda1 = lambda1,mode = 'PENALTY', pos=TRUE)
  }
  return (list(U,V))
}

spams.nnsc <- function(X,return_lasso= FALSE,model= NULL,lambda1= NULL,
                       numThreads = -1,batchsize = -1,K= -1,iter=-1,t0=1e-5,
                       clean=TRUE,rho=1.0,modeParam=0,batch=FALSE) {
  
  if(is.null(lambda1)) {
    stop("ERROR nnsc : lambda1 must be defined\n")
  }
  U = spams.trainDL(X,model = model,numThreads = numThreads,batchsize = batchsize,
                K = K,iter = iter, t0 = t0, clean = clean, rho = rho,verbose=FALSE, 
                modeParam = modeParam,batch = batch, lambda1 = lambda1,
                mode = 'PENALTY', posAlpha=TRUE,posD=TRUE,whiten=FALSE)
  if (! return_lasso) {
    return(U)
  }
  V = spams.lasso(X,D = U,return_reg_path = FALSE, numThreads = numThreads,
                  lambda1 = lambda1,mode = 'PENALTY')
  return (list(U,V))

}

###########  END dictLearn ##############
# utility functions
spams.im2col_sliding <- function(A,m,n,RGB = FALSE) {
  mm = nrow(A)
  nn = ncol(A)
  M = m * n
  N =  (mm - m + 1) * (nn -n + 1)
  B = matrix(c(0),nrow = M,ncol = N)
  im2col_sliding(A,B,m,n,RGB)
  return(B)
}

spams.displayPatches <- function(D) {
  V = 1
  n = nrow(D)
  K = ncol(D)
  sizeEdge = sqrt(n)
  if( as.integer(sizeEdge) != sizeEdge) {
    V = 3
    sizeEdge = sqrt(n/V)
  }
  p = 4.5
  M = max(D)
  m = min(D)
  if (m >= 0) {
    me = 0
    sig = sqrt(mean(D * D))
  } else {
    me = mean(D)
    sig = sqrt(mean((D - me) * (D -me)))
  }
  D = D - me
  D = matrix(mapply(function(x,y) min(x,y),mapply(function(x,y) max(x,y),D,-p * sig),p * sig),nrow = nrow(D),ncol = ncol(D))
  M = max(D)
  m = min(D)
  D = (D - m)/ (M - m)
  nb1 = sqrt(K)
  nBins = as.integer(nb1)
  if (nBins != nb1) {
    nBins = nBins + 1
  }
  tmp = array(c(0),dim = c((sizeEdge+1)*nBins+1,(sizeEdge+1)*nBins+1,V))
  mm = sizeEdge * sizeEdge
  for (ii in 1:nBins) {
    for (jj in 1:nBins) {
      io = ii
      jo = jj
      offsetx = 0
      offsety = 0
      ii = ((ii - 1 + offsetx) %% nBins) + 1
      jj = ((jj - 1 + offsety) %% nBins) + 1
      ix = (io-1)*nBins+jo
      if(ix > K) {
        break
      }
      patchCol = D[1:n,ix]
      patchCol = array(patchCol,c(sizeEdge,sizeEdge,V))
      x1 = (ii-1)*(sizeEdge+1) + 2
      x2 = ii*(sizeEdge+1)
      y1 = (jj-1)*(sizeEdge+1) + 2
      y2 = jj*(sizeEdge+1)
      tmp[x1:x2,y1:y2,] = patchCol
      ii = io
      jj = jo;
    }
  }
#
  return(tmp)
  
}

##################################################
spams.simpleGroupTree <- function(degrees) {
  return ( simpleGroupTree(degrees,length(degrees)))
}
spams.readGroupStruct <- function(file) {
  return ( readGroupStruct(file))
}
spams.groupStructOfString <- function(s) {
  return ( groupStructOfString(s))
}
spams.graphOfGroupStruct <- function(gstruct) {
  x = graphOfGroupStruct(gstruct)
  eta_g = x[[1]]
  lg = x[[2]]
  lv = x[[3]]
  indptr = lg[[1]]
  indices = lg[[2]]
  data = lg[[3]]
  shape = lg[[4]]
  groups = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  indptr = lv[[1]]
  indices = lv[[2]]
  data = lv[[3]]
  shape = lv[[4]]
  groups_var = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  graph = list('eta_g'= x[[1]],'groups' = groups,'groups_var' = groups_var)
  return(graph)
}
spams.treeOfGroupStruct <- function(gstruct) {
  x = treeOfGroupStruct(gstruct)
  nbvars = x[[1]]
  perm = x[[2]]
  eta_g = x[[3]]
  lg = x[[4]]
  own_variables= x[[5]]
  N_own_variables = x[[6]]
  indptr = lg[[1]]
  indices = lg[[2]]
  data = lg[[3]]
  shape = lg[[4]]
  groups = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  tree = list('eta_g'= eta_g,'groups' = groups,'own_variables' = own_variables,
                'N_own_variables' = N_own_variables)
  return (list(perm,tree,nbvars))
}
