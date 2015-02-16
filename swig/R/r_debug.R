spams.tst <- function(A,flag) {
  x = tst(0,flag,A)
  return (x)
}
spams.xtst <- function(A,flag) {
  require('Matrix')
  x = xtst(A,flag)
  # !!!
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  y = sparseMatrix(i = indices, p = indptr, x = data, index1 = FALSE)
  if(flag)
    return (list(y,x[[2]]))
  else
    return (y)
}
spams.ztst <- function(X,ret = FALSE) {
#  x = ztst(A,alpha)
#  x = do.call(ztst,c(list(A,alpha),params))
  x = .mycall('ztst',c('X',0,0,0,'ret'))
  if (ret) {
    iter<- x[[4]]
    x[[4]] <- iter[[1]]
    return(x)
  } else {
    return (x[[1]])
  }
}
