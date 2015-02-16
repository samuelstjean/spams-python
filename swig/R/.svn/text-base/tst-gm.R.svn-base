library("spams")

n=1
m=100
p=20
w0= matrix(rep(1,m),ncol=n)
x1= matrix(rnorm(m*p),m,p)
y1=matrix(sample(c(-1,1),m,replace=TRUE),ncol=n)

paramlasso= list("loss"="logistic","regul"="lasso","lambda"=0.1,"verbose" = TRUE)
Wlasso= spams.FistaFlat(y1,x1,w0,param=paramlasso,return_optim_info=FALSE)
