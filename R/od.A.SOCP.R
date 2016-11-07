od.A.SOCP <-
function(F, b, A, type, kappa, t.max){
# misocp formulation of the A-optimality problem max 1/trace(inv(F' * diag(w) * F))
# with the linear constraints A * w<=b
# w0: starting design

# variables:
# w(1), ...w(n)
# mu(1), ...mu(n)
# w(1) + mu(1), ..., w(n) + mu(n)
# w(1)-mu(1), ..., w(n)-mu(n)
# 2Y(1, 1), ..., 2Y(n, 1), ..., 2Y(1, m), ..., 2Y(n, m)

start <- as.numeric(proc.time()[3])
info <- paste("Running A.SOCP for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

requireNamespace('Matrix',quietly=TRUE)
  
if(requireNamespace('gurobi',quietly=TRUE)){

model <- list()

n <- dim(F)[1]
m <- dim(F)[2]

Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)



da <- dim(A)
nc <- da[1] #number of constraints
M <- 4 * n + m * n #number of parameters


#######################constraints
# Aeq * x=beq
# sum(F(:, i) * Y(:, i)')=diag(m)
# A * w<=b

Aeq <- matrix(0, nrow = m * m + nc + 2 * n, ncol = M) 
A2 <- matrix(0, nrow = nc, ncol = M)
A2[, 1:n] <- A 
k <- 1

for (i in 1 : (m * m)){
  for (j in 1 : n){
    u <- i %% m
    if (u == 0) 
  u <- m
    Aeq[i, (3 + k) * n + j] = 0.5 * Fp[j, u]
  }
  if (i %% m == 0) k <- k + 1
}

A3 <- matrix(0, nrow = 2 * n, ncol = M)
for (i in 1:n){
  A3[i, i] <- 1
  A3[i, n + i] <- 1
  A3[i, 2 * n + i] <- -1
  A3[n + i, i] <- 1
  A3[n + i, n + i] <- -1
  A3[n + i, 3 * n + i] <- 1
}

beq <- rep(0, m * m + nc + 2 * n)
beq[1:(m^2)] <- as.vector(diag(m))

for (i in 1:nc){
  Aeq[m * m + i, ] <- A2[i, ]
  beq[m * m + i] <- b[i]
}

for (i in 1:(2 * n)) 
  Aeq[m * m + nc + i, ] <- A3[i, ]

model$A <- Matrix::Matrix(Aeq, sparse=TRUE)

model$rhs <- beq

#############
sense <- rep("=", m * m + nc + 2 * n)
for (i in (m * m + 1):(m * m + nc)) 
  sense[i] <- "<="
model$sense <- sense


############cones##############
# norm^2(Y(:, i))<=mu(i) * w(i)

model$cones <- vector("list", n)
for (i in 1:n){
 a <- 4:(3 + m)
 model$cones[[i]] <- c(2 * n + i, a * n + i, 3 * n + i)
}

#########lower bound###########
# w>=0
# mu>=0

lb <- rep(-Inf, M)
lb[1:(3 * n)] <- 0
model$lb <- lb

##############objective############
#minimize sum(mu)

c <- rep(0, M)
c[(n + 1):(2 * n)] <- 1
model$obj <- c

##########variable types############

vtypes <- rep("C", M)
if (type=="exact") {
  for (i in 1:n) 
    vtypes[i] <- "I"
}
model$vtypes <- vtypes


#model$start  <-  rep(NA,  dim(model$A)[2])
#model$start[1:n] <-  w0

model$modelsense <- "min"
params  <-  list(TimeLimit=t.max,  MIPGap=0)

result <- gurobi::gurobi(model, params)

status <- result$status
w.best <- NULL; Phi.best <- 0
w.temp <- result$x
if (is.numeric(w.temp) && is.vector(w.temp)) {
  if (length(w.temp) >= dim(A)[2]) {
     w.temp <- w.temp[1:dim(A)[2]]
     w.temp <- pmax(0, w.temp)
     if (type == "exact")      
        w.temp <- round(w.temp)
     if (sum(A %*% w.temp - b <= 1e-4) == dim(A)[1]) {
        w.best <- w.temp
        Phi.best <- od.crit(F, w.best, crit="A")  
     }
  }
}


t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("A.SOCP finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)


list(w.best=w.best, Phi.best=Phi.best, status=status, t.act=t.act)

}else{
 return("Package gurobi is required to run this algorithm")
}
}

