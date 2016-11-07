od.D.SOCP <-
function(F, b, A, type, kappa, t.max){
# (MI)SOCP formulation of the D-optimality problem 
# with the linear constraints A*w<=b

# variables:
# w1, ..., wn /// 1...n
# J11, ..., Jmm /// n + 1...n + m^2
# t11, ...t1m, ..., tn1, ...tnm /// n + m^2 + 1...n + m^2 + n*m
# Z11, ..., Z1m, ..., Zn1, ..., Znm /// n + m^2 + n*m + 1...n + m^2 + 2*n*m
# t11 + w1, ...t1m + w1, ..., tn1 + wn, ...tnm + wn /// n + m^2 + 2*n*m + 1...n + m^2 + 3*n*m
# t11-w1, ...t1m-w1, ..., tn1-wn, ...tnm-wn /// n + m^2 + 3*n*m + 1...n + m^2 + 4*n*m
# newvars


geomeanG <- function(m, n, M){

cones.index <- NULL

if(m == 2){
  cones.index <- vector("list", 1)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  newvars <- 3 #number of new variables
  ncg <- 2 #number of equality constraints
  cn <- 1 #number of cones
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
}

###############################################
if(m == 3){
  cones.index <- vector("list", 3)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  ncg <- 6
  newvars <- 9
  cn <- 3
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, M + 9] <- 1/2;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, M + 9] <- -1/2;  Ag[4, M + 5] <- -1
  Ag[5, M + 3] <- 1/2;  Ag[5, M + 6] <- 1/2;  Ag[5, M + 7] <- -1
  Ag[6, M + 3] <- 1/2;  Ag[6, M + 6] <- -1/2;  Ag[6, M + 8] <- -1
}

###############################################
if(m == 4){
  cones.index <- vector("list", 3)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  ncg <- 6
  newvars <- 9
  cn <- 3
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, M + 3] <- 1/2;  Ag[5, M + 6] <- 1/2;  Ag[5, M + 7] <- -1
  Ag[6, M + 3] <- 1/2;  Ag[6, M + 6] <- -1/2;  Ag[6, M + 8] <- -1
}

################################################
if(m == 5){
  cones.index <- vector("list", 6)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cn <- 6
  newvars <- 18
  ncg <- 12
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, M + 18] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, M + 18] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, M + 3] <- 1/2;  Ag[7, M + 6] <- 1/2;  Ag[7, M + 10] <- -1
  Ag[8, M + 3] <- 1/2;  Ag[8, M + 6] <- -1/2;  Ag[8, M + 11] <- -1
  Ag[9, M + 9] <- 1/2;  Ag[9, M + 18] <- 1/2;  Ag[9, M + 13] <- -1
  Ag[10, M + 9] <- 1/2;  Ag[10, M + 18] <- -1/2;  Ag[10, M + 14] <- -1
  Ag[11, M + 12] <- 1/2;  Ag[11, M + 15] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 12] <- 1/2;  Ag[12, M + 15] <- -1/2;  Ag[12, M + 17] <- -1
}

###############################################
if(m == 6){
  cones.index <- vector("list", 6)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cn <- 6
  newvars <- 18
  ncg <- 12
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, M + 3] <- 1/2;  Ag[7, M + 6] <- 1/2;  Ag[7, M + 10] <- -1
  Ag[8, M + 3] <- 1/2;  Ag[8, M + 6] <- -1/2;  Ag[8, M + 11] <- -1
  Ag[9, M + 9] <- 1/2;  Ag[9, M + 18] <- 1/2;  Ag[9, M + 13] <- -1
  Ag[10, M + 9] <- 1/2;  Ag[10, M + 18] <- -1/2;  Ag[10, M + 14] <- -1
  Ag[11, M + 12] <- 1/2;  Ag[11, M + 15] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 12] <- 1/2;  Ag[12, M + 15] <- -1/2;  Ag[12, M + 17] <- -1
}

###############################################
if(m == 7){
  cones.index <- vector("list", 7)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cones.index[[7]] <- c(M + 19, M + 20, M + 21)
  cn <- 7
  newvars <- 21
  ncg <- 14
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, n + 6 * m + 7] <- 1;  Ag[7, M + 21] <- 1;  Ag[7, M + 10] <- -1
  Ag[8, n + 6 * m + 7] <- 1;  Ag[8, M + 21] <- -1;  Ag[8, M + 11] <- -1
  Ag[9, M + 3] <- 1/2;  Ag[9, M + 6] <- 1/2;  Ag[9, M + 13] <- -1
  Ag[10, M + 3] <- 1/2;  Ag[10, M + 6] <- -1/2;  Ag[10, M + 14] <- -1
  Ag[11, M + 9] <- 1/2;  Ag[11, M + 12] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 9] <- 1/2;  Ag[12, M + 12] <- -1/2;  Ag[12, M + 17] <- -1
  Ag[13, M + 15] <- 1/2;  Ag[13, M + 18] <- 1/2;  Ag[13, M + 19] <- -1
  Ag[14, M + 15] <- 1/2;  Ag[14, M + 18] <- -1/2;  Ag[14, M + 20] <- -1
}

###############################################
if(m == 8){
  cones.index <- vector("list", 7)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cones.index[[7]] <- c(M + 19, M + 20, M + 21)
  cn <- 7
  newvars <- 21
  ncg <- 14
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
  Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
  Ag[9, M + 3] <- 1/2;  Ag[9, M + 6] <- 1/2;  Ag[9, M + 13] <- -1
  Ag[10, M + 3] <- 1/2;  Ag[10, M + 6] <- -1/2;  Ag[10, M + 14] <- -1
  Ag[11, M + 9] <- 1/2;  Ag[11, M + 12] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 9] <- 1/2;  Ag[12, M + 12] <- -1/2;  Ag[12, M + 17] <- -1
  Ag[13, M + 15] <- 1/2;  Ag[13, M + 18] <- 1/2;  Ag[13, M + 19] <- -1
  Ag[14, M + 15] <- 1/2;  Ag[14, M + 18] <- -1/2;  Ag[14, M + 20] <- -1
}

###############################################
if(m == 9){
  cones.index <- vector("list", 11)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cones.index[[7]] <- c(M + 19, M + 20, M + 21)
  cones.index[[8]] <- c(M + 22, M + 23, M + 24)
  cones.index[[9]] <- c(M + 25, M + 26, M + 27)
  cones.index[[10]] <- c(M + 28, M + 29, M + 30)
  cones.index[[11]] <- c(M + 31, M + 32, M + 33)
  cn <- 11
  newvars <- 33
  ncg <- 22
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
  Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
  Ag[9, n + 8 * m + 9] <- 1;  Ag[9, M + 33] <- 1;  Ag[9, M + 13] <- -1
  Ag[10, n + 8 * m + 9] <- 1;  Ag[10, M + 33] <- -1;  Ag[10, M + 14] <- -1
  Ag[11, M + 3] <- 1/2;  Ag[11, M + 6] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 3] <- 1/2;  Ag[12, M + 6] <- -1/2;  Ag[12, M + 17] <- -1
  Ag[13, M + 9] <- 1/2;  Ag[13, M + 12] <- 1/2;  Ag[13, M + 19] <- -1
  Ag[14, M + 9] <- 1/2;  Ag[14, M + 12] <- -1/2;  Ag[14, M + 20] <- -1
  Ag[15, M + 15] <- 1/2;  Ag[15, M + 33] <- 1/2;  Ag[15, M + 22] <- -1
  Ag[16, M + 15] <- 1/2;  Ag[16, M + 33] <- -1/2;  Ag[16, M + 23] <- -1
  Ag[17, M + 18] <- 1/2;  Ag[17, M + 21] <- 1/2;  Ag[17, M + 25] <- -1
  Ag[18, M + 18] <- 1/2;  Ag[18, M + 21] <- -1/2;  Ag[18, M + 26] <- -1
  Ag[19, M + 24] <- 1/2;  Ag[19, M + 33] <- 1/2;  Ag[19, M + 28] <- -1
  Ag[20, M + 24] <- 1/2;  Ag[20, M + 33] <- -1/2;  Ag[20, M + 29] <- -1
  Ag[21, M + 27] <- 1/2;  Ag[21, M + 30] <- 1/2;  Ag[21, M + 31] <- -1
  Ag[22, M + 27] <- 1/2;  Ag[22, M + 30] <- -1/2;  Ag[22, M + 32] <- -1
}

###############################################
if(m == 10){
  cones.index <- vector("list", 11)
  cones.index[[1]] <- c(M + 1, M + 2, M + 3)
  cones.index[[2]] <- c(M + 4, M + 5, M + 6)
  cones.index[[3]] <- c(M + 7, M + 8, M + 9)
  cones.index[[4]] <- c(M + 10, M + 11, M + 12)
  cones.index[[5]] <- c(M + 13, M + 14, M + 15)
  cones.index[[6]] <- c(M + 16, M + 17, M + 18)
  cones.index[[7]] <- c(M + 19, M + 20, M + 21)
  cones.index[[8]] <- c(M + 22, M + 23, M + 24)
  cones.index[[9]] <- c(M + 25, M + 26, M + 27)
  cones.index[[10]] <- c(M + 28, M + 29, M + 30)
  cones.index[[11]] <- c(M + 31, M + 32, M + 33)
  cn <- 11
  newvars <- 33
  ncg <- 22
  Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
  Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
  Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
  Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
  Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
  Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
  Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
  Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
  Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
  Ag[9, n + 8 * m + 9] <- 1;  Ag[9, n + 9 * m + 10] <- 1;  Ag[9, M + 13] <- -1
  Ag[10, n + 8 * m + 9] <- 1;  Ag[10, n + 9 * m + 10] <- -1;  Ag[10, M + 14] <- -1
  Ag[11, M + 3] <- 1/2;  Ag[11, M + 6] <- 1/2;  Ag[11, M + 16] <- -1
  Ag[12, M + 3] <- 1/2;  Ag[12, M + 6] <- -1/2;  Ag[12, M + 17] <- -1
  Ag[13, M + 9] <- 1/2;  Ag[13, M + 12] <- 1/2;  Ag[13, M + 19] <- -1
  Ag[14, M + 9] <- 1/2;  Ag[14, M + 12] <- -1/2;  Ag[14, M + 20] <- -1
  Ag[15, M + 15] <- 1/2;  Ag[15, M + 33] <- 1/2;  Ag[15, M + 22] <- -1
  Ag[16, M + 15] <- 1/2;  Ag[16, M + 33] <- -1/2;  Ag[16, M + 23] <- -1
  Ag[17, M + 18] <- 1/2;  Ag[17, M + 21] <- 1/2;  Ag[17, M + 25] <- -1
  Ag[18, M + 18] <- 1/2;  Ag[18, M + 21] <- -1/2;  Ag[18, M + 26] <- -1
  Ag[19, M + 24] <- 1/2;  Ag[19, M + 33] <- 1/2;  Ag[19, M + 28] <- -1
  Ag[20, M + 24] <- 1/2;  Ag[20, M + 33] <- -1/2;  Ag[20, M + 29] <- -1
  Ag[21, M + 27] <- 1/2;  Ag[21, M + 30] <- 1/2;  Ag[21, M + 31] <- -1
  Ag[22, M + 27] <- 1/2;  Ag[22, M + 30] <- -1/2;  Ag[22, M + 32] <- -1
}

return(list(cones.index = cones.index,  newvars = newvars,  ncg = ncg,  Ag = Ag, cn = cn)) 
}






start <- as.numeric(proc.time()[3])
info <- paste("Running D.SOCP for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

requireNamespace('Matrix',quietly=TRUE)
  
if(requireNamespace('gurobi',quietly=TRUE)){
model  <-  list()

n <- dim(F)[1]
m <- dim(F)[2]

Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)

if(m > 10){
  stop("This procedure is only implemented for m<=10.")
}


model$modelsense <- "max"

da <- dim(A);
nc <- da[1] #number of constraints
M <- n + m^2 + 4*m*n #number of parameters

geo <- geomeanG(m, n, M)
varnum <- geo$newvars + M


###########linear constraints Aeq*x=beq################
# c1: sum(F[, i]*Z[i, ])=J
# c2: J: lower triangular
# c3: sum(t[i, j])<=J[j, j]
# c4: A*w<=b
# c5: tij + wi
# c6: tij-wi

# number of constraints
nC <- c(m * m, 1/2 * m * (m - 1), m, nc, m * n, m * n, geo$ncg)
sC <- sum(nC)
Aeq <- matrix(0, nrow=sum(nC), ncol=varnum)

# c1: sum(Fp[, i]*Z[i, ])=J
for (i in 1:(m*m))
  Aeq[i, i + n] <- -1 

for (i in 1:m){
  for (j in 1:n){
    for (ii in 1:m){
      Aeq[(i-1) * m + ii, (n + m^2 + m * n + (j-1) * m + ii)] <- 1/2 * Fp[j, i] 
    }
  }
}

# c2: J is lower triangular
ind <- 1
for (i in 1:(m-1)){
  for (j in i:(m-1)){
    Aeq[m*m + ind, n + (i-1)*m + j + 1] <- 1
    ind <- ind + 1
  }
}

# c3: sum(t[i, j])<=J[j, j]
for (i in 1:m){
  Aeq[m*m + m*(m-1)/2 + i, n + (i-1)*m + i] <- -1
  for (j in 1:n){
    Aeq[m*m + m*(m-1)/2 + i, n + m*m + (j-1)*m + i] <- 1
  }
}

# c4: A*w<=b
Aeq[(m ^ 2 + 1/2 * m * (m - 1) + m + 1) : (m ^ 2 + 0.5 * m * (m - 1) + m + nc), 1:n] <- A


#c5:tij + wi
S <- sum(nC[1:4])
for (i in 1:(n*m)){
  Aeq[S + i, n + m^2 + i] <- 1
  Aeq[S + i, n + m^2 + 2*n*m + i] <- -1
  Aeq[S + i, floor((i-1)/m) + 1] <- 1
}

#c6:tij-wi
S <- sum(nC[1:5])
for (i in 1:(n * m)){
  Aeq[S + i, n + m ^ 2 + i] <- 1
  Aeq[S + i, n + m ^ 2 + 3 * n * m + i] <- -1
  Aeq[S + i, floor((i-1)/m) + 1] <- -1
}

#c7:geomean
S <- sum(nC[1:6])
if(nC[7] > 0) Aeq[(S + 1):sC, ] <- geo$Ag

#build the model matrix
model$A <- Matrix::Matrix(Aeq, sparse=TRUE)

#right-hand side
rhs <- rep(0, sC)
rhs[(sum(nC[1:3]) + 1):(sum(nC[1:4]))] <- b
model$rhs <- rhs 

#inequalities sense
model$sense <- rep("=", sC)
for(i in (m * m + 1/2 * m * (m-1) + 1) : (m^2 + 1/2 * m * (m-1) + m + nc)) 
  model$sense[i] <- "<="


###########cones#######################
# Z[i, j]^2<=2*u[i, j]*w[i]
# geo$cones.index
model$cones <- vector("list", n + geo$cn)
model$cones <- geo$cones.index
for (i in 1 : (m * n))
  model$cones[[i + geo$cn]] <- c(n + m^2 + 2 * m * n + i, n + m ^ 2 + 3 * m * n + i, n + m ^ 2 + m * n + i)

#########lower bound###########
blx <- rep(-Inf, varnum) 
blx[1:n] <- 0 #w
blx[(n + m^2 + 1):(n + m^2 + n*m)] <- 0 #t
blx[(n + m^2 + 2*m*n + 1):(n + m^2 + 3*m*n)] <- 0 #t + w
# 1st indices of all cones in geomean
if (m == 2) 
  blx[M + 1] <- 0
if (m == 3) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
}
if (m == 4) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
}
if (m == 5) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
}
if (m == 6) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
}
if (m == 7) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0; blx[M + 19] <- 0
}
if (m == 8) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0; 
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0; blx[M + 19] <- 0
}
if (m == 9) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0; blx[M + 19] <- 0 
  blx[M + 22] <- 0; blx[M + 25] <- 0; blx[M + 28] <- 0; blx[M + 31] <- 0
}
if (m==10) {
  blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0 
  blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0; blx[M + 19] <- 0
  blx[M + 22] <- 0; blx[M + 25] <- 0; blx[M + 28] <- 0; blx[M + 31] <- 0
}

#upper bound
bux <- rep(Inf, varnum)
model$lb <- blx
model$ub <- bux


##############objective############
#maximize prod(J[j, j])
c <- rep(0, varnum)
c[varnum] <- 1
model$obj <- c

##########variable types############

vtypes <- rep("C", varnum)
if (type=="exact") {
  for (i in 1:n) vtypes[i] <- "I"
  model$lb[1:n] <- -0.01
}
model$vtypes <- vtypes

params  <-  list(TimeLimit=t.max,  OptimalityTol=1e-9,  FeasibilityTol=1e-9)

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
        Phi.best <- od.crit(F, w.best, crit="D")  
     }
  }
}

t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("D.SOCP finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)


list(w.best=w.best, Phi.best=Phi.best, status=status, t.act=t.act)

}else{
 return("Package gurobi is required to run this algorithm")
}
}

