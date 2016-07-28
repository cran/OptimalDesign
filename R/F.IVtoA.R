F.IVtoA <- function (F, R = NULL) 
{
    n <- dim(F)[1]
    if (is.null(R)) 
      R <- 1:n
    k <- length(R)
    FR <- F[R, ]
    L <- od.infmat(FR, rep(1, k))
    print(paste("Condition number of L is", rcond(L)), quote = FALSE)
    print("Note that if L is singular, the problem of IV-optimality cannot be transformed to a problem of A-optimality; see the manual.", quote=FALSE)
    C.inv <- solve(chol(L))
    return(F %*% C.inv)
}
