F.simplex <-
function(formula,n.factors,n.levels)
{
# The rules for creating the formula are the same as in the lm function but:
# 1) the formula must not contain the dependent variable;
# 2) the factors (independent variables) must be labeled x1,x2,...
# n.factors is an integer number >=2.
# n.levels is an integer number >=2.

# verify input
if(!is.numeric(n.factors) || !is.finite(n.factors) || (!length(n.factors)==1) || (n.factors <= 1)) stop("n.factors must be an integer greater than 1.")
if(!is.numeric(n.levels) || !is.finite(n.levels) || (!length(n.levels)==1) || (n.levels <= 1)) stop("n.levels must be an integer greater than 1.")



d <- n.factors
lst <- c()
labs <- c()
cmb <- t(combn(n.levels + n.factors - 2, n.factors - 1))
data <- as.data.frame((cbind(cmb, n.levels + n.factors - 1) - cbind(0, cmb) - 1) / (n.levels - 1))
for (i in 1:d){
    labs <- c(labs, paste("x", as.character(i), sep=""))
}
names(data) <- labs[1:n.factors]
mf <- model.frame(formula, data)
F <- as.data.frame(model.matrix(attr(mf, "terms"), mf))
return(as.matrix(F))
}
