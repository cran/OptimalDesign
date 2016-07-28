F.cube <-
function (formula, lower, upper, n.levels)
{
# The rules for creating the formula are the same as in the lm function but:
# 1) the formula must not contain the dependent variable;
# 2) the d factors (independent variables) must be labeled x1,x2,...
# lower is a d-dimensional vector of the smallest values of factors.
# upper is a d-dimensional vector of the largest values of factors.
# n.levels is a d-dimensional vector of the numbers of levels of each factor.

# verify input
if (!is.numeric(lower) || !is.finite(lower) || !is.vector(lower)) stop("Lower must be a real number or a real vector.")
if (!is.numeric(upper) || !is.finite(upper) || !is.vector(upper)) stop("Upper must be a real number or a real vector.")
if (!is.numeric(n.levels) || !is.finite(n.levels) || !is.vector(n.levels) || !all(n.levels==round(n.levels))) stop("n.levels must be integer vector.")
if (length(lower) != length(upper)) stop ("Lower and upper must be of the same length.")
if (any(lower >= upper)) stop ("lower[i] has to be smaller than upper[i] for all i.")
if (length(lower) != length(n.levels)) stop("n.levels must be of the same length as lower and upper.")



d <- length(n.levels)
lst <- c(); labs <- c()
for (i in 1:d){
    lst <- c(lst,list(seq(lower[i], upper[i], length=n.levels[i])))
    labs <- c(labs,paste("x", as.character(i), sep=""))
}
data <- expand.grid(lst)
names(data) <- labs
mf <- model.frame(formula, data)
F <- as.data.frame(model.matrix(attr(mf, "terms"), mf))
return(as.matrix(F))
}
