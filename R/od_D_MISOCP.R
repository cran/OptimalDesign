od_D_MISOCP <- function(Fx, b1, A1, b2, A2, b3, A3, bin  = FALSE,
                        type = c("exact", "approximate"), gap  = NULL,
                        t.max = Inf)
{
  if (!requireNamespace("gurobi", quietly = TRUE))
    stop("(MI)SOCP requires package gurobi.")

  type <- match.arg(type)

  n <- nrow(Fx)   # number of candidate points
  m <- ncol(Fx)   # number of parameters

  # ---- linear constraints on w ----
  k1 <- length(b1)
  k2 <- length(b2)
  k3 <- length(b3)
  b  <- c(b1, b2, b3)
  A  <- rbind(A1, A2, A3)
  sense <- c(rep("<=", k1), rep(">=", k2), rep("=", k3))
  nc <- nrow(A)

  
  # variables:
  #  w (n)
  #  lower-tri L (m^2)
  #  the 4*m*n block for Z's
  M  <- n + m^2 + 4 * m * n
  
  
  # generic SOC representation of geometric mean of m nonnegative numbers
  .geomeanG_any <- function(m, n, M) {
  if (m < 2L)
    stop("m must be >= 2 for the D-criterion formulation")

  # ---------- build the list of "conceptual" cones ----------
  cones <- list()
  L <- integer(0)     
  k <- 1L
  i <- 1L
  while (i + 1L <= m) {
    # cone k: pair of two REAL diagonals
    cones[[k]] <- list(type = "stage1", x = i, y = i + 1L)
    L <- c(L, k)
    k <- k + 1L
    i <- i + 2L
  }
  if (i <= m) {
    # odd leftover: (dm, t)
    cones[[k]] <- list(type = "stage1_last_with_t", x = i)
    L <- c(L, k)
    k <- k + 1L
  }

  # reduce the list L by pairwise combining,
  # if length(L) is odd, the last element is combined with t
  while (length(L) > 1L) {
    newL <- integer(0)
    j <- 1L
    while (j + 1L <= length(L)) {
      cones[[k]] <- list(type = "later", xcone = L[j], ycone = L[j + 1L])
      newL <- c(newL, k)
      k <- k + 1L
      j <- j + 2L
    }
    if (j <= length(L)) {
      cones[[k]] <- list(type = "later_with_t", xcone = L[j])
      newL <- c(newL, k)
      k <- k + 1L
    }
    L <- newL
  }

  cn <- length(cones)          # number of cones
  newvars <- 3L * cn           # each cone adds 3 vars
  ncg <- 2L * cn               # and 2 linear rows
  totvars <- M + newvars
  Ag <- matrix(0, nrow = ncg, ncol = totvars)
  t.idx <- M + newvars         # *last* new var is the global t
  row <- 1L

  # helper: output var of cone "k" is always its 3rd new var
  cone_out <- function(k) M + (k - 1L) * 3L + 3L

  for (idx in seq_len(cn)) {
    cone <- cones[[idx]]
    A.idx <- M + (idx - 1L) * 3L + 1L
    B.idx <- M + (idx - 1L) * 3L + 2L

    if (cone$type == "stage1") {
      # inputs are REAL diagonals
      x.idx <- n + (cone$x - 1L) * m + cone$x
      y.idx <- n + (cone$y - 1L) * m + cone$y
      # z^2 <= 4 x y  <=>  ||(x - y, 2z)|| <= x + y
      Ag[row,   x.idx] <- 1
      Ag[row,   y.idx] <- 1
      Ag[row,   A.idx] <- -1
      Ag[row+1, x.idx] <- 1
      Ag[row+1, y.idx] <- -1
      Ag[row+1, B.idx] <- -1

    } else if (cone$type == "stage1_last_with_t") {
      # inputs are REAL diagonal and final t
      x.idx <- n + (cone$x - 1L) * m + cone$x
      if (m == 3L) {
        # tiny special case: the original code did this for m = 3
        Ag[row,   x.idx] <- 1
        Ag[row,   t.idx] <- 1/2
        Ag[row,   A.idx] <- -1
        Ag[row+1, x.idx] <- 1
        Ag[row+1, t.idx] <- -1/2
        Ag[row+1, B.idx] <- -1
      } else {
        Ag[row,   x.idx] <- 1
        Ag[row,   t.idx] <- 1
        Ag[row,   A.idx] <- -1
        Ag[row+1, x.idx] <- 1
        Ag[row+1, t.idx] <- -1
        Ag[row+1, B.idx] <- -1
      }

    } else if (cone$type == "later") {
      # later levels: always use 1/2 factors
      x.idx <- cone_out(cone$xcone)
      y.idx <- cone_out(cone$ycone)
      Ag[row,   x.idx] <- 1/2
      Ag[row,   y.idx] <- 1/2
      Ag[row,   A.idx] <- -1
      Ag[row+1, x.idx] <- 1/2
      Ag[row+1, y.idx] <- -1/2
      Ag[row+1, B.idx] <- -1

    } else if (cone$type == "later_with_t") {
      x.idx <- cone_out(cone$xcone)
      Ag[row,   x.idx] <- 1/2
      Ag[row,   t.idx] <- 1/2
      Ag[row,   A.idx] <- -1
      Ag[row+1, x.idx] <- 1/2
      Ag[row+1, t.idx] <- -1/2
      Ag[row+1, B.idx] <- -1
    }

    row <- row + 2L
  }

  cones.index <- lapply(seq_len(cn), function(k)
    c(M + (k - 1L) * 3L + 1L,
      M + (k - 1L) * 3L + 2L,
      M + (k - 1L) * 3L + 3L))

  list(cones.index = cones.index,
       newvars     = newvars,
       ncg         = ncg,
       Ag          = Ag,
       cn          = cn)
}



  # generic geometric-mean cones 
  geo <- .geomeanG_any(m, n, M)

  varnum <- M + geo$newvars

  # counts of equality blocks
  nC <- c(m * m, 0.5 * m * (m - 1),  m,  nc,  m * n,  m * n, geo$ncg)
  sC <- sum(nC)

  Aeq <- matrix(0, nrow = sC, ncol = varnum)

  ## tie L_ij to Fx and Z (first m^2 rows)
  # diag and lower-tri structure
  for (i in 1:(m * m)) {
    Aeq[i, i + n] <- -1
  }

  # constraints: L = sum_j Fx[j, ] * Z_j  
  for (i in 1:m) {
    for (j in 1:n) {
      for (ii in 1:m) {
        Aeq[(i - 1) * m + ii,
             (n + m^2 + m * n + (j - 1) * m + ii)] <- 0.5 * Fx[j, i]
      }
    }
  }

  ## enforce lower-triangular structure (next m(m-1)/2 rows)
  ind <- 1
  for (i in 1:(m - 1)) {
    for (j in i:(m - 1)) {
      Aeq[m * m + ind, n + (i - 1) * m + j + 1] <- 1
      ind <- ind + 1
    }
  }

  ## link diag(L) with "sum of Z" (next m rows)
  base <- m * m + m * (m - 1) / 2
  for (i in 1:m) {
    Aeq[base + i, n + (i - 1) * m + i] <- -1
    for (j in 1:n) {
      Aeq[base + i,
           n + m * m + (j - 1) * m + i] <- 1
    }
  }

  ## put user constraints on w (next nc rows)
  start <- m^2 + 0.5 * m * (m - 1) + m + 1
  Aeq[start:(start + nc - 1), 1:n] <- A

  ## two big groups tying aux vars 
  S <- sum(nC[1:4])
  for (i in 1:(n * m)) {
    Aeq[S + i, n + m^2 + i] <- 1
    Aeq[S + i, n + m^2 + 2 * n * m + i] <- -1
    Aeq[S + i, floor((i - 1)/m) + 1] <- 1
  }
  S <- sum(nC[1:5])
  for (i in 1:(n * m)) {
    Aeq[S + i, n + m^2 + i] <- 1
    Aeq[S + i, n + m^2 + 3 * n * m + i] <- -1
    Aeq[S + i, floor((i - 1)/m) + 1] <- -1
  }

  ## add the geometric-mean rows
  S <- sum(nC[1:6])
  if (nC[7] > 0) {
    Aeq[(S + 1):sC, ] <- geo$Ag
  }

  # ---- model object ----
  model <- list()
  model$modelsense <- "max"
  model$A   <- Matrix::Matrix(Aeq, sparse = TRUE)

  rhs <- rep(0, sC)
  rhs[(sum(nC[1:3]) + 1):(sum(nC[1:4]))] <- b
  model$rhs   <- rhs
  model$sense <- rep("=", sC)

  # change some equalities to <= or user-defined sense
  for (i in (m^2 + 0.5 * m * (m - 1) + 1):(m^2 + 0.5 * m * (m - 1) + m)) {
    model$sense[i] <- "<="
  }
  for (i in (m^2 + 0.5 * m * (m - 1) + m + 1):
           (m^2 + 0.5 * m * (m - 1) + m + nc)) {
    model$sense[i] <- sense[i - (m^2 + 0.5 * m * (m - 1) + m)]
  }

  # ---- quadratic cones ----
  model$quadcon <- vector("list", n + geo$cn)

  # 1) geometric-mean cones
  for (i in seq_len(geo$cn)) {
    ind <- geo$cones.index[[i]]
    qc <- list()
    qc$Qc  <- Matrix::spMatrix(varnum, varnum, ind, ind, c(-1, 1, 1))
    qc$rhs <- 0
    model$quadcon[[i]] <- qc
  }

  # 2) the n*m little cones
  for (i in 1:(m * n)) {
    qc <- list()
    idxs <- c(n + m^2 + 2 * m * n + i,
              n + m^2 + 3 * m * n + i,
              n + m^2 + m * n + i)
    qc$Qc  <- Matrix::spMatrix(varnum, varnum,
                               idxs, idxs, c(-1, 1, 1))
    qc$rhs <- 0
    model$quadcon[[i + geo$cn]] <- qc
  }

  # ---- bounds ----
  blx <- rep(-Inf, varnum)
  blx[1:n] <- 0
  blx[(n + m^2 + 1):(n + m^2 + n * m)] <- 0
  blx[(n + m^2 + 2 * m * n + 1):(n + m^2 + 3 * m * n)] <- 0

  # IMPORTANT: first var of each geometric-mean cone must be >= 0
  first_cone_vars <- vapply(geo$cones.index, function(v) v[1], numeric(1))
  blx[first_cone_vars] <- 0

  bux <- rep(Inf, varnum)
  model$lb <- blx
  model$ub <- bux

  # ---- objective: maximize the final t ----
  cc <- rep(0, varnum)
  cc[varnum] <- 1
  model$obj <- cc

  # ---- variable types ----
  vtypes <- rep("C", varnum)
  if (type == "exact") {
    vtypes[1:n] <- if (bin) "B" else "I"
  } else if (type == "approximate" && bin) {
    model$ub[1:n] <- 1
  }
  model$vtype <- vtypes

  # ---- gurobi params ----
  if (!is.null(gap)) {
    params <- list(TimeLimit = t.max, MIPFocus = 1, MIPGap = gap)
  } else {
    params <- list(TimeLimit = t.max, MIPFocus = 1)
  }

  result <- gurobi::gurobi(model, params)
  gc()

  res <- result$x
  if (!is.vector(res) || !is.numeric(res) || !all(is.finite(res)) || length(res) < n) {
    w.best <- NULL
  } else {
    w.best <- res[1:n]
    if (type == "exact")
      w.best <- round(w.best)
  }

  list(w.best = w.best, status = result$status)
}