# scs.R — Sparse Cluster Selection: R API over the non-Gurobi C++ solver (src_scs.cpp).
# Idiom mirrors glm()/predict(): scs_fit() returns a fit object; predict()/stability() are layers on top.
#   fit  <- scs_fit(X, y, cluster, lambda = 5)                 # symmetric L0 budget
#   fit2 <- scs_fit(X, y, cluster, lambda = c(worst = 4, best = 2))  # asymmetric tail budgets
#   predict(fit, newdata = Xnew)                               # new-cluster prediction (CART soft map)
#   stability(fit, B = 100)                                    # per-cluster selection freq + FDR bound
suppressMessages({ library(Rcpp); library(RcppEigen) })
`%||%` <- function(a, b) if (is.null(a)) b else a
# Compile the C++ solver. Override SCS_CPP to point elsewhere; defaults to ./src_scs.cpp.
sourceCpp(getOption("scs.cpp", "src_scs.cpp"))

# D1 (resolved 2026-06-19): selection uses a small SCALE-RELATIVE ridge mu = eps*median(diag(Xfull'Xfull));
# reported coefficients come from an UNPENALIZED LS refit on the selected support (refit=TRUE) so the
# paper's "unshrunk" claim is literally true. Pass an explicit numeric `mu` to override the relative rule
# (back-compat: mu=1 reproduces the old ridge-everything behaviour). `obj` stays the ridged selection
# objective; `fixef`/`gamma` are the refit. See scs_eps_sensitivity() for the ε-sensitivity report.
scs_refit_unshrunk <- function(Xf, y, A, sel_i, lev) {
  pf <- ncol(Xf)
  if (length(sel_i) == 0) {
    cf <- stats::lm.fit(Xf, y)$coefficients
    return(list(fixef = setNames(cf, colnames(Xf)), gamma = setNames(numeric(0), character(0))))
  }
  D <- A[, sel_i, drop = FALSE]
  cf <- stats::lm.fit(cbind(Xf, D), y)$coefficients
  cf[is.na(cf)] <- 0  # guard against incidental rank deficiency
  list(fixef = setNames(cf[1:pf], colnames(Xf)),
       gamma = setNames(cf[(pf + 1):(pf + length(sel_i))], lev[sel_i]))
}

# Build augmented design [intercept | X | one-hot(cluster)] and fit.
scs_fit <- function(X, y, cluster, lambda, mu = NULL, eps = 1e-4, refit = TRUE,
                    intercept = TRUE, n_restart = 4L, seed = 1L) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  cluster <- factor(cluster); lev <- levels(cluster); K <- length(lev)
  n <- length(y); stopifnot(nrow(X) == n, length(cluster) == n)
  Xf <- if (intercept) cbind(`(Intercept)` = 1, X) else X
  p_fix <- ncol(Xf)
  A <- model.matrix(~ cluster - 1); colnames(A) <- lev
  Xfull <- cbind(Xf, A)
  if (is.null(mu)) mu <- eps * stats::median(colSums(Xfull^2))  # scale-relative ridge (D1)

  asym <- length(lambda) == 2
  if (asym) {
    lp <- as.integer(lambda[[1]]); ln <- as.integer(lambda[[2]])
    r <- scs_solve_dir_cpp(Xfull, as.double(y), p_fix, K, lp, ln, mu, as.integer(n_restart), as.integer(seed))
    sel_i <- c(r$selected_pos, r$selected_neg)
  } else {
    lam <- as.integer(lambda)
    r <- scs_solve_turbo_cpp(Xfull, as.double(y), p_fix, K, lam, mu, as.integer(n_restart), as.integer(seed))
    sel_i <- r$selected
  }
  sel <- lev[sel_i]

  if (refit) {                                   # D1: unshrunk LS on the selected support
    rf <- scs_refit_unshrunk(Xf, y, A, sel_i, lev)
    fixef <- rf$fixef
    gamma <- setNames(numeric(K), lev); gamma[names(rf$gamma)] <- rf$gamma
  } else {
    beta <- r$beta
    gamma <- beta[(p_fix + 1):(p_fix + K)]; names(gamma) <- lev
    fixef <- beta[1:p_fix]; names(fixef) <- colnames(Xf)
  }
  # tail assignment from the reported (refit) gamma so signs match what we report
  if (asym) { sel_pos <- names(which(gamma[sel] > 0)); sel_neg <- names(which(gamma[sel] < 0)) }
  else      { sel_pos <- sel_neg <- NULL }

  structure(list(
    call = match.call(), fixef = fixef, gamma = gamma,
    selected = sel, selected_pos = sel_pos, selected_neg = sel_neg,
    obj = r$obj, lambda = lambda, mu = mu, eps = eps, refit = refit, asymmetric = asym,
    levels = lev, p_fix = p_fix, intercept = intercept,
    data = list(X = X, y = y, cluster = cluster)
  ), class = "scs")
}

print.scs <- function(x, ...) {
  cat("Sparse Cluster Selection (L0 MIO LMM), non-Gurobi\n")
  cat(sprintf("  clusters K=%d, fixed effects p=%d, mu=%g, objective=%.6g\n",
              length(x$levels), x$p_fix, x$mu, x$obj))
  if (x$asymmetric) cat(sprintf("  asymmetric budgets: worst=%d, best=%d\n", x$lambda[[1]], x$lambda[[2]]))
  else cat(sprintf("  budget lambda=%d\n", x$lambda))
  cat(sprintf("  selected %d clusters: %s\n", length(x$selected),
              paste(utils::head(x$selected, 12), collapse = ", ")))
  invisible(x)
}
coef.scs <- function(object, ...) list(fixef = object$fixef, gamma = object$gamma)

# predict(): supervised CART map from cluster-level features to gamma, soft assignment (paper Sec 2.5).
# Default features = cluster means of X (uses X; pass Zsummary to add Z features).
predict.scs <- function(object, newdata, newcluster = NULL, type = c("response", "gamma"), ...) {
  type <- match.arg(type)
  newdata <- as.matrix(newdata)
  if (!is.null(newcluster)) {                      # known cluster -> use its gamma
    g <- object$gamma[as.character(newcluster)]; g[is.na(g)] <- 0
  } else {
    if (!requireNamespace("rpart", quietly = TRUE)) stop("install 'rpart' for new-cluster prediction")
    Xtr <- object$data$X; cl <- object$data$cluster
    fn <- paste0("f", seq_len(ncol(Xtr)))
    cmeans <- aggregate(Xtr, list(cl = cl), mean)
    feat <- cmeans[, -1, drop = FALSE]; colnames(feat) <- fn
    tr_df <- data.frame(y = object$gamma[as.character(cmeans$cl)], feat)
    tr <- rpart::rpart(y ~ ., data = tr_df, method = "anova",
                       control = rpart::rpart.control(maxdepth = max(1, ceiling(log2(length(object$levels))))))
    nd <- as.data.frame(newdata); colnames(nd) <- fn
    g <- as.numeric(predict(tr, newdata = nd))
  }
  if (type == "gamma") return(g)
  Xf <- if (object$intercept) cbind(1, newdata) else newdata
  as.numeric(Xf %*% object$fixef) + g
}

# stability(): subsample observations B times, refit, record per-cluster selection frequency.
# Meinshausen-Buhlmann (2010): for threshold pi_thr, E[# false flags] <= q^2 / ((2*pi_thr - 1) * K),
# where q = average #selected. Returns freq, the flagged set at pi_thr, and the implied error bound.
stability <- function(object, ...) UseMethod("stability")
stability.scs <- function(object, B = 100L, frac = 0.5, pi_thr = 0.6, seed = 1L, ...) {
  set.seed(seed)
  X <- object$data$X; y <- object$data$y; cl <- object$data$cluster
  lev <- object$levels; K <- length(lev)
  counts <- setNames(numeric(K), lev); qsum <- 0; used <- 0
  n <- length(y)
  for (b in 1:B) {
    idx <- sample.int(n, floor(frac * n))
    cb <- droplevels(cl[idx])
    if (nlevels(cb) < 2) next
    lamb <- if (object$asymmetric) object$lambda else min(object$lambda, nlevels(cb) - 1)
    fitb <- tryCatch(scs_fit(X[idx, , drop = FALSE], y[idx], cb, lamb,
                             mu = object$mu, intercept = object$intercept,
                             n_restart = 2L, seed = b), error = function(e) NULL)
    if (is.null(fitb)) next
    counts[fitb$selected] <- counts[fitb$selected] + 1
    qsum <- qsum + length(fitb$selected); used <- used + 1
  }
  freq <- counts / max(used, 1)
  q <- qsum / max(used, 1)
  Ev_bound <- if (pi_thr > 0.5) q^2 / ((2 * pi_thr - 1) * K) else NA_real_
  flagged <- names(sort(freq[freq >= pi_thr], decreasing = TRUE))
  structure(list(freq = sort(freq, decreasing = TRUE), flagged = flagged,
                 pi_thr = pi_thr, q = q, K = K, B_used = used,
                 expected_false_flags_bound = Ev_bound),
            class = "scs_stability")
}
print.scs_stability <- function(x, ...) {
  cat(sprintf("Stability selection: B_used=%d, avg selected q=%.2f, pi_thr=%.2f\n", x$B_used, x$q, x$pi_thr))
  cat(sprintf("  flagged (freq>=pi_thr): %s\n", paste(x$flagged, collapse = ", ")))
  cat(sprintf("  Meinshausen-Buhlmann bound on expected false flags: %.3f\n", x$expected_false_flags_bound))
  invisible(x)
}

# ============================================================================
# scs_path(): solve the symmetric L0 problem along the budget path lambda=0..lambda_max.
# Shared engine for deliverable C (audit budget) and D (cluster ranking). Returns the
# objective at each budget and the K x (lambda_max+1) membership matrix (cluster x budget).
# ============================================================================
scs_path <- function(X, y, cluster, lambda_max = NULL, mu = NULL, eps = 1e-4,
                     refit = TRUE, intercept = TRUE, n_restart = 4L, seed = 1L) {
  cluster <- factor(cluster); lev <- levels(cluster); K <- length(lev)
  if (is.null(lambda_max)) lambda_max <- min(K - 1L, max(2L, ceiling(K / 2)))
  budgets <- 0:lambda_max
  # fix mu once (from the full design) so the ridge is constant along the path -> comparable objectives
  if (is.null(mu)) {
    Xf <- if (intercept) cbind(1, as.matrix(X)) else as.matrix(X)
    A <- model.matrix(~ cluster - 1)
    mu <- eps * stats::median(colSums(cbind(Xf, A)^2))
  }
  member <- matrix(0L, nrow = K, ncol = length(budgets), dimnames = list(lev, paste0("l", budgets)))
  obj <- numeric(length(budgets))
  fits <- vector("list", length(budgets))
  for (i in seq_along(budgets)) {
    f <- scs_fit(X, y, cluster, lambda = budgets[i], mu = mu, refit = refit,
                 intercept = intercept, n_restart = n_restart, seed = seed)
    obj[i] <- f$obj
    if (length(f$selected)) member[f$selected, i] <- 1L
    fits[[i]] <- f
  }
  structure(list(budgets = budgets, obj = obj, member = member, levels = lev,
                 lambda_max = lambda_max, mu = mu, fits = fits),
            class = "scs_path")
}
print.scs_path <- function(x, ...) {
  cat(sprintf("SCS budget path: lambda=0..%d, K=%d, mu=%g\n", x$lambda_max, length(x$levels), x$mu))
  cat(sprintf("  objective: %s\n", paste(sprintf("%.4g", x$obj), collapse = " -> ")))
  invisible(x)
}

# ---- Deliverable C: lambda-as-audit-budget ---------------------------------
# Decision-theoretic framing: with capacity to review only lambda units, the L0 solve returns the
# lambda clusters that maximize detected deviation. detected(lambda)=obj(0)-obj(lambda) is the total
# deviation captured (RSS units); marginal(lambda) is the value of the lambda-th audit slot.
# Suggested budget = smallest lambda capturing >= target (default 90%) of the max detectable deviation.
# Two budget recommendations, both reported:
#  - target-based: smallest lambda capturing >= `target` of the max detectable deviation (default 90%);
#  - marginal-gain elbow: largest lambda whose review slot still buys >= `elbow_frac` of the FIRST slot's
#    gain (default 10%) -> the point of diminishing returns. (marginal generally decreasing but not
#    guaranteed monotone under local search, so we take the last budget above the elbow threshold.)
audit_budget <- function(object, ...) UseMethod("audit_budget")
audit_budget.scs_path <- function(object, target = 0.90, elbow_frac = 0.10, ...) {
  obj <- object$obj; b <- object$budgets
  detected <- obj[1] - obj                         # cumulative deviation captured (>=0, nondecreasing*)
  total <- detected[length(detected)]
  frac <- if (total > 0) detected / total else rep(0, length(detected))
  marginal <- c(NA_real_, -diff(obj))              # value of each additional audit slot (row i = budget b[i])
  first_slot <- if (length(marginal) >= 2) marginal[2] else NA_real_  # gain of the 1st reviewed cluster
  hit <- which(frac >= target)
  suggested <- if (length(hit)) b[min(hit)] else object$lambda_max
  worth <- which(!is.na(marginal) & marginal >= elbow_frac * first_slot)
  suggested_elbow <- if (length(worth)) b[max(worth)] else 0L
  structure(list(
    table = data.frame(lambda = b, objective = obj, detected = detected,
                       frac_of_max = frac, marginal_gain = marginal,
                       marginal_frac = marginal / first_slot),
    suggested_budget = suggested, suggested_budget_elbow = suggested_elbow,
    target = target, elbow_frac = elbow_frac, total_detectable = total),
    class = "scs_audit")
}
audit_budget.scs <- function(object, target = 0.90, elbow_frac = 0.10, lambda_max = NULL, ...) {
  p <- scs_path(object$data$X, object$data$y, object$data$cluster,
                lambda_max = lambda_max, mu = object$mu, refit = object$refit,
                intercept = object$intercept, ...)
  audit_budget(p, target = target, elbow_frac = elbow_frac)
}
print.scs_audit <- function(x, ...) {
  cat("Audit-budget profile (lambda = number of clusters you can review):\n")
  tb <- x$table; tb$objective <- signif(tb$objective, 5)
  tb$detected <- signif(tb$detected, 4); tb$frac_of_max <- round(tb$frac_of_max, 3)
  tb$marginal_gain <- signif(tb$marginal_gain, 4); tb$marginal_frac <- round(tb$marginal_frac, 3)
  print(tb, row.names = FALSE)
  cat(sprintf("Suggested budget  [target %.0f%% of detectable deviation]: lambda = %d\n",
              100 * x$target, x$suggested_budget))
  cat(sprintf("Suggested budget  [marginal-gain elbow, slot >= %.0f%% of 1st]: lambda = %d\n",
              100 * x$elbow_frac, x$suggested_budget_elbow))
  invisible(x)
}

# ---- Deliverable D: lambda-path cluster ranking (LARS-like) -----------------
# Rank clusters by the budget at which they FIRST enter the active set (entry_lambda): clusters that
# enter at small budgets are the most important. Ties broken by selection frequency across the path
# (times_selected) then |gamma| at first entry. A stable importance ordering, more robust than point gamma-hat.
cluster_ranking <- function(object, ...) UseMethod("cluster_ranking")
cluster_ranking.scs_path <- function(object, ...) {
  M <- object$member; b <- object$budgets; lev <- object$levels
  entry_lambda <- apply(M, 1, function(r) { w <- which(r == 1L); if (length(w)) b[min(w)] else NA_integer_ })
  times_selected <- rowSums(M)
  # |gamma| at the budget = entry_lambda (first appearance), 0 if never selected
  gamma_at_entry <- vapply(seq_along(lev), function(k) {
    el <- entry_lambda[k]; if (is.na(el)) return(0)
    abs(object$fits[[which(b == el)]]$gamma[lev[k]])
  }, numeric(1))
  df <- data.frame(cluster = lev, entry_lambda = entry_lambda,
                   times_selected = times_selected, gamma_at_entry = gamma_at_entry,
                   stringsAsFactors = FALSE)
  ord <- order(ifelse(is.na(df$entry_lambda), Inf, df$entry_lambda),
               -df$times_selected, -df$gamma_at_entry)
  df <- df[ord, ]; df$rank <- ifelse(is.na(df$entry_lambda), NA_integer_, seq_len(nrow(df)))
  rownames(df) <- NULL
  structure(list(ranking = df, lambda_max = object$lambda_max), class = "scs_ranking")
}
cluster_ranking.scs <- function(object, lambda_max = NULL, ...) {
  p <- scs_path(object$data$X, object$data$y, object$data$cluster,
                lambda_max = lambda_max, mu = object$mu, refit = object$refit,
                intercept = object$intercept, ...)
  cluster_ranking(p)
}
print.scs_ranking <- function(x, ...) {
  cat("Cluster importance ranking (lambda-path order of entry):\n")
  print(utils::head(x$ranking, 20), row.names = FALSE)
  invisible(x)
}

# ---- D1 support: epsilon-sensitivity of the selected set --------------------
# Report how the selection moves as the relative ridge eps varies; Jaccard vs the smallest-eps fit.
scs_eps_sensitivity <- function(X, y, cluster, lambda,
                                eps_grid = 10^(-6:-1), intercept = TRUE, ...) {
  fits <- lapply(eps_grid, function(e)
    scs_fit(X, y, cluster, lambda, eps = e, intercept = intercept, ...))
  sets <- lapply(fits, function(f) f$selected)
  ref <- sets[[1]]
  jacc <- vapply(sets, function(s) {
    u <- length(union(s, ref)); if (u == 0) 1 else length(intersect(s, ref)) / u
  }, numeric(1))
  data.frame(eps = eps_grid, mu = vapply(fits, `[[`, numeric(1), "mu"),
             n_selected = lengths(sets),
             jaccard_vs_smallest_eps = round(jacc, 3),
             selected = vapply(sets, paste, character(1), collapse = ","),
             stringsAsFactors = FALSE)
}
