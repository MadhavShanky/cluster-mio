# scs.R — Sparse Cluster Selection: R API over the non-Gurobi C++ solver (src/src_scs.cpp).
# Idiom mirrors glm()/predict(): scs_fit() returns a fit object; predict()/stability() are layers on top.
#   fit  <- scs_fit(X, y, cluster, lambda = 5)                 # symmetric L0 budget
#   fit2 <- scs_fit(X, y, cluster, lambda = c(worst = 4, best = 2))  # asymmetric tail budgets
#   predict(fit, newdata = Xnew)                               # new-cluster prediction (CART soft map)
#   stability(fit, B = 100)                                    # per-cluster selection freq + FDR bound

#' clusterMIO: Sparse Cluster Selection for Linear Mixed Models
#'
#' An \eqn{\ell_0} (cardinality-constrained) selection of which cluster random
#' effects in a linear mixed model are held nonzero, with the rest set to exactly
#' zero. The solver profiles out the fixed effects once and exploits the disjoint
#' support of the one-hot cluster dummies, so each candidate swap has a closed-form
#' gain costing \eqn{O(p^2)} to evaluate, independent of the sample size. Depends
#' only on RcppEigen (no commercial optimizer). Scalable companion to the Julia
#' prototype (\code{clust_mio.jl}, \code{simul.jl}) in the repository root.
#'
#' @keywords internal
#' @useDynLib clusterMIO, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate coef cor lm.fit median model.matrix predict setNames
#' @importFrom utils head
"_PACKAGE"

# D1 (resolved 2026-06-19): selection uses a small SCALE-RELATIVE ridge mu = eps*median(diag(Xfull'Xfull));
# reported coefficients come from an UNPENALIZED LS refit on the selected support (refit=TRUE) so the
# paper's "unshrunk" claim is literally true. Pass an explicit numeric `mu` to override the relative rule
# (back-compat: mu=1 reproduces the old ridge-everything behaviour). `obj` stays the ridged selection
# objective; `fixef`/`gamma` are the refit. See scs_eps_sensitivity() for the ε-sensitivity report.

#' Unpenalized least-squares refit on a selected cluster support
#'
#' Internal helper for the D1 "unshrunk" reporting rule: given the augmented fixed
#' design and the selected cluster columns, return ordinary (ridge-free) LS
#' coefficients so the reported \code{fixef}/\code{gamma} are literally unshrunk.
#' @keywords internal
#' @noRd
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

#' Fit a sparse cluster-selection (L0 MIO LMM) model
#'
#' Builds the augmented design \code{[intercept | X | one-hot(cluster)]} and solves
#' the \eqn{\ell_0} cluster-selection problem: choose at most \code{lambda} clusters
#' whose random effect is nonzero (or, with a two-vector \code{lambda}, at most
#' \code{lambda[1]} worst-tail and \code{lambda[2]} best-tail clusters), minimizing
#' the ridge-conditioned profiled objective via best-improvement local search with
#' random restarts.
#'
#' @param X numeric matrix (or coercible) of fixed-effect covariates, \code{n} rows.
#' @param y numeric response vector of length \code{n}.
#' @param cluster factor (or coercible) of cluster labels, length \code{n}.
#' @param lambda integer L0 budget. A length-1 value gives the symmetric budget
#'   (at most \code{lambda} selected clusters); a length-2 value \code{c(worst, best)}
#'   gives asymmetric tail budgets (directional solver).
#' @param mu optional explicit ridge. Default \code{NULL} uses the scale-relative
#'   rule \code{mu = eps * median(colSums(Xfull^2))} (D1). \code{mu = 1} reproduces
#'   the legacy ridge-everything behaviour.
#' @param eps relative-ridge scale used when \code{mu} is \code{NULL} (default 1e-4).
#' @param refit logical; if \code{TRUE} (default) report an unpenalized LS refit on
#'   the selected support so \code{fixef}/\code{gamma} are unshrunk (the \code{obj}
#'   value is always the ridged selection objective).
#' @param intercept logical; prepend an intercept column to \code{X} (default TRUE).
#' @param n_restart integer number of local-search restarts (default 4).
#' @param seed integer RNG seed; the C++ restarts draw from R's RNG via
#'   \code{set.seed(seed)} so results are reproducible.
#' @return An object of class \code{"scs"}: a list with \code{fixef}, \code{gamma}
#'   (length-K, named by cluster), \code{selected} (and \code{selected_pos}/
#'   \code{selected_neg} when asymmetric), \code{obj}, and the inputs needed by
#'   \code{predict}, \code{stability}, \code{audit_budget} and \code{cluster_ranking}.
#' @examples
#' set.seed(1)
#' K <- 20; nk <- 15; n <- K * nk
#' cluster <- factor(rep(seq_len(K), each = nk))
#' X <- matrix(rnorm(n * 2), n, 2)
#' g <- numeric(K); g[c(3, 11)] <- c(4, -4)         # two deviating clusters
#' y <- X %*% c(1, -1) + g[as.integer(cluster)] + rnorm(n)
#' fit <- scs_fit(X, y, cluster, lambda = 2)
#' fit$selected
#' @export
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

#' @export
#' @noRd
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

#' @export
#' @noRd
coef.scs <- function(object, ...) list(fixef = object$fixef, gamma = object$gamma)

#' Predict cluster effects / responses for new data
#'
#' For a known cluster, returns its fitted \code{gamma}. For an unknown cluster,
#' applies the supervised CART map from cluster-level features (cluster means of
#' \code{X}) to \code{gamma} (paper Sec 2.5) for a soft assignment.
#' @param object an \code{"scs"} fit.
#' @param newdata numeric matrix of new covariate rows.
#' @param newcluster optional vector of known cluster labels (one per row of
#'   \code{newdata}); if supplied, the fitted \code{gamma} is used directly.
#' @param type \code{"response"} (default; linear predictor + cluster effect) or
#'   \code{"gamma"} (cluster effect only).
#' @param ... unused.
#' @return numeric vector of predictions.
#' @export
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

#' Stability selection over subsamples (Meinshausen-Buhlmann / CPSS)
#'
#' Subsamples the observations \code{B} times, refits the selector, and records the
#' per-cluster selection frequency, the set flagged at threshold \code{pi_thr}, and
#' the implied bound on the expected number of false flags.
#' @param object an \code{"scs"} fit.
#' @param ... passed to methods.
#' @return An object of class \code{"scs_stability"}.
#' @export
stability <- function(object, ...) UseMethod("stability")

#' @rdname stability
#' @param B integer number of subsample fits (for \code{pairing = "complementary"}
#'   this is \code{2 * #pairs} and must be a positive even number).
#' @param frac subsample fraction (forced to 0.5 under complementary pairing).
#' @param pi_thr selection-frequency threshold for flagging.
#' @param seed integer RNG seed.
#' @param scheme \code{"obs"} (observation-level subsampling, default) or
#'   \code{"cluster_strat"} (within-cluster subsampling so every candidate cluster
#'   is eligible in every subsample).
#' @param pairing \code{"independent"} (default; B independent half-subsamples) or
#'   \code{"complementary"} (Shah-Samworth complementary pairs: equal floor(n/2)
#'   halves; the CPSS sampling structure that makes the reference a genuine
#'   false-inclusion bound).
#' @param n_restart integer local-search restarts inside each subsample fit.
# scheme: "obs" = observation-level subsampling (default, original MB-style); "cluster_strat" =
#   cluster-stratified subsampling -- draw floor(frac*n_k) (>=1) observations WITHIN each cluster, so every
#   candidate cluster is eligible in every subsample. Removes the small-cluster selection-frequency cap
#   1-(1-frac)^{n_k} that otherwise makes pi_thr>=0.8 unreachable for small clusters, and supplies the
#   "all candidates eligible" precondition the MB/CPSS denominator K needs (see reviews/strengthen_synthesis.md).
# n_restart: local-search restarts inside each subsample fit (was hard-coded 2L); exposed to study
#   restart-jitter / robust-core sensitivity (sim/restart_robustcore_sim.R). Defaults preserve old behaviour.
# pairing: "independent" (default, original) draws B independent half-subsamples; "complementary" draws
#   B/2 Shah-Samworth COMPLEMENTARY PAIRS -- each pair splits a random permutation into two EQUAL halves A
#   and A^c of size floor(n/2) (one observation left out per pair when n is odd), so frac is forced to 1/2
#   and the two halves are an exact equal-size CPSS pair. Complementary pairs are the sampling structure the CPSS bound
#   (Shah & Samworth 2013, Thm 1) requires: it controls E|hatS_tau intersect L_theta| <= theta/(2tau-1) *
#   E|hatS_{1/2} intersect L_theta| with NO null-cluster exchangeability assumption, so under N subset
#   L_{q/K} the reference q^2/((2*pi_thr-1)*K) becomes a genuine false-inclusion bound (Supplement Prop CPSS).
#   The numeric reference is the SAME formula as MB; only the justification + the realized frequencies change.
#' @export
stability.scs <- function(object, B = 100L, frac = 0.5, pi_thr = 0.6, seed = 1L,
                          scheme = c("obs", "cluster_strat"),
                          pairing = c("independent", "complementary"), n_restart = 2L, ...) {
  scheme <- match.arg(scheme); pairing <- match.arg(pairing)
  # Fail loud on an unmatched named argument: a stale deployment whose signature lacked `pairing`/`scheme`
  # would otherwise swallow it into `...` and silently run a DIFFERENT procedure (this happened once on the
  # cluster). Refuse rather than downgrade.
  dots <- list(...)
  if (length(dots))
    stop("stability.scs(): unexpected argument(s) ", paste(sQuote(names(dots)), collapse = ", "),
         " -- refusing to run so a silently dropped 'pairing'/'scheme' cannot change the procedure.",
         call. = FALSE)
  if (pairing == "complementary" && scheme == "cluster_strat")
    stop("pairing='complementary' requires scheme='obs': within-stratum complements give unequal halves / ",
         "singleton strata, so they are not Shah-Samworth complementary pairs and the CPSS bound fails.",
         call. = FALSE)
  if (pairing == "complementary") {
    # CPSS fixes the subsample fraction at 1/2 (each pair = two equal halves) and needs an even number of
    # half-fits (B = 2 x #pairs); refuse anything that would silently reinterpret the design.
    if (!is.numeric(frac) || length(frac) != 1L || abs(frac - 0.5) > 1e-9)
      stop("pairing='complementary' fixes frac=0.5 (each pair is two equal halves); got frac=", frac, ".",
           call. = FALSE)
    if (B < 2L || B %% 2L != 0L)
      stop("pairing='complementary' requires a positive even B (B = number of half-fits = 2 x #pairs); got B=",
           B, ".", call. = FALSE)
  }
  set.seed(seed)
  X <- object$data$X; y <- object$data$y; cl <- object$data$cluster
  lev <- object$levels; K <- length(lev)
  counts <- setNames(numeric(K), lev); qsum <- 0; used <- 0
  n <- length(y)
  by_cluster <- if (scheme == "cluster_strat") split(seq_len(n), cl) else NULL
  # fit the selector on one subsample (index set idx) and accumulate per-cluster selection counts;
  # returns TRUE if the fit was usable (>=2 clusters and solver succeeded), FALSE otherwise.
  do_fit <- function(idx, sd) {
    cb <- droplevels(cl[idx])
    if (nlevels(cb) < 2) return(FALSE)
    lamb <- if (object$asymmetric) object$lambda else min(object$lambda, nlevels(cb) - 1)
    fitb <- tryCatch(scs_fit(X[idx, , drop = FALSE], y[idx], cb, lamb,
                             mu = object$mu, intercept = object$intercept,
                             n_restart = as.integer(n_restart), seed = sd), error = function(e) NULL)
    if (is.null(fitb)) return(FALSE)
    counts[fitb$selected] <<- counts[fitb$selected] + 1
    qsum <<- qsum + length(fitb$selected); used <<- used + 1
    TRUE
  }
  # draw the 'A' half of a subsample (observation-level or cluster-stratified)
  draw_half <- function(half_frac) {
    if (scheme == "cluster_strat") {                      # within-cluster subsample, keep every cluster
      unlist(lapply(by_cluster, function(ii)
        if (length(ii) <= 1L) ii else sample(ii, max(1L, floor(half_frac * length(ii))))), use.names = FALSE)
    } else {
      sample.int(n, floor(half_frac * n))
    }
  }
  if (pairing == "complementary") {
    # Shah-Samworth complementary pairs. Precompute ALL half-sample index sets under the master seed
    # BEFORE any fitting: scs_fit() calls R's set.seed() internally (the C++ restarts draw from R's RNG),
    # so interleaving fits would drive every subsequent subsample draw off the fit seeds rather than `seed`.
    # Each pair is a random permutation split into two EQUAL halves of size floor(n/2); for odd n the single
    # leftover observation is left out of that pair, so A and A^c are a valid equal-size CPSS pair.
    m <- floor(n / 2L)
    npairs <- max(1L, B %/% 2L)
    halves <- vector("list", 2L * npairs)
    for (p in seq_len(npairs)) {
      perm <- sample.int(n, n)
      halves[[2L * p - 1L]] <- perm[1:m]
      halves[[2L * p]]      <- perm[(m + 1L):(2L * m)]    # disjoint from, and equal size to, its partner
    }
    planned <- length(halves)                             # CPSS denominator = planned halves (a failed fit
    for (j in seq_along(halves)) do_fit(halves[[j]], j)   #   counts as an empty selection, still in denom)
    denom <- planned
  } else {
    for (b in 1:B) do_fit(draw_half(frac), b)             # original independent subsampling (unchanged)
    denom <- used
  }
  freq <- counts / max(denom, 1)
  q <- qsum / max(denom, 1)
  Ev_bound <- if (pi_thr > 0.5) q^2 / ((2 * pi_thr - 1) * K) else NA_real_
  flagged <- names(sort(freq[freq >= pi_thr], decreasing = TRUE))
  structure(list(freq = sort(freq, decreasing = TRUE), flagged = flagged,
                 pi_thr = pi_thr, q = q, K = K, B_used = used, B_planned = denom,
                 scheme = scheme, pairing = pairing,
                 n_restart = n_restart, expected_false_flags_bound = Ev_bound),
            class = "scs_stability")
}

#' @export
#' @noRd
print.scs_stability <- function(x, ...) {
  cat(sprintf("Stability selection: B_used=%d, avg selected q=%.2f, pi_thr=%.2f\n", x$B_used, x$q, x$pi_thr))
  cat(sprintf("  flagged (freq>=pi_thr): %s\n", paste(x$flagged, collapse = ", ")))
  ref_lab <- if (!is.null(x$pairing) && x$pairing == "complementary") "CPSS (complementary-pairs)" else "Meinshausen-Buhlmann"
  cat(sprintf("  %s reference on expected false flags: %.3f\n", ref_lab, x$expected_false_flags_bound))
  invisible(x)
}

# ============================================================================
# scs_path(): solve the symmetric L0 problem along the budget path lambda=0..lambda_max.
# Shared engine for deliverable C (audit budget) and D (cluster ranking). Returns the
# objective at each budget and the K x (lambda_max+1) membership matrix (cluster x budget).
# ============================================================================

#' Solve the symmetric L0 problem along the budget path
#'
#' Fits \code{scs_fit} at every budget \code{lambda = 0, 1, ..., lambda_max} with a
#' single fixed ridge \code{mu} (so the objectives are comparable across budgets).
#' Shared engine for \code{audit_budget} and \code{cluster_ranking}.
#' @param X,y,cluster as in \code{scs_fit}.
#' @param lambda_max maximum budget (default \code{min(K-1, max(2, ceiling(K/2)))}).
#' @param mu,eps,refit,intercept,n_restart,seed as in \code{scs_fit} (mu fixed once
#'   from the full design so the ridge is constant along the path).
#' @return An object of class \code{"scs_path"} with the per-budget objective and the
#'   \code{K x (lambda_max+1)} membership matrix.
#' @export
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

#' @export
#' @noRd
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

#' Audit-budget profile (deliverable C)
#'
#' Frames \code{lambda} as the number of clusters one has capacity to review and
#' reports how much detectable deviation each additional audit slot captures, with
#' two suggested budgets (a target-coverage budget and a marginal-gain elbow).
#' @param object an \code{"scs_path"} or \code{"scs"} object.
#' @param ... passed to methods.
#' @return An object of class \code{"scs_audit"}.
#' @export
audit_budget <- function(object, ...) UseMethod("audit_budget")

#' @rdname audit_budget
#' @param target coverage target (default 0.90 of max detectable deviation).
#' @param elbow_frac marginal-gain elbow threshold (default 0.10 of the first slot).
#' @export
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

#' @rdname audit_budget
#' @param lambda_max maximum budget for the internally-computed path.
#' @export
audit_budget.scs <- function(object, target = 0.90, elbow_frac = 0.10, lambda_max = NULL, ...) {
  p <- scs_path(object$data$X, object$data$y, object$data$cluster,
                lambda_max = lambda_max, mu = object$mu, refit = object$refit,
                intercept = object$intercept, ...)
  audit_budget(p, target = target, elbow_frac = elbow_frac)
}

#' @export
#' @noRd
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

#' Cluster importance ranking along the budget path (deliverable D)
#'
#' Ranks clusters by the budget at which they first enter the active set (earlier
#' entry = more important), breaking ties by selection frequency across the path and
#' then by \code{|gamma|} at first entry.
#' @param object an \code{"scs_path"} or \code{"scs"} object.
#' @param ... passed to methods.
#' @return An object of class \code{"scs_ranking"}.
#' @export
cluster_ranking <- function(object, ...) UseMethod("cluster_ranking")

#' @rdname cluster_ranking
#' @export
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

#' @rdname cluster_ranking
#' @param lambda_max maximum budget for the internally-computed path.
#' @export
cluster_ranking.scs <- function(object, lambda_max = NULL, ...) {
  p <- scs_path(object$data$X, object$data$y, object$data$cluster,
                lambda_max = lambda_max, mu = object$mu, refit = object$refit,
                intercept = object$intercept, ...)
  cluster_ranking(p)
}

#' @export
#' @noRd
print.scs_ranking <- function(x, ...) {
  cat("Cluster importance ranking (lambda-path order of entry):\n")
  print(utils::head(x$ranking, 20), row.names = FALSE)
  invisible(x)
}

# ---- D1 support: epsilon-sensitivity of the selected set --------------------
# Report how the selection moves as the relative ridge eps varies; Jaccard vs the smallest-eps fit.

#' Epsilon-sensitivity of the selected set (D1 support)
#'
#' Reports how the selected cluster set moves as the relative-ridge scale \code{eps}
#' varies over a grid, via the Jaccard overlap against the smallest-\code{eps} fit.
#' @param X,y,cluster,lambda as in \code{scs_fit}.
#' @param eps_grid numeric grid of relative-ridge scales (default \code{10^(-6:-1)}).
#' @param intercept logical; prepend an intercept (default TRUE).
#' @param ... passed to \code{scs_fit}.
#' @return A data frame: \code{eps}, realized \code{mu}, \code{n_selected}, the
#'   Jaccard overlap vs the smallest-\code{eps} fit, and the selected set.
#' @export
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
