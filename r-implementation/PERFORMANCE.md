# Performance roadmap

The released solver is the **validated** core: on a fixed `(data, lambda, seed)` it
reproduces the paper's selected sets bit-for-bit, and the R package wrapper is
bit-identical to the original `Rcpp::sourceCpp` path (verified). The optimizations
below are **deferred, not rejected**, for one reason:

> **Reproducibility gate.** The local search accepts a swap on a strict
> `gain > best + 1e-12` margin. Any change that reorders floating-point summation
> (grouped reductions, GEMM vs. dot-product, Cholesky-solve vs. explicit inverse)
> can move the last bits and, at a near-tie, **flip which cluster is selected** —
> which would change the flagged-hospital sets reported in the (frozen) manuscript.
> Each item below must therefore be landed behind a regression battery that asserts
> the selected set and objective are **unchanged across a grid of seeds and DGPs**
> before it is accepted. Items are tagged `[bit-safe]` (no FP reorder) or
> `[gated]` (needs the regression battery).

Findings are ordered by expected win. Source: efficiency council (Codex GPT-5.5,
high effort) + synthesis. Line references are to `src/src_scs.cpp` and `R/scs.R`.

## Tier 1 — architectural (largest wins)

1. **Grouped `FastCtx` construction; drop the dense one-hot.** `[gated]`
   `scs_fit()` builds `Xfull = [Xf | one-hot(cluster)]` (dense `n×(p+K)`) and
   `build_ctx()` forms `G`, `nvec`, `yvec` column-by-column (`F'a_k`), costing
   `O(n·p·K)` time and `O(nK)` memory. Pass `Xf` + integer cluster ids instead and
   accumulate in one row pass: `k=cid[i]; G.col(k)+=F.row(i); nvec[k]++; yvec[k]+=Y[i]`.
   → `O(n·p)` time, no dense dummy. Biggest single win, especially for large `K`.
   *Gate:* grouped sums are algebraically `F'A`, `A'Y`, `diag(A'A)` but reorder the
   additions — validate selected sets at tie-like supports.

2. **Path solver that builds the context once.** `[gated]`
   `scs_path()` refits via `scs_fit()` at every budget, rebuilding `model.matrix`,
   `cbind`, and `build_ctx` (`O(n·p² + n·K·p)`) `lambda_max+1` times. Add
   `scs_path_turbo_cpp(Xf, y, cid, lambda_max, mu, n_restart, seed)` that builds one
   `FastCtx` and sweeps budgets. *Gate:* must replay the exact
   `set.seed(seed); sample.int(K, lambda)` sequence per budget so RNG draws are
   identical. Same hot path benefits `audit_budget()` and `cluster_ranking()`.

3. **C++ batch worker for `stability()`.** `[gated]`
   The `B`-fit subsample loop constructs a full `scs` object and a dense dummy
   matrix per fit. A label-based C++ worker returning only selection counts removes
   `B×` R-object + dummy overhead (with item 1, each subsample context drops the `K`
   factor). *Gate:* preserve the current sampling order exactly — for `independent`
   sampling draw-then-fit-with-seed-`b`; for CPSS pass the precomputed `halves` as
   is. Do **not** reorder or precompute independent subsamples (changes the stream).

## Tier 2 — local, mostly bit-safe

4. **Selected-only dummy columns for the unshrunk refit.** `[bit-safe]`
   `scs_refit_unshrunk()` slices `D = A[, sel_i]` from the full one-hot. Build `D`
   (or its `n×lambda` block) directly from cluster ids — identical 0/1 matrix, so
   `lm.fit` input and output are unchanged; saves `O(nK)` construction. (Only pays
   off bundled with item 1, which removes the full `A`.)

5. **`partial_sort` for top-`lambda` warm starts.** `[bit-safe]`
   Turbo/fast warm starts `std::sort` all `K` score pairs; `std::partial_sort` to
   `lambda` gives the identical prefix under the same `(score, index)` comparator.
   `O(K log K) → O(K log lambda)`. Pure comparisons, no FP arithmetic.

6. **`fast_obj()` gain-only when coefficients aren't requested.** `[gated]`
   Every candidate eval reconstructs `beta = u - B·G_S·gamma` and the full
   objective. When `beta_out == nullptr`, return `(fixed0_2n - r·gamma)/(2n)` and
   skip the `Gg`, `B·Gg` work. Saves `O(p·lambda + p²)` per call (hot in
   `scs_solve_dir_cpp` and the legacy `fast_localsearch`). *Gate:* changes the
   floating summation path of the objective comparison.

## Tier 3 — constants / numerics

7. **Avoid the full `K×K` `W`.** `[gated]` Store `G`, `BG`; compute `W_ij` only for
   selected pairs and `diagW` columnwise. Context `O(p·K²)→O(p·K)`, memory
   `O(K²)→O(K)`. GEMM-vs-dot reorders sums.
8. **Cholesky solves instead of `FtF.inverse()`** (`build_ctx`, l.129). `[gated]`
   `LLT(FtF)` once; `u=llt.solve(FtY)`, `BG=llt.solve(G)`. Same big-O, lower
   constants and better conditioning; `.inverse()` shifts low bits.
9. **Remove allocation from the turbo candidate loop** (`turbo_localsearch`). Reuse
   work vectors; optionally precompute `Q*G`, `G'z` before the `j` scan.
   Allocation-only cleanup is `[bit-safe]`; batched GEMM is `[gated]`.
10. **Cache `XtX`, `XtY`, `YtY` in the legacy general-`q` solver** (`scs_solve_cpp`,
    `inner_obj`, `swap_pass`). Evaluate supports from submatrices instead of
    rebuilding `Xs` and multiplying by `X` each candidate: candidate eval
    `O(n·k² + n·p) → O(k³ + p·k)`, removing `n` from the inner loop. `[gated]`.

**Hard constraint (all tiers):** do not parallelize restarts, budgets, or stability
fits unless every `sample.int` draw is precomputed in the current order and
replayed — otherwise the RNG stream changes even when the objective does not.
