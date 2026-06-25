// Sparse Cluster Selection (SCS) — non-Gurobi, pure-C++ solver.
// Reproduces the cluster-mio profiled objective EXACTLY:
//   support ind = fixed cols (always in) + selected cluster cols
//   alpha = Y - Xs (mu I + Xs'Xs)^{-1} Xs'Y      (ridge mu on all selected cols, as in clust_mio.jl)
//   c(s)  = <Y, alpha> / (2n)
//   grad_j = -(x_j' alpha)^2 / (2 n mu)   over ALL columns
//   beta_support = (mu I + Xs'Xs)^{-1} Xs'Y  ( == (1/mu) Xs' alpha )
// Engine: gradient-ranking warm start + best-improvement swap local search with random restarts (L1).
// (OA branch-and-cut certificate L2 added separately.)
//
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::LLT<MatrixXd> LLTd;

// ---- inner solve on a given support (0-based column indices) ----
static double inner_obj(const MatrixXd& X, const VectorXd& Y,
                        const std::vector<int>& ind, double mu,
                        VectorXd* grad_out = nullptr,
                        VectorXd* coef_out = nullptr) {
  const int n = X.rows();
  const int k = (int)ind.size();
  MatrixXd Xs(n, k);
  for (int j = 0; j < k; ++j) Xs.col(j) = X.col(ind[j]);
  MatrixXd G = Xs.transpose() * Xs;
  G.diagonal().array() += mu;
  VectorXd XtY = Xs.transpose() * Y;
  LLTd llt(G);
  VectorXd coef = llt.solve(XtY);
  VectorXd alpha = Y - Xs * coef;
  double obj = Y.dot(alpha) / (2.0 * n);
  if (grad_out) {
    VectorXd t = X.transpose() * alpha;          // length p
    *grad_out = -(t.array().square()) / (2.0 * n * mu);
  }
  if (coef_out) *coef_out = coef;
  return obj;
}

// build support = [0..p_fix-1] + active cluster cols (global indices)
static std::vector<int> build_support(int p_fix, int K, int q,
                                      const std::vector<std::vector<char>>& active) {
  std::vector<int> ind;
  ind.reserve(p_fix + q * K);
  for (int j = 0; j < p_fix; ++j) ind.push_back(j);
  for (int l = 0; l < q; ++l)
    for (int kk = 0; kk < K; ++kk)
      if (active[l][kk]) ind.push_back(p_fix + l * K + kk);
  return ind;
}

// one best-improvement swap pass over all blocks; returns true if improved
static bool swap_pass(const MatrixXd& X, const VectorXd& Y, int p_fix, int K, int q,
                      double mu, std::vector<std::vector<char>>& active, double& cur_obj) {
  bool any = false;
  for (int l = 0; l < q; ++l) {
    bool improved = true;
    while (improved) {
      improved = false;
      // candidate in/out lists for block l
      std::vector<int> in_idx, out_idx;
      for (int kk = 0; kk < K; ++kk) (active[l][kk] ? in_idx : out_idx).push_back(kk);
      double best = cur_obj; int bi = -1, bo = -1;
      for (size_t a = 0; a < in_idx.size(); ++a) {
        for (size_t b = 0; b < out_idx.size(); ++b) {
          active[l][in_idx[a]] = 0; active[l][out_idx[b]] = 1;
          std::vector<int> ind = build_support(p_fix, K, q, active);
          double o = inner_obj(X, Y, ind, mu);
          active[l][in_idx[a]] = 1; active[l][out_idx[b]] = 0;  // revert
          if (o < best - 1e-12) { best = o; bi = in_idx[a]; bo = out_idx[b]; }
        }
      }
      if (bi >= 0) {
        active[l][bi] = 0; active[l][bo] = 1;
        cur_obj = best; improved = true; any = true;
      }
    }
  }
  return any;
}

// gradient-ranking warm start for one block configuration
static void rank_warmstart(const MatrixXd& X, const VectorXd& Y, int p_fix, int K, int q,
                           const std::vector<int>& lambda, double mu,
                           std::vector<std::vector<char>>& active) {
  // fixed-only fit residual gradient
  std::vector<int> ind0;
  for (int j = 0; j < p_fix; ++j) ind0.push_back(j);
  VectorXd grad;
  inner_obj(X, Y, ind0, mu, &grad);
  for (int l = 0; l < q; ++l) {
    std::vector<std::pair<double,int>> sc(K);
    for (int kk = 0; kk < K; ++kk) sc[kk] = { grad(p_fix + l * K + kk), kk };  // most negative = best
    std::sort(sc.begin(), sc.end());
    std::fill(active[l].begin(), active[l].end(), 0);
    for (int r = 0; r < lambda[l] && r < K; ++r) active[l][sc[r].second] = 1;
  }
}

// ================= FAST PATH (q=1): Schur-complement elimination of the fixed block =================
// Exact same profiled objective c(S), but each support evaluation is a |S|x|S| solve assembled from
// precomputed scalars (cost independent of n after an O(n) precompute). Scales to large K.
struct FastCtx {
  int n, p_fix, K;
  double mu, YtY;
  Eigen::MatrixXd B;        // (F'F+mu I)^{-1}, p_fix x p_fix
  Eigen::MatrixXd G;        // p_fix x K, columns g_k = F' a_k
  Eigen::VectorXd FtY, u;   // p_fix
  Eigen::MatrixXd W;        // K x K, W = G' B G
  Eigen::MatrixXd BG;       // p_fix x K, = B * G   (columns B g_k)
  Eigen::VectorXd nvec, yvec, gu;  // length K: n_k, sum Y_k, g_k' u
  Eigen::VectorXd rhs;      // length K: yvec - gu  (also: cluster warm-start score = rhs^2)
  Eigen::VectorXd diagW, dvec; // length K: W_kk ; d_k = n_k + mu - W_kk
  double fixed0_2n;         // fixed-only objective * 2n = YtY - FtY' u ; gain(S)=r_S'H_S^{-1}r_S, c*2n = fixed0_2n - gain
};

static FastCtx build_ctx(const MatrixXd& X, const VectorXd& Y, int p_fix, int K, double mu) {
  FastCtx c; c.n = X.rows(); c.p_fix = p_fix; c.K = K; c.mu = mu;
  MatrixXd F = X.leftCols(p_fix);
  // Small ridge mu on the FULL augmented design (fixed + selected), as in best-subset with ridge
  // (Bertsimas-King-Mazumder). This conditions the selection objective; the reported coefficients come
  // from an unpenalized refit, so the estimand is unchanged. Removing the fixed-block ridge is NOT
  // cosmetic -- it shifts the local-search selection and changes the stability watch-list (NY 14 -> 4).
  MatrixXd FtF = F.transpose() * F; FtF.diagonal().array() += mu;
  c.B = FtF.inverse();
  c.FtY = F.transpose() * Y;
  c.u = c.B * c.FtY;
  c.YtY = Y.dot(Y);
  c.G = MatrixXd::Zero(p_fix, K);
  c.nvec = VectorXd::Zero(K); c.yvec = VectorXd::Zero(K);
  for (int k = 0; k < K; ++k) {
    VectorXd a = X.col(p_fix + k);           // one-hot cluster column
    c.G.col(k) = F.transpose() * a;          // g_k
    c.nvec(k) = a.sum();
    c.yvec(k) = a.dot(Y);
  }
  c.BG = c.B * c.G;
  c.W = c.G.transpose() * c.BG;
  c.gu = c.G.transpose() * c.u;
  c.rhs = c.yvec - c.gu;
  c.diagW = c.W.diagonal();
  c.dvec = c.nvec.array() + mu - c.diagW.array();
  c.fixed0_2n = c.YtY - c.FtY.dot(c.u);
  return c;
}

// forward decl (defined below)
static double fast_obj(const FastCtx& c, const std::vector<int>& sel,
                       VectorXd* gamma_out, VectorXd* beta_out);

// build H_S for a cluster set sel (0-based)
static MatrixXd build_H(const FastCtx& c, const std::vector<int>& sel) {
  const int s = (int)sel.size();
  MatrixXd H(s, s);
  for (int i = 0; i < s; ++i) {
    H(i, i) = c.dvec(sel[i]);
    for (int j = i + 1; j < s; ++j) { double w = -c.W(sel[i], sel[j]); H(i, j) = w; H(j, i) = w; }
  }
  return H;
}
static double gain_of(const FastCtx& c, const std::vector<int>& sel) {
  if (sel.empty()) return 0.0;
  const int s = (int)sel.size();
  VectorXd r(s); for (int i = 0; i < s; ++i) r(i) = c.rhs(sel[i]);
  MatrixXd H = build_H(c, sel);
  Eigen::LLT<MatrixXd> llt(H);
  return r.dot(llt.solve(r));
}

// TURBO best-improvement swap local search using closed-form Schur swap gains.
// Maximizes gain(S) = r_S' H_S^{-1} r_S  (objective c*2n = fixed0_2n - gain). Pure swaps: |S| stays = lambda.
//
// Turbo #3 (rank-1 downdate + closed-form leave corrections): factor H_S ONCE per pass to get M = H_S^{-1}
// and the full-set precomputes (w_full = M r_S, R = BG_S M, Q_S = R BG_S', z_full = R r_S). Each leave-
// candidate core C = S\{a} (delete index a) is then a rank-1 correction of these, so the quantities the
// candidate scan needs reduce to closed forms with NO per-a m-sized solve:
//   gain_C = gain_S - w_full(a)^2 / M_aa                                         (O(1))
//   p_a    = R.col(a) - bg_a M_aa  (= (BG)_C M_{C,a})                            (O(p_fix))
//   z_C    = z_full - bg_a w_full(a) - p_a w_full(a)/M_aa                         (O(p_fix))
//   Q_C    = Q_S - bg_a p_a' - p_a bg_a' - M_aa bg_a bg_a' - p_a p_a'/M_aa        (O(p_fix^2))
// Collapses the pass from O(lambda^4) (refactor per a) to O(lambda^3 + lambda*K*p_fix^2). Algebraically
// identical to the explicit-H_C^{-1} gains -> still EXACT (validated == brute force at machine eps).
static void turbo_localsearch(const FastCtx& c, int lambda, std::vector<int>& sel, double& gain) {
  const int pf = c.p_fix;
  bool improved = true;
  while (improved) {
    improved = false;
    const int m = (int)sel.size();
    if (m == 0) break;                            // lambda=0: no swaps
    double best_gain = gain; int best_a = -1, best_j = -1;
    std::vector<char> inset(c.K, 0); for (int x : sel) inset[x] = 1;

    // --- full-set precompute (once per accepted move) ---
    MatrixXd H = build_H(c, sel);
    Eigen::LLT<MatrixXd> llt(H);
    MatrixXd M = llt.solve(MatrixXd::Identity(m, m));     // H_S^{-1}, O(m^3)
    VectorXd rS(m); MatrixXd BGS(pf, m);
    for (int i = 0; i < m; ++i) { rS(i) = c.rhs(sel[i]); BGS.col(i) = c.BG.col(sel[i]); }
    VectorXd w_full = M * rS;                              // O(m^2)
    double   gain_S = rS.dot(w_full);
    MatrixXd R = BGS * M;                                  // p_fix x m, O(p_fix m^2)
    MatrixXd Q_S = R * BGS.transpose();                    // p_fix x p_fix, O(p_fix^2 m)
    VectorXd z_full = R * rS;                              // p_fix, O(p_fix m)

    for (int a = 0; a < m; ++a) {
      double gain_C; VectorXd z; MatrixXd Q;
      if (m == 1) { gain_C = 0.0; z = VectorXd::Zero(pf); Q = MatrixXd::Zero(pf, pf); }
      else {
        const double Maa = M(a, a), wa = w_full(a);
        const VectorXd bg_a = BGS.col(a);
        VectorXd p_a = R.col(a) - bg_a * Maa;              // (BG)_C M_{C,a}
        gain_C = gain_S - wa * wa / Maa;
        z = z_full - bg_a * wa - p_a * (wa / Maa);
        Q = Q_S - bg_a * p_a.transpose() - p_a * bg_a.transpose()
              - Maa * (bg_a * bg_a.transpose()) - (p_a * p_a.transpose()) / Maa;
      }
      // candidate adds j (j not in current sel)
      for (int j = 0; j < c.K; ++j) {
        if (inset[j]) continue;
        const VectorXd& gj = c.G.col(j);
        double hMh = gj.dot(Q * gj);            // h_j' H_C^{-1} h_j ; O(p_fix^2)
        double s_j = c.dvec(j) - hMh;
        if (s_j <= 1e-12) continue;             // numerical guard
        double hw = -gj.dot(z);                 // h_j' w_C  (h_j = -(BG)_C' g_j)
        double num = c.rhs(j) - hw;
        double g_new = gain_C + num * num / s_j;
        if (g_new > best_gain + 1e-12) { best_gain = g_new; best_a = a; best_j = j; }
      }
    }
    if (best_a >= 0) { sel[best_a] = best_j; gain = best_gain; improved = true; }
  }
}

// [[Rcpp::export]]
List scs_solve_turbo_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                         const Eigen::Map<Eigen::VectorXd> Y,
                         int p_fix, int K, int lambda, double mu,
                         int n_restart = 4, int seed = 1) {
  FastCtx c = build_ctx(X, Y, p_fix, K, mu);
  // warm start: top-lambda by rhs^2
  std::vector<std::pair<double,int>> sc(K);
  for (int k = 0; k < K; ++k) sc[k] = { -(c.rhs(k) * c.rhs(k)), k };
  std::sort(sc.begin(), sc.end());
  std::vector<int> sel; for (int r = 0; r < lambda && r < K; ++r) sel.push_back(sc[r].second);
  double gain = gain_of(c, sel);
  turbo_localsearch(c, lambda, sel, gain);
  std::vector<int> best = sel; double best_gain = gain;

  Function set_seed("set.seed"); set_seed(seed);
  Function sample_int("sample.int");
  for (int r = 1; r < n_restart; ++r) {
    IntegerVector pick = sample_int(K, lambda);
    std::vector<int> s2; for (int t = 0; t < pick.size(); ++t) s2.push_back(pick[t] - 1);
    double g2 = gain_of(c, s2);
    turbo_localsearch(c, lambda, s2, g2);
    if (g2 > best_gain + 1e-12) { best_gain = g2; best = s2; }
  }
  std::sort(best.begin(), best.end());
  // recover beta via fast_obj (exact)
  VectorXd gamma, beta_fix;
  double obj = fast_obj(c, best, &gamma, &beta_fix);
  VectorXd beta = VectorXd::Zero(p_fix + K);
  beta.head(p_fix) = beta_fix;
  for (size_t i = 0; i < best.size(); ++i) beta(p_fix + best[i]) = gamma(i);
  IntegerVector sel1(best.size()); for (size_t i = 0; i < best.size(); ++i) sel1[i] = best[i] + 1;
  return List::create(_["obj"] = obj, _["gain"] = best_gain,
                      _["obj_from_gain"] = (c.fixed0_2n - best_gain) / (2.0 * c.n),
                      _["beta"] = beta, _["selected"] = sel1);
}

// objective c(S) for a selected cluster set sel (0-based), via the reduced lambda x lambda system
static double fast_obj(const FastCtx& c, const std::vector<int>& sel,
                       VectorXd* gamma_out = nullptr, VectorXd* beta_out = nullptr) {
  const int s = (int)sel.size();
  if (s == 0) {
    double obj2n = c.YtY - c.FtY.dot(c.u);
    if (beta_out) *beta_out = c.u;
    if (gamma_out) *gamma_out = VectorXd();
    return obj2n / (2.0 * c.n);
  }
  MatrixXd H(s, s); VectorXd r(s);
  for (int i = 0; i < s; ++i) {
    r(i) = c.rhs(sel[i]);
    H(i, i) = (c.nvec(sel[i]) + c.mu) - c.W(sel[i], sel[i]);
    for (int j = i + 1; j < s; ++j) { double w = -c.W(sel[i], sel[j]); H(i, j) = w; H(j, i) = w; }
  }
  Eigen::LLT<MatrixXd> llt(H);
  VectorXd gamma = llt.solve(r);
  // G_S gamma  (p_fix)
  VectorXd Gg = VectorXd::Zero(c.p_fix);
  for (int i = 0; i < s; ++i) Gg += c.G.col(sel[i]) * gamma(i);
  VectorXd beta = c.u - c.B * Gg;
  double obj2n = c.YtY - c.FtY.dot(beta);
  for (int i = 0; i < s; ++i) obj2n -= c.yvec(sel[i]) * gamma(i);
  if (gamma_out) *gamma_out = gamma;
  if (beta_out) *beta_out = beta;
  return obj2n / (2.0 * c.n);
}

static void fast_localsearch(const FastCtx& c, int lambda, std::vector<int>& sel, double& obj) {
  bool improved = true;
  while (improved) {
    improved = false;
    std::vector<char> inset(c.K, 0); for (int x : sel) inset[x] = 1;
    std::vector<int> out; for (int k = 0; k < c.K; ++k) if (!inset[k]) out.push_back(k);
    double best = obj; int bi = -1, bo = -1;
    for (size_t a = 0; a < sel.size(); ++a) {
      int saved = sel[a];
      for (size_t b = 0; b < out.size(); ++b) {
        sel[a] = out[b];
        double o = fast_obj(c, sel);
        if (o < best - 1e-12) { best = o; bi = (int)a; bo = out[b]; }
      }
      sel[a] = saved;
    }
    if (bi >= 0) { sel[bi] = bo; obj = best; improved = true; }
  }
}

// [[Rcpp::export]]
List scs_solve_fast_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                        const Eigen::Map<Eigen::VectorXd> Y,
                        int p_fix, int K, int lambda, double mu,
                        int n_restart = 4, int seed = 1) {
  FastCtx c = build_ctx(X, Y, p_fix, K, mu);

  // warm start: top-lambda clusters by rhs^2 (= (a_k' residual_at_fixed_fit)^2)
  std::vector<std::pair<double,int>> sc(K);
  for (int k = 0; k < K; ++k) sc[k] = { -(c.rhs(k) * c.rhs(k)), k };
  std::sort(sc.begin(), sc.end());
  std::vector<int> sel; for (int r = 0; r < lambda && r < K; ++r) sel.push_back(sc[r].second);
  double obj = fast_obj(c, sel);
  fast_localsearch(c, lambda, sel, obj);
  std::vector<int> best_sel = sel; double best_obj = obj;

  Function set_seed("set.seed"); set_seed(seed);
  Function sample_int("sample.int");
  for (int r = 1; r < n_restart; ++r) {
    IntegerVector pick = sample_int(K, lambda);
    std::vector<int> s2; for (int t = 0; t < pick.size(); ++t) s2.push_back(pick[t] - 1);
    double o2 = fast_obj(c, s2);
    fast_localsearch(c, lambda, s2, o2);
    if (o2 < best_obj - 1e-12) { best_obj = o2; best_sel = s2; }
  }

  std::sort(best_sel.begin(), best_sel.end());
  VectorXd gamma, beta_fix;
  double fobj = fast_obj(c, best_sel, &gamma, &beta_fix);
  VectorXd beta = VectorXd::Zero(p_fix + K);
  beta.head(p_fix) = beta_fix;
  for (size_t i = 0; i < best_sel.size(); ++i) beta(p_fix + best_sel[i]) = gamma(i);
  IntegerVector sel1(best_sel.size()); for (size_t i = 0; i < best_sel.size(); ++i) sel1[i] = best_sel[i] + 1;
  return List::create(_["obj"] = fobj, _["beta"] = beta, _["selected"] = sel1);
}

// ================= ASYMMETRIC / DIRECTIONAL selection =================
// Tail-count budgets on the FITTED deviations: choose a support S so that the unconstrained fitted
// gamma has at most lambda_pos positive and lambda_neg negative entries, minimizing c(S).
// (Definition: budgets on each tail of the fitted random effects -> "flag <=L+ worst & <=L- best".)
// Engine: forward-greedy respecting tail caps + best-improvement swaps + random restarts. Uses the
// fast (Schur) inner solve so we get fitted gamma each eval. Exact-global is not claimed for the
// asymmetric variant (the symmetric solver is the certified one); validated vs brute force below.
static bool dir_feasible(const VectorXd& g, int lpos, int lneg, double tol, int& np, int& nn) {
  np = 0; nn = 0;
  for (int i = 0; i < g.size(); ++i) { if (g(i) > tol) ++np; else if (g(i) < -tol) ++nn; }
  return np <= lpos && nn <= lneg;
}

// [[Rcpp::export]]
List scs_solve_dir_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                       const Eigen::Map<Eigen::VectorXd> Y,
                       int p_fix, int K, int lambda_pos, int lambda_neg, double mu,
                       int n_restart = 4, int seed = 1) {
  FastCtx c = build_ctx(X, Y, p_fix, K, mu);
  const int cap = std::min(K, lambda_pos + lambda_neg);
  const double tol = 1e-9;
  auto obj_of = [&](const std::vector<int>& sel, VectorXd& g) {
    VectorXd b; return fast_obj(c, sel, &g, &b);
  };
  auto polish = [&](std::vector<int> sel, double obj) {
    // best-improvement swaps + add/drop respecting tail caps
    bool improved = true;
    while (improved) {
      improved = false;
      std::vector<char> inset(K, 0); for (int x : sel) inset[x] = 1;
      double best = obj; std::vector<int> bestsel = sel; bool found = false;
      // swaps
      for (size_t a = 0; a < sel.size(); ++a) for (int j = 0; j < K; ++j) {
        if (inset[j]) continue;
        std::vector<int> cand = sel; cand[a] = j;
        VectorXd g; double o = obj_of(cand, g); int np, nn;
        if (o < best - 1e-12 && dir_feasible(g, lambda_pos, lambda_neg, tol, np, nn)) { best = o; bestsel = cand; found = true; }
      }
      // drops
      for (size_t a = 0; a < sel.size(); ++a) {
        std::vector<int> cand; for (size_t t = 0; t < sel.size(); ++t) if (t != a) cand.push_back(sel[t]);
        VectorXd g; double o = obj_of(cand, g); int np, nn;
        if (o < best - 1e-12 && dir_feasible(g, lambda_pos, lambda_neg, tol, np, nn)) { best = o; bestsel = cand; found = true; }
      }
      // adds (if room)
      if ((int)sel.size() < cap) for (int j = 0; j < K; ++j) {
        if (inset[j]) continue;
        std::vector<int> cand = sel; cand.push_back(j);
        VectorXd g; double o = obj_of(cand, g); int np, nn;
        if (o < best - 1e-12 && dir_feasible(g, lambda_pos, lambda_neg, tol, np, nn)) { best = o; bestsel = cand; found = true; }
      }
      if (found) { sel = bestsel; obj = best; improved = true; }
    }
    return std::make_pair(obj, sel);
  };

  // forward greedy start
  std::vector<int> sel; double obj;
  { VectorXd g0; obj = obj_of(sel, g0); }
  bool grow = true;
  while ((int)sel.size() < cap && grow) {
    grow = false; std::vector<char> inset(K, 0); for (int x : sel) inset[x] = 1;
    double best = obj; int bj = -1;
    for (int j = 0; j < K; ++j) {
      if (inset[j]) continue;
      std::vector<int> cand = sel; cand.push_back(j);
      VectorXd g; double o = obj_of(cand, g); int np, nn;
      if (o < best - 1e-12 && dir_feasible(g, lambda_pos, lambda_neg, tol, np, nn)) { best = o; bj = j; }
    }
    if (bj >= 0) { sel.push_back(bj); obj = best; grow = true; }
  }
  auto bestpair = polish(sel, obj);

  // random restarts
  Function set_seed("set.seed"); set_seed(seed);
  Function sample_int("sample.int");
  for (int r = 1; r < n_restart; ++r) {
    int sz = std::max(1, cap);
    IntegerVector pick = sample_int(K, std::min(sz, K));
    std::vector<int> s2; for (int t = 0; t < pick.size(); ++t) s2.push_back(pick[t] - 1);
    VectorXd g; int np, nn; double o2 = obj_of(s2, g);
    // trim to feasibility by dropping smallest-|gamma| until feasible
    while (!dir_feasible(g, lambda_pos, lambda_neg, tol, np, nn) && !s2.empty()) {
      int drop = 0; double mn = 1e300;
      for (size_t i = 0; i < s2.size(); ++i) if (std::abs(g(i)) < mn) { mn = std::abs(g(i)); drop = (int)i; }
      s2.erase(s2.begin() + drop); o2 = obj_of(s2, g);
    }
    auto pr = polish(s2, o2);
    if (pr.first < bestpair.first - 1e-12) bestpair = pr;
  }

  std::vector<int> best_sel = bestpair.second;
  std::sort(best_sel.begin(), best_sel.end());
  VectorXd gamma, beta_fix; double fobj = fast_obj(c, best_sel, &gamma, &beta_fix);
  VectorXd beta = VectorXd::Zero(p_fix + K);
  beta.head(p_fix) = beta_fix;
  std::vector<int> pos, neg;
  for (size_t i = 0; i < best_sel.size(); ++i) {
    beta(p_fix + best_sel[i]) = gamma(i);
    if (gamma(i) > tol) pos.push_back(best_sel[i] + 1);
    else if (gamma(i) < -tol) neg.push_back(best_sel[i] + 1);
  }
  return List::create(_["obj"] = fobj, _["beta"] = beta,
                      _["selected_pos"] = wrap(pos), _["selected_neg"] = wrap(neg));
}

// [[Rcpp::export]]
List scs_solve_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                   const Eigen::Map<Eigen::VectorXd> Y,
                   int p_fix, int K, IntegerVector lambda, double mu,
                   int n_restart = 4, int seed = 1) {
  const int q = lambda.size();
  std::vector<int> lam(q);
  for (int l = 0; l < q; ++l) lam[l] = lambda[l];

  auto run_from = [&](std::vector<std::vector<char>> active) {
    double obj = inner_obj(X, Y, build_support(p_fix, K, q, active), mu);
    while (swap_pass(X, Y, p_fix, K, q, mu, active, obj)) {}
    return std::make_pair(obj, active);
  };

  // start 1: gradient ranking
  std::vector<std::vector<char>> active(q, std::vector<char>(K, 0));
  rank_warmstart(X, Y, p_fix, K, q, lam, mu, active);
  auto best = run_from(active);

  // random restarts (R's RNG so set.seed controls it)
  Function set_seed("set.seed"); set_seed(seed);
  Function sample_int("sample.int");
  for (int r = 1; r < n_restart; ++r) {
    std::vector<std::vector<char>> a2(q, std::vector<char>(K, 0));
    for (int l = 0; l < q; ++l) {
      IntegerVector pick = sample_int(K, lam[l]);
      for (int t = 0; t < pick.size(); ++t) a2[l][pick[t] - 1] = 1;
    }
    auto res = run_from(a2);
    if (res.first < best.first - 1e-12) best = res;
  }

  // assemble beta (length p = p_fix + q*K) and selected sets
  std::vector<int> ind = build_support(p_fix, K, q, best.second);
  VectorXd grad, coef;
  double obj = inner_obj(X, Y, ind, mu, &grad, &coef);
  VectorXd beta = VectorXd::Zero(p_fix + q * K);
  for (size_t j = 0; j < ind.size(); ++j) beta(ind[j]) = coef(j);

  List sel(q);
  for (int l = 0; l < q; ++l) {
    std::vector<int> s;
    for (int kk = 0; kk < K; ++kk) if (best.second[l][kk]) s.push_back(kk + 1); // 1-based
    sel[l] = wrap(s);
  }
  return List::create(_["obj"] = obj, _["beta"] = beta, _["selected"] = sel);
}
