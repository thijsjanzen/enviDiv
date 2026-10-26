// envidiv likelihood, version 4: hybrid of closed forms and Boost odeint.
//
// Closed forms are used for everything that follows constant-rate
// birth-death dynamics: E_H, E_S (logistic solution) and log D_H, log D_S
// along branches. The three quantities without a closed form are integrated
// with an adaptive odeint stepper (as in secsse / TRAISIE), each as a
// one-dimensional linear ODE in which E_S and log D_S come from the closed
// forms:
//   * E_U across each lowstand (once per lowstand):
//       E_U' = 2 mu E_S + 2 alpha E_S^2 + 2 lambda_L E_S^3 - 2 theta E_U
//   * rho = D_U / D_S along branch segments inside a lowstand:
//       rho' = h - rho (2 theta + g),  h = 2 mu + 4 alpha E_S + 6 lambda_L E_S^2,
//       g = d log D_S / ds = -(lambda_L + mu) + 2 lambda_L E_S
//   * the pending pair below a pump node, Qt = Q / (D_A D_B), from tau_min to
//     the drop age a:
//       Qt' = (2 alpha + 4 lambda_L E_S) - Qt (2 theta + 2 g),
//       Qt(tau_min) = 1 if the older first daughter node lies inside the
//       lowstand (point mass), 0 otherwise;
//     then log D_H(a) = log D_A(tau_min) + log D_B(tau_min)
//                       + 2 [log D_S(a) - log D_S(tau_min)] + log Qt(a).
// Ages are cut into segments at all rises and drops; per segment the start
// values of F = 1 - E, LD (cumulative log D increment) and E_U are stored.
// Includes the dating error sigma on pump nodes (mixture over node labels),
// the sampling fraction f, and consistent treatment of node ages within BTOL
// of a lowstand boundary. The time-scale factor is applied in R.
//

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include "boost/numeric/odeint.hpp"
#include <array>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <functional>
#include <tuple>
using namespace Rcpp;
namespace odeint = boost::numeric::odeint;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

inline double gfun(double r, double dt) {
  return std::fabs(r * dt) < 1e-12 ? dt : std::expm1(r * dt) / r;
}
inline double prop_F(double F0, double l, double m, double dt) {
  double r = l - m;
  return F0 * std::exp(r * dt) / (1.0 + l * F0 * gfun(r, dt));
}
inline double prop_logD(double F0, double l, double m, double dt) {
  double r = l - m;
  return r * dt - 2.0 * std::log1p(l * F0 * gfun(r, dt));
}

// ---- odeint for one-dimensional linear ODEs
using state1 = std::array<double, 1>;
enum class Method { RKCK54, RKF78, RKD5, BS };

struct Solver {
  Method method;
  double atol, rtol;
  Solver(const std::string& m = "odeint::runge_kutta_cash_karp54",
         double a = 1e-12, double r = 1e-10) : atol(a), rtol(r) {
    if (m == "odeint::runge_kutta_cash_karp54")     method = Method::RKCK54;
    else if (m == "odeint::runge_kutta_fehlberg78") method = Method::RKF78;
    else if (m == "odeint::runge_kutta_dopri5")     method = Method::RKD5;
    else if (m == "odeint::bulirsch_stoer")         method = Method::BS;
    else Rcpp::stop("unknown method: " + m);
  }
  // integrate x' = f(x, t) from t0 to t1; returns x(t1)
  template <class F>
  double run(F f, double x0, double t0, double t1) const {
    if (!(t1 > t0)) return x0;
    state1 x = {x0};
    auto sys = [&f](const state1& y, state1& dy, double t) { dy[0] = f(y[0], t); };
    double dt = 0.1 * (t1 - t0);
    switch (method) {
      case Method::RKCK54:
        odeint::integrate_adaptive(odeint::make_controlled<odeint::runge_kutta_cash_karp54<state1>>(atol, rtol),
                                   sys, x, t0, t1, dt); break;
      case Method::RKF78:
        odeint::integrate_adaptive(odeint::make_controlled<odeint::runge_kutta_fehlberg78<state1>>(atol, rtol),
                                   sys, x, t0, t1, dt); break;
      case Method::RKD5:
        odeint::integrate_adaptive(odeint::make_controlled<odeint::runge_kutta_dopri5<state1>>(atol, rtol),
                                   sys, x, t0, t1, dt); break;
      case Method::BS: {
        odeint::bulirsch_stoer<state1> st(atol, rtol);
        odeint::integrate_adaptive(st, sys, x, t0, t1, dt); break;
      }
    }
    return x[0];
  }
};

struct State {
  double logD = 0.0;
  double rho = 0.0;
  bool low = false;
  double s = 0.0;
};

struct Seg {
  double b, e;      // start and end age
  bool low;
  int stand;        // -1 in high water
  double lam;       // speciation rate in this segment
  double th2;       // 2 * theta (low water only)
  double F, LD, TH, EU;  // values at the start of the segment
};

struct Env {
  const double mu, lH, lL, al, f;
  Env(double mu_, double lH_, double lL_, double al_, double f_) :
    mu(mu_), lH(lH_), lL(lL_), al(al_), f(f_) {}

  std::vector<double> rise, drop;
  Solver solver;
  std::vector<Seg> seg;
  std::vector<int> first_seg, last_seg;   // per stand

  int stand_of(double s, double tol = 1e-9) const {
    for (size_t j = 0; j < rise.size(); ++j)
      if (s > rise[j] + tol && s < drop[j] - tol) return (int)j;
    return -1;
  }
  int seg_at(double s) const {         // largest segment with start <= s
    int lo = 0, hi = (int)seg.size() - 1;
    while (lo < hi) { int mid = (lo + hi + 1) / 2; if (seg[mid].b <= s) lo = mid; else hi = mid - 1; }
    return lo;
  }
  int seg_in_stand(int j, double s) const {
    int i = first_seg[j];
    while (i < last_seg[j] && seg[i + 1].b <= s) ++i;
    return i;
  }
  double F_in(int i,  double s) const { return prop_F(seg[i].F, seg[i].lam, mu, s - seg[i].b); }
  double LD_in(int i, double s) const { return seg[i].LD + prop_logD(seg[i].F, seg[i].lam, mu, s - seg[i].b); }
  double TH_in(int i, double s) const { return seg[i].TH + (seg[i].low ? seg[i].th2 * (s - seg[i].b) : 0.0); }

  double FH_at(double s)      const { return F_in(seg_at(s), s); }
  double LD_at(double s)      const { return LD_in(seg_at(s), s); }
  double FS(int j, double s)  const { return F_in(seg_in_stand(j, s), s); }
  double LDl(int j, double s) const { return LD_in(seg_in_stand(j, s), s); }
  double THl(int j, double s) const { return TH_in(seg_in_stand(j, s), s); }
  // speciation rate at a node: water state from stand_of() (with tolerance,
  // as everywhere else)
  double lam_node(double s) const { return stand_of(s) >= 0 ? lL : lH; }
};

Env make_env(double mu, double lH, double lL, double al, double f,
             const std::vector<double>& rise, const std::vector<double>& drop,
             const Solver& solver) {
  Env e(mu, lH, lL, al, f);
  e.solver = solver;
  e.rise = rise; e.drop = drop;
  std::vector<double> B = {0.0};
  for (double x : rise) B.push_back(x);
  for (double x : drop) B.push_back(x);
  std::sort(B.begin(), B.end());
  B.erase(std::unique(B.begin(), B.end(), [](double a, double b) { return std::fabs(a - b) < 1e-12; }), B.end());
  size_t n = B.size();
  e.seg.resize(n);
  for (size_t i = 0; i < n; ++i) {
    Seg& g = e.seg[i];
    g.b = B[i];
    g.e = (i + 1 < n) ? B[i + 1] : std::numeric_limits<double>::infinity();
    double mid = std::isfinite(g.e) ? 0.5 * (g.b + g.e) : g.b + 1.0;
    g.stand = e.stand_of(mid, 0.0);
    g.low = g.stand >= 0;
    g.lam = g.low ? lL : lH;
    g.th2 = g.low ? 2.0 * (al + lL + mu) : 0.0;
  }
  e.first_seg.assign(rise.size(), -1); e.last_seg.assign(rise.size(), -1);
  for (size_t i = 0; i < n; ++i) if (e.seg[i].low) {
    int j = e.seg[i].stand;
    if (e.first_seg[j] < 0) e.first_seg[j] = (int)i;
    e.last_seg[j] = (int)i;
  }
  // pass over segments upwards in age
  e.seg[0].F = f; e.seg[0].LD = 0.0; e.seg[0].TH = 0.0; e.seg[0].EU = 1.0 - f;
  for (size_t i = 0; i + 1 < n; ++i) {
    Seg& g = e.seg[i];
    double D = g.e - g.b;
    double Fend = prop_F(g.F, g.lam, mu, D);
    double LDend = g.LD + prop_logD(g.F, g.lam, mu, D);
    double THend = g.TH + (g.low ? g.th2 * D : 0.0);
    double EUend = 0.0;
    if (g.low) {
      double F0 = g.F, lam = g.lam, th2 = g.th2, b = g.b, en = g.e;
      auto dEU = [&](double EU, double u) {
        double ES = 1.0 - prop_F(F0, lam, mu, u - b);
        return 2 * mu * ES + 2 * al * ES * ES + 2 * lam * ES * ES * ES - th2 * EU;
      };
      EUend = solver.run(dEU, g.EU, b, en);
    }
    Seg& h = e.seg[i + 1];
    h.LD = LDend; h.TH = THend;
    if (!g.low && h.low) {            // rise: S and U start from H
      h.F = Fend; h.EU = 1.0 - Fend;
    } else if (g.low && !h.low) {     // drop: H starts from U
      h.F = 1.0 - EUend; h.EU = 0.0;
    } else {                          // (adjacent segments in the same state)
      h.F = Fend; h.EU = g.low ? EUend : 0.0;
    }
  }
  return e;
}

// D_U / D_S propagation within stand j from s0 to s1 (odeint)
State prop_low(const Env& e, int j, double s0, double s1, State st) {
  if (s1 <= s0) return st;
  double LD1 = e.LDl(j, s1), LD0 = e.LDl(j, s0);
  double rho = st.rho;
  int i0 = e.seg_in_stand(j, s0), i1 = e.seg_in_stand(j, s1);
  for (int i = i0; i <= i1; ++i) {
    double a = std::max(s0, e.seg[i].b), b = std::min(s1, e.seg[i].e);
    if (b <= a) continue;
    double lam = e.seg[i].lam, th2 = e.seg[i].th2;
    auto drho = [&](double r, double u) {
      double ES = 1.0 - e.F_in(i, u);
      double h = 2 * e.mu + 4 * e.al * ES + 6 * lam * ES * ES;
      double g = -(lam + e.mu) + 2 * lam * ES;
      return h - r * (th2 + g);
    };
    rho = e.solver.run(drho, rho, a, b);
  }
  st.rho = rho;
  st.logD += LD1 - LD0;
  st.low = true;
  st.s = s1;
  return st;
}

// Propagate a branch state from age s_lo to s_hi. Ages within BTOL of a
// lowstand boundary are treated as lying on the boundary, consistent with
// stand_of(): a parent within BTOL of a drop or rise is in high water, so the
// state is converted to high water there; a branch starting within BTOL above
// a rise still passes through that lowstand.
const double BTOL = 1e-9;
State prop_branch(const Env& e, double s_lo, double s_hi, State st, int stop_j = -1) {
  double s = s_lo;
  int j = e.stand_of(s);
  while (true) {
    if (j < 0) {
      double nxt = s_hi;
      int jn = -1;
      for (size_t k = 0; k < e.rise.size(); ++k) {
        double r = e.rise[k];
        if (r >= s - BTOL && r < s_hi - BTOL && r < nxt) { nxt = r; jn = (int)k; }
      }
      if (jn < 0 || nxt > s) st.logD += e.LD_at(std::max(nxt, s)) - e.LD_at(s);
      if (jn < 0) { st.s = s_hi; return st; }
      s = nxt; j = jn;
      st.rho = 1.0; st.low = true;
      if (stop_j == j) { st.s = s; return st; }
    } else {
      if (stop_j == j) { st.s = s; return st; }
      double top = std::min(e.drop[j], s_hi);
      st = prop_low(e, j, s, top, st);
      s = top;
      if (s_hi < e.drop[j] - BTOL) { st.s = s; return st; }   // parent inside the lowstand
      st.logD += std::log(st.rho);                             // reached the drop: U -> H
      st.low = false;
      j = -1;
      s = top;
      if (s >= s_hi - 1e-12) { st.s = s_hi; return st; }
    }
  }
}

double lse(const std::vector<double>& x) {
  double m = NEG_INF;
  for (double v : x) m = std::max(m, v);
  if (!std::isfinite(m)) return m;
  double s = 0;
  for (double v : x) s += std::exp(v - m);
  return m + std::log(s);
}

bool sum_states(const std::vector<State>& v, State& out) {
  if (v.empty()) return false;
  std::vector<double> l; for (auto& s : v) l.push_back(s.logD);
  double tot = lse(l);
  out = v[0]; out.logD = tot;
  if (v[0].low) {
    double r = 0; for (auto& s : v) r += std::exp(s.logD - tot) * s.rho;
    out.rho = r;
  }
  return true;
}

struct Hyp { double age; int j; State st; };

}  // namespace

// [[Rcpp::export]]
double envidiv_loglik_cpp(NumericVector pars,
                          IntegerMatrix edge,
                          NumericVector ages,
                          int n_tip,
                          NumericVector rise,
                          NumericVector drop,
                          double sigma = 0.0,
                          double max_z = 5.0,
                          double tol = 1e-6,
                          double f = 1.0,
                          std::string method = "odeint::runge_kutta_cash_karp54",
                          double atol = 1e-12,
                          double rtol = 1e-10) {
  Solver solver(method, atol, rtol);
  Env e = make_env(pars["mu"], pars["lambda_H"], pars["lambda_L"], pars["alpha"],
                   f,
                   as<std::vector<double>>(rise),
                   as<std::vector<double>>(drop),
                   solver);
  int nn = ages.size();
  std::vector<double> age = as<std::vector<double>>(ages);
  std::vector<std::vector<int>> kids(nn);
  std::vector<int> parent(nn, -1);
  for (int r = 0; r < edge.nrow(); ++r) {
    kids[edge(r, 0) - 1].push_back(edge(r, 1) - 1);
    parent[edge(r, 1) - 1] = edge(r, 0) - 1;
  }
  int crown = n_tip;
  double T = age[crown];
  if (e.stand_of(T) >= 0) stop("crown inside a lowstand");
  const std::vector<double>& D = e.drop;
  const double logf = std::log(f);

  std::function<std::vector<Hyp>(int)> node_hyps = [&](int v) -> std::vector<Hyp> {
    std::vector<Hyp> out;
    if (v < n_tip) { State t; t.logD = logf; out.push_back({0.0, -1, t}); return out; }
    double x = age[v];
    std::vector<std::tuple<double, int, double>> H;
    if (sigma <= 0) {
      int jj = -1;
      for (size_t j = 0; j < D.size(); ++j) if (std::fabs(D[j] - x) < tol) { jj = (int)j; break; }
      H.emplace_back(jj >= 0 ? D[jj] : x, jj, 0.0);
    } else {
      H.emplace_back(x, -1, 0.0);
      double kid_max = 0;
      for (int k : kids[v]) kid_max = std::max(kid_max, age[k]);
      for (size_t j = 0; j < D.size(); ++j) {
        double a = D[j];
        if (std::fabs(x - a) >= max_z * sigma) continue;
        double tr = std::min(age[parent[v]] - a, a - kid_max);
        if (tr <= 0 || std::fabs(x - a) > tr) continue;
        double lw = R::dnorm(x, a, sigma, 1) -
          std::log(R::pnorm(tr / sigma, 0, 1, 1, 0) - R::pnorm(-tr / sigma, 0, 1, 1, 0));
        H.emplace_back(a, (int)j, lw);
      }
    }
    std::vector<std::vector<Hyp>> kh;
    for (int k : kids[v]) kh.push_back(node_hyps(k));
    for (auto& h : H) {
      double sv = std::get<0>(h); int jp = std::get<1>(h); double lw = std::get<2>(h);
      State st;
      if (jp < 0) {
        std::vector<State> ds(2); bool ok = true;
        for (int c = 0; c < 2; ++c) {
          std::vector<State> v2;
          for (auto& hc : kh[c]) if (hc.age < sv - 1e-12) v2.push_back(prop_branch(e, hc.age, sv, hc.st));
          if (!sum_states(v2, ds[c])) ok = false;
        }
        if (!ok) continue;
        int js = e.stand_of(sv);
        st.logD = std::log(e.lam_node(sv)) + ds[0].logD + ds[1].logD;
        if (js < 0) {
          st.low = false;
        } else {
          st.rho = 2.0 * (1.0 - e.FS(js, sv)); st.low = true;
        }
      } else {
        double a1 = e.rise[jp], a2 = e.drop[jp];
        std::vector<std::vector<State>> starts(2);
        for (int c = 0; c < 2; ++c)
          for (auto& hc : kh[c]) if (hc.age < a2 - 1e-12)
            starts[c].push_back(prop_branch(e, hc.age, a2, hc.st, jp));
        if (starts[0].empty() || starts[1].empty()) continue;
        std::vector<double> terms;
        for (const State& A0 : starts[0]) for (const State& B0 : starts[1]) {
          double tmin = std::max(A0.s, B0.s);
          State A = (A0.s < tmin) ? prop_low(e, jp, A0.s, tmin, A0) : A0;
          State B = (B0.s < tmin) ? prop_low(e, jp, B0.s, tmin, B0) : B0;
          double Qt = (tmin > a1 + 1e-12) ? 1.0 : 0.0;          // point mass
          int i0 = e.seg_in_stand(jp, tmin), i1 = e.seg_in_stand(jp, a2);
          for (int i = i0; i <= i1; ++i) {
            double lo = std::max(tmin, e.seg[i].b), hi = std::min(a2, e.seg[i].e);
            if (hi <= lo) continue;
            double lam = e.seg[i].lam, th2 = e.seg[i].th2;
            auto dQ = [&](double q, double tau) {
              double ES = 1.0 - e.F_in(i, tau);
              double g = -(lam + e.mu) + 2 * lam * ES;
              return (2 * e.al + 4 * lam * ES) - q * (th2 + 2 * g);
            };
            Qt = e.solver.run(dQ, Qt, lo, hi);
          }
          double logQt = (Qt > 0) ? std::log(Qt) : NEG_INF;
          double G = e.LDl(jp, a2) - e.LDl(jp, tmin);
          terms.push_back(A.logD + B.logD + 2.0 * G + logQt);
        }
        st.logD = lse(terms); st.low = false;
      }
      st.logD += lw;
      out.push_back({sv, jp, st});
    }
    return out;
  };

  double ll = 0;
  for (int k : kids[crown]) {
    std::vector<State> v2;
    for (auto& h : node_hyps(k)) if (h.age < T - 1e-12) v2.push_back(prop_branch(e, h.age, T, h.st));
    State s;
    if (!sum_states(v2, s)) return NEG_INF;
    ll += s.logD;
  }
  return ll - 2.0 * std::log(e.FH_at(T));
}
