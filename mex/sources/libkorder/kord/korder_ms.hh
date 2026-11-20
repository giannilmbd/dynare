/*
 * Copyright © 2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Perturbation of order k with an exogenous Markov chain s_t ∈ {1,…,S} with transition matrix P.

   Switching parameters are endogenous variables a with a dummy equation whose residual is a_t − σ â(s_t):
   the regime enters the system only through the σ¹ forcing −â(s) in those equation rows (‘forcing’).
   The policy is regime dependent, y_t = g(y*_{t−1}, u_t, σ; s_t). Next period the composite is
   G(·; s, s′) = g**(g*(·; s), u′, σ; s′), so each regime pair (s, s′) gets its own KOrder object: in pair (s, s′)
   the containers g and g* hold regime s (current and inner function), g** holds regime s′ (outer function), and G
   and the Z stack are those of the pair. The equations of regime s are Σ_{s′} P[s,s′] F(s, s′) = 0.

   Derivatives without σ are regime free and are solved by the stock routines in every pair. Derivatives with σ:
     y^i σ^j:   A X(s) + B Σ_{s′}P[s,s′] X(s′) (g*_y)^{⊗i} = R(s)  (coupled Sylvester),
     y^i u^j σ^k (j ≥ 1): A X(s) = R(s)                             (per regime),
     σ^k:       A X(s) + [0 0 f_{y**₊}] Σ_{s′}P[s,s′] X(s′) = R(s)   (coupled).
   The coupled systems are solved jointly for the S regimes: the σ^k system as one LU of order S·ny,
   (I_S ⊗ A₀ + P ⊗ C) x = r, and the y^iσ^j system as one generalised Sylvester equation of order S·ny,
   Ā x̃ + B̄ x̃ (g*_y)^{⊗i} = r, with B̄ = P ⊗ B and the unknowns ordered so that the columns of B̄ act on the last
   S·nyss rows (MatrixSMS, MatrixAMS, MatrixBMS). No decomposition of P is needed, so P may have complex eigenvalues
   or be defective. See MS_korder_algebra.md for the derivation.

   Differences from the stock KOrder step: g_σ(s) ≠ 0 at first order (inserted, with G_σ); the σ-order gates
   hasMoment(k) and the fillG gate hasMoment(k−m) are removed (odd σ derivatives are non-zero with regimes); after
   solving y^iσ^j or σ^k, G is updated with the new outer tensor only (the stock update with the whole g**
   container would count twice the terms g**_{y^{i+m}σ^{j−m}}(g*_σ)^m that the conditional G already holds). */

#ifndef KORDER_MS_H
#define KORDER_MS_H

#include "korder.hh"

#include <memory>
#include <vector>

/* I_S ⊗ A₀ + P ⊗ C, A₀ = [f_y] + [0 [f_y**₊]·[g**_y*] 0], C = [0 0 f_{y**₊}]: the stacked matrix of the coupled σ^k
   equations; row and column s·ny + r is variable r of regime s. With S = 1, P = 1 it is MatrixS. */
class MatrixSMS : public PLUMatrix
{
public:
  MatrixSMS(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
            const PartitionY& ypart, const TwoDMatrix& P);
};

class KOrderMS
{
  const int S;
  const PartitionY ypart;
  const int ny, nu;
  TwoDMatrix P;       // S × S, P(s, s′)
  TwoDMatrix forcing; // ny × S, added to F_σ(s) (equation rows)
  /* One journal per pair, not written on disk: the pairs of a σ symmetry are computed in parallel
     (thread_pool jobs) and Journal is not thread safe. Declared before ‘pairs’, which refer to them. */
  std::vector<std::unique_ptr<Journal>> pair_journals;
  std::vector<std::unique_ptr<KOrder>> pairs;
  std::unique_ptr<MatrixSMS> matS; // stacked σ^k matrix
  /* Stacked y^iσ^j Sylvester: Ā (S·ny × S·ny) and B̄ = P ⊗ B (S·ny × S·nyss). Unknown x̃: rows of regime s
     0…nstat+npred−1 at s·(nstat+npred) + r, rows nstat+npred… at S·(nstat+npred) + s·nyss + r − nstat − npred. */
  TwoDMatrix Abar, Bbar;
  Journal& journal;

  KOrder&
  pr(int s, int sp)
  {
    return *pairs[s * S + sp];
  }
  /* Pair (s, s′) enters F(s) = Σ_{s′} P[s,s′] F(s,s′) only if P[s,s′] ≠ 0; the G of a pair feeds only its
     own F. Diagonal pairs hold the derivatives of their regime (g** views of other pairs point into them)
     and are always computed. */
  [[nodiscard]] bool
  active(int s, int sp) const
  {
    return s == sp || P.get(s, sp) != 0.0;
  }

  // σ-free symmetry: g and G do not depend on the regime; copy those of pair (0, 0) into the other pairs
  void shareSigmaFree(const Symmetry& sym);
  void solveOrder1();
  void performStep(int order);
  // Common part for a symmetry with σ: fills G, computes F(s) = Σ P F(s,s′); returns the conditional G per pair
  std::vector<FGSTensor*> rhs(int i, int j, int k, std::vector<std::unique_ptr<FGSTensor>>& Fs);
  void insertMS(std::vector<std::unique_ptr<FGSTensor>>& X);
  void updateG(const std::vector<FGSTensor*>& Gc, const std::vector<std::unique_ptr<FGSTensor>>& X,
               int dim_outer);
  void recover_ys(int i, int j);
  void recover_yus(int i, int j, int k);
  void recover_s(int i);
  // X(s) ↦ x (S·ny rows, regime s in rows s·ny…), solve the stacked σ^k system, and back
  void solveStackedS(std::vector<std::unique_ptr<FGSTensor>>& X) const;

public:
  KOrderMS(int num_stat, int num_pred, int num_both, int num_forw,
           const TensorContainer<FSSparseTensor>& fcont, const TwoDMatrix& gy, const TwoDMatrix& gu,
           const TwoDMatrix& v, const UNormalMoments& moments, const TwoDMatrix& P_,
           const TwoDMatrix& forcing_, Journal& jr);
  void solve(int maxorder);
  [[nodiscard]] const FGSContainer&
  getFoldDers(int s) const
  {
    return pairs[s * S]->getFoldDers();
  }
  [[nodiscard]] int
  nRegimes() const
  {
    return S;
  }
};

#endif
