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
     σ^k:       A X(s) + [0 0 f_{y**₊}] Σ_{s′}P[s,s′] X(s′) = R(s)   (coupled),
   decoupled with the eigen-decomposition P′ W = W Λ (real eigenvalues; supplied by the caller):
   X̃_m = Σ_s X(s) W[s,m], X(s) = Σ_m X̃_m W⁻¹[m,s]. See MS_korder_algebra.md for the derivation.

   Differences from the stock KOrder step: g_σ(s) ≠ 0 at first order (inserted, with G_σ); the σ-order gates
   hasMoment(k) and the fillG gate hasMoment(k−m) are removed (odd σ derivatives are non-zero with regimes); after
   solving y^iσ^j or σ^k, G is updated with the new outer tensor only (the stock update with the whole g**
   container would count twice the terms g**_{y^{i+m}σ^{j−m}}(g*_σ)^m that the conditional G already holds). */

#ifndef KORDER_MS_H
#define KORDER_MS_H

#include "korder.hh"

#include <memory>
#include <vector>

/* A + λ [0 0 f_{y**₊}]: the matrix of the coupled σ^k equations for eigenvalue λ (λ = 1 is MatrixS). */
class MatrixSL : public PLUMatrix
{
public:
  MatrixSL(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
           const PartitionY& ypart, double lambda);
};

class KOrderMS
{
  const int S;
  const PartitionY ypart;
  const int ny, nu;
  TwoDMatrix P;     // S × S, P(s, s′)
  Vector lam;       // eigenvalues of P′
  TwoDMatrix W;     // P′ W = W diag(lam)
  TwoDMatrix Winv;  // W⁻¹
  TwoDMatrix forcing; // ny × S, added to F_σ(s) (equation rows)
  std::vector<std::unique_ptr<KOrder>> pairs;
  std::vector<std::unique_ptr<MatrixSL>> matSL;
  Journal& journal;

  KOrder&
  pr(int s, int sp)
  {
    return *pairs[s * S + sp];
  }

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
  // X̃_m = Σ_s X(s) W[s,m] and back
  std::vector<std::unique_ptr<FGSTensor>> toEigen(const std::vector<std::unique_ptr<FGSTensor>>& X) const;
  std::vector<std::unique_ptr<FGSTensor>> fromEigen(const std::vector<std::unique_ptr<FGSTensor>>& Xt) const;

public:
  KOrderMS(int num_stat, int num_pred, int num_both, int num_forw,
           const TensorContainer<FSSparseTensor>& fcont, const TwoDMatrix& gy, const TwoDMatrix& gu,
           const TwoDMatrix& v, const UNormalMoments& moments, const TwoDMatrix& P_,
           const Vector& lam_, const TwoDMatrix& W_, const TwoDMatrix& Winv_,
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
