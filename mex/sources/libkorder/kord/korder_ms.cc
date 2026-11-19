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

#include "korder_ms.hh"
#include "kord_exception.hh"

/* A + λ [0 0 f_{y**₊}], built as MatrixS with the last summand scaled by λ. */
MatrixSL::MatrixSL(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
                   const PartitionY& ypart, double lambda) :
    PLUMatrix(ypart.ny())
{
  zeros();
  IntSequence c {1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  ConstTwoDMatrix gss_ys(ypart.nstat + ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);
  TwoDMatrix sub2(*this, ypart.nstat + ypart.npred, ypart.nyss());
  sub2.add(lambda, f_yss);
  calcPLU();
}

KOrderMS::KOrderMS(int num_stat, int num_pred, int num_both, int num_forw,
                   const TensorContainer<FSSparseTensor>& fcont, const TwoDMatrix& gy,
                   const TwoDMatrix& gu, const TwoDMatrix& v, const UNormalMoments& moments,
                   const TwoDMatrix& P_, const Vector& lam_, const TwoDMatrix& W_,
                   const TwoDMatrix& Winv_, const TwoDMatrix& forcing_, Journal& jr) :
    S(P_.nrows()),
    ypart(num_stat, num_pred, num_both, num_forw),
    ny(ypart.ny()),
    nu(gu.ncols()),
    P(P_),
    lam(lam_),
    W(W_),
    Winv(Winv_),
    forcing(forcing_),
    journal(jr)
{
  KORD_RAISE_IF(P.ncols() != S, "Transition matrix is not square in KOrderMS");
  KORD_RAISE_IF(lam.length() != S || W.nrows() != S || W.ncols() != S || Winv.nrows() != S
                    || Winv.ncols() != S,
                "Wrong dimensions of the eigen-decomposition in KOrderMS");
  KORD_RAISE_IF(forcing.nrows() != ny || forcing.ncols() != S,
                "Forcing must be ny × S in KOrderMS");
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      {
        auto p = std::make_unique<KOrder>(num_stat, num_pred, num_both, num_forw, fcont, gy, gu, v,
                                          moments, jr);
        p->switchToFolded();
        p->_fGstack.inner_sigma_zero = false;
        p->_uGstack.inner_sigma_zero = false;
        pairs.push_back(std::move(p));
      }
  for (int m = 0; m < S; m++)
    matSL.push_back(std::make_unique<MatrixSL>(fcont.get(Symmetry {1}),
                                               pairs[0]->_uZstack.getStackSizes(), gy, ypart,
                                               lam[m]));
}

std::vector<std::unique_ptr<FGSTensor>>
KOrderMS::toEigen(const std::vector<std::unique_ptr<FGSTensor>>& X) const
{
  std::vector<std::unique_ptr<FGSTensor>> Xt;
  for (int m = 0; m < S; m++)
    {
      auto t = std::make_unique<FGSTensor>(*X[0]);
      t->zeros();
      for (int s = 0; s < S; s++)
        t->add(W.get(s, m), *X[s]);
      Xt.push_back(std::move(t));
    }
  return Xt;
}

std::vector<std::unique_ptr<FGSTensor>>
KOrderMS::fromEigen(const std::vector<std::unique_ptr<FGSTensor>>& Xt) const
{
  std::vector<std::unique_ptr<FGSTensor>> X;
  for (int s = 0; s < S; s++)
    {
      auto t = std::make_unique<FGSTensor>(*Xt[0]);
      t->zeros();
      for (int m = 0; m < S; m++)
        t->add(Winv.get(m, s), *Xt[m]);
      X.push_back(std::move(t));
    }
  return X;
}

/* Insert X(s) as the derivative of regime s: an owned copy in g of pairs (s, ·); g* of pairs (s, ·) and g** of
   pairs (·, s) are row views (FGSTensor(first_row, num, t) shares the storage of t, as in KOrder::insertDerivative)
   into the g tensor of pair (s, s), which lives as long as the solver. */
void
KOrderMS::insertMS(std::vector<std::unique_ptr<FGSTensor>>& X)
{
  const Symmetry& sym = X[0]->getSym();
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      pr(s, sp).g<Storage::fold>().insert(std::make_unique<FGSTensor>(*X[s]));
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      {
        KOrder& p = pr(s, sp);
        FGSTensor& own = p.g<Storage::fold>().get(sym);
        FGSTensor& outer = pr(sp, sp).g<Storage::fold>().get(sym);
        p.gs<Storage::fold>().insert(std::make_unique<FGSTensor>(ypart.nstat, ypart.nys(), own));
        p.gss<Storage::fold>().insert(
            std::make_unique<FGSTensor>(ypart.nstat + ypart.npred, ypart.nyss(), outer));
      }
}

/* Complete the conditional G of each pair: the term with the new inner tensor (outer dimension 1) and, if
   dim_outer > 0, the term with the new outer tensor X(s′) only. */
void
KOrderMS::updateG(const std::vector<FGSTensor*>& Gc,
                  const std::vector<std::unique_ptr<FGSTensor>>& X, int dim_outer)
{
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      {
        KOrder& p = pr(s, sp);
        FGSTensor& G = *Gc[s * S + sp];
        p.Gstack<Storage::fold>().multAndAdd(1, p.gss<Storage::fold>(), G);
        if (dim_outer > 0)
          {
            FGSContainer tmp(4);
            tmp.insert(std::make_unique<FGSTensor>(ypart.nstat + ypart.npred, ypart.nyss(),
                                                   pr(sp, sp).g<Storage::fold>().get(X[sp]->getSym())));
            p.Gstack<Storage::fold>().multAndAdd(dim_outer, tmp, G);
          }
      }
}

/* For symmetry yⁱuʲσᵏ (k ≥ 1): in every pair, insert G_yⁱuʲu′ᵐσᵏ⁻ᵐ (m = 1…k, no moment gate), the conditional
   G_yⁱuʲσᵏ, and compute F = F_yⁱuʲσᵏ + Dᵢⱼₖ + Eᵢⱼₖ. Returns the right-hand sides −Σ_{s′} P[s,s′] F(s, s′) in Fs and
   pointers to the conditional G of each pair. */
std::vector<FGSTensor*>
KOrderMS::rhs(int i, int j, int k, std::vector<std::unique_ptr<FGSTensor>>& Fs)
{
  Symmetry sym {i, j, 0, k};
  JournalRecordPair pa(journal);
  pa << "MS: recovering symmetry " << sym << endrec;
  std::vector<FGSTensor*> Gc(S * S);
  std::vector<std::unique_ptr<FGSTensor>> Fp(S * S);
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      {
        KOrder& p = pr(s, sp);
        for (int m = 1; m <= k; m++)
          p.G<Storage::fold>().insert(p.faaDiBrunoG<Storage::fold>(Symmetry {i, j, m, k - m}));
        auto G = p.faaDiBrunoG<Storage::fold>(sym);
        Gc[s * S + sp] = G.get();
        p.G<Storage::fold>().insert(std::move(G));
        auto F = p.faaDiBrunoZ<Storage::fold>(sym);
        F->add(1.0, p.calcD_ijk<Storage::fold>(i, j, k));
        F->add(1.0, p.calcE_ijk<Storage::fold>(i, j, k));
        Fp[s * S + sp] = std::move(F);
      }
  Fs.clear();
  for (int s = 0; s < S; s++)
    {
      auto t = std::make_unique<FGSTensor>(*Fp[s * S]);
      t->zeros();
      for (int sp = 0; sp < S; sp++)
        t->add(P.get(s, sp), *Fp[s * S + sp]);
      t->mult(-1.0);
      Fs.push_back(std::move(t));
    }
  return Gc;
}

/* Order 1: g_y, g_u are those of the constant model (already in every pair). g_σ(s) solves
   A c(s) + [0 0 f_{y**₊}] Σ_{s′} P[s,s′] c(s′) = −forcing(s); then G_σ per pair. */
void
KOrderMS::solveOrder1()
{
  JournalRecordPair pa(journal);
  pa << "MS: first order sigma derivatives" << endrec;
  std::vector<std::unique_ptr<FGSTensor>> X;
  for (int s = 0; s < S; s++)
    {
      auto t = std::make_unique<FGSTensor>(ny,
                                           TensorDimens(Symmetry {0, 0, 0, 1}, pairs[0]->nvs));
      for (int r = 0; r < ny; r++)
        t->get(r, 0) = -forcing.get(r, s);
      X.push_back(std::move(t));
    }
  auto Xt = toEigen(X);
  for (int m = 0; m < S; m++)
    matSL[m]->multInv(*Xt[m]);
  X = fromEigen(Xt);
  insertMS(X);
  for (auto& p : pairs)
    p->G<Storage::fold>().insert(p->faaDiBrunoG<Storage::fold>(Symmetry {0, 0, 0, 1}));
}

void
KOrderMS::recover_ys(int i, int j)
{
  std::vector<std::unique_ptr<FGSTensor>> R;
  auto Gc = rhs(i, 0, j, R);
  auto Rt = toEigen(R);
  KOrder& p0 = *pairs[0];
  for (int m = 0; m < S; m++)
    {
      FGSTensor& der = *Rt[m];
      if (ypart.nys() > 0 && ypart.nyss() > 0)
        {
          UGSTensor tmp(der);
          TwoDMatrix gs_y(p0.gs<Storage::unfold>().get(Symmetry {1, 0, 0, 0}));
          TwoDMatrix Bl(p0.matB);
          Bl.mult(lam[m]);
          GeneralSylvester sylv(i, ny, ypart.nys(), ypart.nstat + ypart.npred, p0.matA.getData(),
                                Bl.getData(), gs_y.getData(), tmp.getData());
          sylv.solve();
          FGSTensor ftmp(tmp);
          der.getData() = const_cast<const Vector&>(ftmp.getData());
        }
      else if (ypart.nys() > 0 && ypart.nyss() == 0)
        p0.matA.multInv(der);
    }
  auto X = fromEigen(Rt);
  insertMS(X);
  updateG(Gc, X, i + j);
}

void
KOrderMS::recover_yus(int i, int j, int k)
{
  std::vector<std::unique_ptr<FGSTensor>> R;
  auto Gc = rhs(i, j, k, R);
  for (int s = 0; s < S; s++)
    pairs[0]->matA.multInv(*R[s]);
  insertMS(R);
  updateG(Gc, R, 0);
}

void
KOrderMS::recover_s(int i)
{
  std::vector<std::unique_ptr<FGSTensor>> R;
  auto Gc = rhs(0, 0, i, R);
  auto Rt = toEigen(R);
  for (int m = 0; m < S; m++)
    matSL[m]->multInv(*Rt[m]);
  auto X = fromEigen(Rt);
  insertMS(X);
  updateG(Gc, X, i);
}

/* The order of KOrder::performStep(); σ-free symmetries by the stock routines in every pair. */
void
KOrderMS::performStep(int order)
{
  JournalRecordPair pa(journal);
  pa << "MS: performing step for order = " << order << endrec;
  for (auto& p : pairs)
    p->recover_y<Storage::fold>(order);
  for (int i = 0; i < order; i++)
    for (auto& p : pairs)
      p->recover_yu<Storage::fold>(i, order - i);
  for (int j = 1; j < order; j++)
    {
      for (int i = j - 1; i >= 1; i--)
        recover_yus(order - j, i, j - i);
      recover_ys(order - j, j);
    }
  for (int i = order - 1; i >= 1; i--)
    recover_yus(0, i, order - i);
  recover_s(order);
}

void
KOrderMS::solve(int maxorder)
{
  solveOrder1();
  for (int k = 2; k <= maxorder; k++)
    performStep(k);
}
