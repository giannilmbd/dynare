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

#include "GeneralSylvester.hh"
#include "thread_pool.hh"

#include <functional>
#include <mutex>

/* I_S ⊗ A₀ + P ⊗ C. A₀ and C are assembled as in MatrixS (MatrixS = A₀ + C). */
MatrixSMS::MatrixSMS(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
                     const PartitionY& ypart, const TwoDMatrix& P) :
    PLUMatrix(P.nrows() * ypart.ny())
{
  const int ny = ypart.ny(), S = P.nrows();
  TwoDMatrix A0(ny, ny), C(ny, ny);
  A0.zeros();
  C.zeros();
  IntSequence c {1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  A0.add(1.0, f_y);
  ConstTwoDMatrix gss_ys(ypart.nstat + ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(A0, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);
  TwoDMatrix sub2(C, ypart.nstat + ypart.npred, ypart.nyss());
  sub2.add(1.0, f_yss);
  zeros();
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      for (int col = 0; col < ny; col++)
        for (int r = 0; r < ny; r++)
          get(s * ny + r, sp * ny + col)
              = (s == sp ? A0.get(r, col) : 0.0) + P.get(s, sp) * C.get(r, col);
  calcPLU();
}

KOrderMS::KOrderMS(int num_stat, int num_pred, int num_both, int num_forw,
                   const TensorContainer<FSSparseTensor>& fcont, const TwoDMatrix& gy,
                   const TwoDMatrix& gu, const TwoDMatrix& v, const UNormalMoments& moments,
                   const TwoDMatrix& P_, const TwoDMatrix& forcing_, Journal& jr) :
    S(P_.nrows()),
    ypart(num_stat, num_pred, num_both, num_forw),
    ny(ypart.ny()),
    nu(gu.ncols()),
    P(P_),
    forcing(forcing_),
    Abar(S * ypart.ny(), S * ypart.ny()),
    Bbar(S * ypart.ny(), S * ypart.nyss()),
    journal(jr)
{
  KORD_RAISE_IF(P.ncols() != S, "Transition matrix is not square in KOrderMS");
  KORD_RAISE_IF(forcing.nrows() != ny || forcing.ncols() != S,
                "Forcing must be ny × S in KOrderMS");
  /* Pair (0, 0) writes to the solver journal (it is the only pair whose records are kept: in a parallel
     group it is the only job that writes there, and the calling thread waits); the others to their own. */
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      {
        Journal* pj = &jr;
        if (s + sp > 0)
          {
            pair_journals.push_back(std::make_unique<Journal>());
            pj = pair_journals.back().get();
          }
        auto p = std::make_unique<KOrder>(num_stat, num_pred, num_both, num_forw, fcont, gy, gu, v,
                                          moments, *pj);
        p->switchToFolded();
        p->_fGstack.inner_sigma_zero = false;
        p->_uGstack.inner_sigma_zero = false;
        pairs.push_back(std::move(p));
      }
  matS = std::make_unique<MatrixSMS>(fcont.get(Symmetry {1}), pairs[0]->_uZstack.getStackSizes(), gy,
                                     ypart, P);
  const int zc = ypart.nstat + ypart.npred, nyss = ypart.nyss();
  const TwoDMatrix& A = pairs[0]->matA;
  const TwoDMatrix& B = pairs[0]->matB;
  Abar.zeros();
  Bbar.zeros();
  for (int s = 0; s < S; s++)
    for (int r = 0; r < ny; r++)
      {
        for (int col = 0; col < zc; col++)
          Abar.get(s * ny + r, s * zc + col) = A.get(r, col);
        for (int col = 0; col < nyss; col++)
          Abar.get(s * ny + r, S * zc + s * nyss + col) = A.get(r, zc + col);
        for (int sp = 0; sp < S; sp++)
          for (int col = 0; col < nyss; col++)
            Bbar.get(s * ny + r, sp * nyss + col) = P.get(s, sp) * B.get(r, col);
      }
}

void
KOrderMS::solveStackedS(std::vector<std::unique_ptr<FGSTensor>>& X) const
{
  const int nc = X[0]->ncols();
  TwoDMatrix x(S * ny, nc);
  for (int s = 0; s < S; s++)
    for (int j = 0; j < nc; j++)
      for (int r = 0; r < ny; r++)
        x.get(s * ny + r, j) = X[s]->get(r, j);
  matS->multInv(x);
  for (int s = 0; s < S; s++)
    for (int j = 0; j < nc; j++)
      for (int r = 0; r < ny; r++)
        X[s]->get(r, j) = x.get(s * ny + r, j);
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
        if (!active(s, sp))
          continue;
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

namespace
{
// A thread_pool job running a closure
class ClosureJob : public thread_pool::job
{
  std::function<void()> f;

public:
  explicit ClosureJob(std::function<void()> f_) : f(std::move(f_))
  {
  }
  void
  operator()([[maybe_unused]] std::mutex& mut) override
  {
    f();
  }
};
}

/* For symmetry yⁱuʲσᵏ (k ≥ 1): in every active pair, insert G_yⁱuʲu′ᵐσᵏ⁻ᵐ (m = 1…k, no moment gate), the
   conditional G_yⁱuʲσᵏ, and compute F = F_yⁱuʲσᵏ + Dᵢⱼₖ + Eᵢⱼₖ. Returns the right-hand sides −Σ_{s′} P[s,s′] F(s, s′)
   in Fs and pointers to the conditional G of each pair (null for inactive pairs).
   The pairs are independent here (each writes only its own G container and F; g** of other pairs is read), so
   they run as parallel jobs of the thread pool; inside a job, the tensor library's own thread_pool::run()
   calls run sequentially. */
std::vector<FGSTensor*>
KOrderMS::rhs(int i, int j, int k, std::vector<std::unique_ptr<FGSTensor>>& Fs)
{
  Symmetry sym {i, j, 0, k};
  JournalRecordPair pa(journal);
  pa << "MS: recovering symmetry " << sym << endrec;
  std::vector<FGSTensor*> Gc(S * S, nullptr);
  std::vector<std::unique_ptr<FGSTensor>> Fp(S * S);
  thread_pool::job_group_t jobs;
  for (int s = 0; s < S; s++)
    for (int sp = 0; sp < S; sp++)
      if (active(s, sp))
        jobs.push_back(std::make_unique<ClosureJob>([this, &Gc, &Fp, sym, i, j, k, s, sp] {
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
        }));
  thread_pool::run(jobs);
  Fs.clear();
  for (int s = 0; s < S; s++)
    {
      auto t = std::make_unique<FGSTensor>(*Fp[s * S + s]);
      t->zeros();
      for (int sp = 0; sp < S; sp++)
        if (active(s, sp))
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
  solveStackedS(X);
  insertMS(X);
  for (auto& p : pairs)
    p->G<Storage::fold>().insert(p->faaDiBrunoG<Storage::fold>(Symmetry {0, 0, 0, 1}));
}

void
KOrderMS::recover_ys(int i, int j)
{
  std::vector<std::unique_ptr<FGSTensor>> R;
  auto Gc = rhs(i, 0, j, R);
  KOrder& p0 = *pairs[0];
  if (ypart.nys() > 0 && ypart.nyss() > 0)
    {
      const int zc = ypart.nstat + ypart.npred, nyss = ypart.nyss();
      auto row = [&](int s, int r) { return r < zc ? s * zc + r : S * zc + s * nyss + r - zc; };
      std::vector<UGSTensor> ur;
      for (int s = 0; s < S; s++)
        ur.emplace_back(*R[s]);
      const int nc = ur[0].ncols();
      TwoDMatrix d(S * ny, nc);
      for (int s = 0; s < S; s++)
        for (int c = 0; c < nc; c++)
          for (int r = 0; r < ny; r++)
            d.get(s * ny + r, c) = ur[s].get(r, c);
      TwoDMatrix gs_y(p0.gs<Storage::unfold>().get(Symmetry {1, 0, 0, 0}));
      GeneralSylvester sylv(i, S * ny, ypart.nys(), S * zc, Abar.getData(), Bbar.getData(),
                            gs_y.getData(), d.getData());
      sylv.solve();
      for (int s = 0; s < S; s++)
        {
          for (int c = 0; c < nc; c++)
            for (int r = 0; r < ny; r++)
              ur[s].get(r, c) = d.get(row(s, r), c);
          FGSTensor ftmp(ur[s]);
          R[s]->getData() = const_cast<const Vector&>(ftmp.getData());
        }
    }
  else if (ypart.nys() > 0 && ypart.nyss() == 0)
    for (int s = 0; s < S; s++)
      p0.matA.multInv(*R[s]);
  insertMS(R);
  updateG(Gc, R, i + j);
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
  solveStackedS(R);
  insertMS(R);
  updateG(Gc, R, i);
}

/* g_yⁱuʲ and G_yⁱuʲ (no σ) involve only σ-free derivatives, which are the same in all pairs: insert copies of
   those of pair (0, 0) (g with its g* and g** row views, as KOrder::insertDerivative does, and G). */
void
KOrderMS::shareSigmaFree(const Symmetry& sym)
{
  KOrder& p0 = *pairs[0];
  const FGSTensor& g0 = p0.g<Storage::fold>().get(sym);
  const FGSTensor& G0 = p0.G<Storage::fold>().get(sym);
  for (int q = 1; q < S * S; q++)
    {
      KOrder& p = *pairs[q];
      p.insertDerivative<Storage::fold>(std::make_unique<FGSTensor>(g0));
      p.G<Storage::fold>().insert(std::make_unique<FGSTensor>(G0));
    }
}

/* The order of KOrder::performStep(); σ-free symmetries by the stock routines in pair (0, 0), copied to the
   other pairs. */
void
KOrderMS::performStep(int order)
{
  JournalRecordPair pa(journal);
  pa << "MS: performing step for order = " << order << endrec;
  pairs[0]->recover_y<Storage::fold>(order);
  shareSigmaFree(Symmetry {order, 0, 0, 0});
  for (int i = 0; i < order; i++)
    {
      pairs[0]->recover_yu<Storage::fold>(i, order - i);
      shareSigmaFree(Symmetry {i, order - i, 0, 0});
    }
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
