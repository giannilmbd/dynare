/*
 * Copyright © 2005 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Higher order at stochastic steady

/* This file defines a number of classes of which KOrderStoch is the main
   purpose. Basically, KOrderStoch calculates first and higher order Taylor
   expansion of a policy rule at σ>0 with explicit forward g**. More formally,
   we have to solve a policy rule g from

    𝔼ₜ[f(g**(g*(y*ₜ,uₜ,σ),uₜ₊₁,σ),g(y*,uₜ,σ),y*,uₜ)]

   As the introduction in approximation.hh argues, g** at tine t+1 must be
   given from outside. Let the explicit 𝔼ₜ(g**(y*,uₜ₊₁,σ) be equal to h(y*,σ).
   Then we have to solve f(h(g*(y*,u,σ),σ),g(y,u,σ),y,u), which is much easier
   than fully implicit system for σ=0.

   Besides the class KOrderStoch, we declare here also classes for the new
   containers corresponding to f(h(g*(y*,u,σ),σ),g(y,u,σ),y,u). Further, we
   declare IntegDerivs and StochForwardDerivs classes which basically calculate
   h as an extrapolation based on an approximation to g at lower σ.
*/

#ifndef KORDER_STOCH_HH
#define KORDER_STOCH_HH

#include <memory>

#include "faa_di_bruno.hh"
#include "journal.hh"
#include "korder.hh"
#include "pascal_triangle.hh"

/* This class is a container, which has a specialized constructor integrating
   the policy rule at given σ. */

template<Storage t>
class IntegDerivs : public ctraits<t>::Tgss
{
public:
  IntegDerivs(int r, const IntSequence& nvs, const typename ctraits<t>::Tgss& g,
              const typename ctraits<t>::Tm& mom, double at_sigma);
};

/* This constructor integrates a rule (namely its g** part) with respect to
   u=σ~·η, and stores to the object the derivatives of this integral h at
   (y*,u,σ)=(ỹ*,0,σ~). The original container of g**, the moments of the
   stochastic shocks ‘mom’ and the σ~ are input.

   The code follows the following derivation

   h(y,σ) = 𝔼ₜ[g(y,u′,σ)]

              1         ⎛  d  ⎞
    = ỹ + ∑  ──    ∑    ⎝i,j,k⎠ [g_yⁱuʲσᵏ] (y*-ỹ*)ⁱ σʲ Σʲ (σ-σ~)ᵏ
         ᵈ⁼¹ d! ⁱ⁺ʲ⁺ᵏ⁼ᵈ

              1          ⎛   d   ⎞                       ⎛m+n⎞
    = ỹ + ∑  ──     ∑    ⎝i,m+n,k⎠ [g_yⁱuᵐ⁺ⁿσᵏ] ŷ*ⁱ Σᵐ⁺ⁿ ⎝m,n⎠ σ~ᵐ σ^ᵏ⁺ⁿ
         ᵈ⁼¹ d! ⁱ⁺ᵐ⁺ⁿ⁺ᵏ⁼ᵈ

              1          ⎛   d   ⎞
    = ỹ + ∑  ──     ∑    ⎝i,m,n,k⎠ [g_yⁱuᵐ⁺ⁿσᵏ] Σᵐ⁺ⁿ σ~ᵐ ŷ*ⁱ σ^ᵏ⁺ⁿ
         ᵈ⁼¹ d! ⁱ⁺ᵐ⁺ⁿ⁺ᵏ⁼ᵈ

              1             ⎛   d   ⎞
    = ỹ + ∑  ──   ∑     ∑   ⎝i,m,n,k⎠ [g_yⁱuᵐ⁺ⁿσᵏ] Σᵐ⁺ⁿ σ~ᵐ ŷ*ⁱ σ^ᵏ⁺ⁿ
         ᵈ⁼¹ d! ⁱ⁺ᵖ⁼ᵈ  ᵐ⁼⁰
                      ⁿ⁺ᵏ⁼ᵖ

              1       ⎛ d ⎞ ⎡      ⎛ p ⎞  1                      ⎤
    = ỹ + ∑  ──   ∑   ⎝i,p⎠ ⎢  ∑   ⎝n,k⎠ ── [g_yⁱuᵐ⁺ⁿσᵏ] Σᵐ⁺ⁿ σ~ᵐ⎥ ŷ*ⁱ σ^ᵏ⁺ⁿ
         ᵈ⁼¹ d! ⁱ⁺ᵖ⁼ᵈ       ⎢ ᵐ⁼⁰        m!                      ⎥
                            ⎣ⁿ⁺ᵏ⁼ᵖ                               ⎦

         ⎛   a   ⎞
   where ⎝b₁,…,bₙ⎠ is a generalized combination number, p=k+n, σ^=σ-σ~,
   ŷ*=y*-ỹ, and we gave up writing the multidimensional indexes in Einstein
   summation.

   This implies that:

                   ⎛ p ⎞  1
    h_yⁱσᵖ =   ∑   ⎝n,k⎠ ── [g_yⁱuᵐ⁺ⁿσᵏ] Σᵐ⁺ⁿ σ~ᵐ
              ᵐ⁼⁰        m!
             ⁿ⁺ᵏ⁼ᵖ

   and this is exactly what the code does.
*/
template<Storage t>
IntegDerivs<t>::IntegDerivs(int r, const IntSequence& nvs, const typename ctraits<t>::Tgss& g,
                            const typename ctraits<t>::Tm& mom, double at_sigma) :
    ctraits<t>::Tgss(4)
{
  int maxd = g.getMaxDim();
  for (int d = 1; d <= maxd; d++)
    {
      for (int i = 0; i <= d; i++)
        {
          int p = d - i;
          Symmetry sym {i, 0, 0, p};
          auto ten = std::make_unique<typename ctraits<t>::Ttensor>(r, TensorDimens(sym, nvs));

          // Calculate derivative h_yⁱσᵖ
          /* This code calculates:

                             ⎛ p ⎞  1
              h_yⁱσᵖ =   ∑   ⎝n,k⎠ ── [g_yⁱuᵐ⁺ⁿσᵏ] Σᵐ⁺ⁿ σ~ᵐ
                        ᵐ⁼⁰        m!
                       ⁿ⁺ᵏ⁼ᵖ

             and stores it in ‘ten’. */
          ten->zeros();
          for (int n = 0; n <= p; n++)
            {
              int k = p - n;
              int povern = PascalTriangle::noverk(p, n);
              int mfac = 1;
              for (int m = 0; i + m + n + k <= maxd; m++, mfac *= m)
                {
                  double mult = (pow(at_sigma, m) * povern) / mfac;
                  Symmetry sym_mn {i, m + n, 0, k};
                  if (m + n == 0 && g.check(sym_mn))
                    ten->add(mult, g.get(sym_mn));
                  if (m + n > 0 && mom.check(Symmetry {m + n}) && g.check(sym_mn))
                    {
                      typename ctraits<t>::Ttensor gtmp(g.get(sym_mn));
                      gtmp.mult(mult);
                      gtmp.contractAndAdd(1, *ten, mom.get(Symmetry {m + n}));
                    }
                }
            }

          this->insert(std::move(ten));
        }
    }
}

/* This class calculates an extrapolation of expectation of forward
   derivatives. It is a container, all calculations are done in a constructor.

   The class calculates derivatives of E[g(y*,u,σ)] at (ȳ*,σ¯). The derivatives
   are extrapolated based on derivatives at (ỹ*,σ~). */

template<Storage t>
class StochForwardDerivs : public ctraits<t>::Tgss
{
public:
  StochForwardDerivs(const PartitionY& ypart, int nu, const typename ctraits<t>::Tgss& g,
                     const typename ctraits<t>::Tm& m, const Vector& ydelta, double sdelta,
                     double at_sigma);
};

/* This is the constructor which performs the integration and the
   extrapolation. Its parameters are: ‘g’ is the container of derivatives at
   (ỹ,σ~); ‘m’ are the moments of stochastic shocks; ‘ydelta’ is a difference
   of the steady states ȳ−ỹ; ‘sdelta’ is the difference between new sigma and
   old sigma σ¯−σ~, and ‘at_sigma’ is σ~. There is no need of inputing the ỹ.

   We do the operation in four steps:
   — Integrate g**, the derivatives are at (ỹ,σ~)
   — Form the (full symmetric) polynomial from the derivatives stacking
     ⎡y*⎤
     ⎣σ ⎦
   — Centralize this polynomial about (ȳ,σ¯)
   — Recover general symmetry tensors from the (full symmetric) polynomial
*/
template<Storage t>
StochForwardDerivs<t>::StochForwardDerivs(const PartitionY& ypart, int nu,
                                          const typename ctraits<t>::Tgss& g,
                                          const typename ctraits<t>::Tm& m, const Vector& ydelta,
                                          double sdelta, double at_sigma) :
    ctraits<t>::Tgss(4)
{
  int maxd = g.getMaxDim();
  int r = ypart.nyss();

  // Make ‘g_int’ be integral of g** at (ỹ,σ~)
  /* This simply constructs IntegDerivs class. Note that the ‘nvs’ of
     the tensors has zero dimensions for shocks, this is because we need to
                                    ⎡y*⎤
     make easily stacks of the form ⎣σ ⎦ in the next step. */
  IntSequence nvs {ypart.nys(), 0, 0, 1};
  IntegDerivs<t> g_int(r, nvs, g, m, at_sigma);

  // Make ‘g_int_sym’ be full symmetric polynomial from ‘g_int’
  /* Here we just form a polynomial whose unique variable corresponds to
     ⎡y*⎤
     ⎣σ ⎦ stack. */
  typename ctraits<t>::Tpol g_int_sym(r, ypart.nys() + 1);
  for (int d = 1; d <= maxd; d++)
    {
      auto ten = std::make_unique<typename ctraits<t>::Ttensym>(r, ypart.nys() + 1, d);
      ten->zeros();
      for (int i = 0; i <= d; i++)
        if (int k {d - i}; g_int.check(Symmetry {i, 0, 0, k}))
          ten->addSubTensor(g_int.get(Symmetry {i, 0, 0, k}));
      g_int_sym.insert(std::move(ten));
    }

  // Make ‘g_int_cent’ the centralized polynomial about (ȳ,σ¯)
  /* Here we centralize the polynomial to (ȳ,σ¯) knowing that the polynomial
     was centralized about (ỹ,σ~). This is done by derivating and evaluating
     the derivated polynomial at (ȳ−ỹ,σ¯-σ~). The stack of this vector is
     ‘delta’ in the code. */
  Vector delta(ypart.nys() + 1);
  Vector dy(delta, 0, ypart.nys());
  ConstVector dy_in(ydelta, ypart.nstat, ypart.nys());
  dy = dy_in;
  delta[ypart.nys()] = sdelta;
  typename ctraits<t>::Tpol g_int_cent(r, ypart.nys() + 1);
  for (int d = 1; d <= maxd; d++)
    {
      g_int_sym.derivative(d - 1);
      auto der = g_int_sym.evalPartially(d, delta);
      g_int_cent.insert(std::move(der));
    }

  // Pull out general symmetry tensors from ‘g_int_cent’
  /* Here we only recover the general symmetry derivatives from the full
     symmetric polynomial. Note that the derivative get the true ‘nvs’. */
  IntSequence ss {ypart.nys(), 0, 0, 1};
  IntSequence pp {0, 1, 2, 3};
  IntSequence true_nvs(nvs);
  true_nvs[1] = nu;
  true_nvs[2] = nu;
  for (int d = 1; d <= maxd; d++)
    if (g_int_cent.check(Symmetry {d}))
      for (int i = 0; i <= d; i++)
        {
          Symmetry sym {i, 0, 0, d - i};
          IntSequence coor(pp.unfold(sym));
          auto ten = std::make_unique<typename ctraits<t>::Ttensor>(
              g_int_cent.get(Symmetry {d}), ss, coor, TensorDimens(sym, true_nvs));
          this->insert(std::move(ten));
        }
}

/* This container corresponds to h(g*(y,u,σ),σ). Note that in our application,
   the σ as a second argument to h will be its fourth variable in symmetry, so
   we have to do four member stack having the second and third stack dummy. */

template<class _Ttype>
class GXContainer : public GContainer<_Ttype>
{
public:
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename StackContainer<_Ttype>::_Ctype;
  using itype = typename StackContainer<_Ttype>::itype;
  GXContainer(const _Ctype* gs, int ngs, int nu) : GContainer<_Ttype>(gs, ngs, nu)
  {
  }
  itype getType(int i, const Symmetry& s) const override;
};

/* This routine corresponds to this stack:
   ⎡ g*(y,u,σ) ⎤
   ⎢   dummy   ⎥
   ⎢   dummy   ⎥
   ⎣     σ     ⎦
*/
template<class _Ttype>
typename GXContainer<_Ttype>::itype
GXContainer<_Ttype>::getType(int i, const Symmetry& s) const
{
  if (i == 0)
    {
      if (s[2] > 0)
        return itype::zero;
      else
        return itype::matrix;
    }
  if (i == 1)
    return itype::zero;
  if (i == 2)
    return itype::zero;
  if (i == 3)
    {
      if (s == Symmetry {0, 0, 0, 1})
        return itype::unit;
      else
        return itype::zero;
    }

  KORD_RAISE("Wrong stack index in GXContainer::getType");
}

/* This container corresponds to f(H(y,u,σ),g(y,u,sigma),y,u), where the H has
   the size (number of rows) as g**. Since it is very simmilar to ZContainer,
   we inherit form it and override only getType() method. */

template<class _Ttype>
class ZXContainer : public ZContainer<_Ttype>
{
public:
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename StackContainer<_Ttype>::_Ctype;
  using itype = typename StackContainer<_Ttype>::itype;
  ZXContainer(const _Ctype* gss, int ngss, const _Ctype* g, int ng, int ny, int nu) :
      ZContainer<_Ttype>(gss, ngss, g, ng, ny, nu)
  {
  }
  itype getType(int i, const Symmetry& s) const override;
};

/* This getType() method corresponds to this stack:
   ⎡ H(y,u,σ) ⎤
   ⎢ g(y,u,σ) ⎥
   ⎢    y     ⎥
   ⎣    u     ⎦
*/
template<class _Ttype>
typename ZXContainer<_Ttype>::itype
ZXContainer<_Ttype>::getType(int i, const Symmetry& s) const
{
  if (i == 0)
    {
      if (s[2] > 0)
        return itype::zero;
      else
        return itype::matrix;
    }
  if (i == 1)
    {
      if (s[2] > 0)
        return itype::zero;
      else
        return itype::matrix;
    }
  if (i == 2)
    {
      if (s == Symmetry {1, 0, 0, 0})
        return itype::unit;
      else
        return itype::zero;
    }
  if (i == 3)
    {
      if (s == Symmetry {0, 1, 0, 0})
        return itype::unit;
      else
        return itype::zero;
    }

  KORD_RAISE("Wrong stack index in ZXContainer::getType");
}

class UnfoldedGXContainer : public GXContainer<UGSTensor>, public UnfoldedStackContainer
{
public:
  using _Ctype = TensorContainer<UGSTensor>;
  UnfoldedGXContainer(const _Ctype* gs, int ngs, int nu) : GXContainer<UGSTensor>(gs, ngs, nu)
  {
  }
};

class FoldedGXContainer : public GXContainer<FGSTensor>, public FoldedStackContainer
{
public:
  using _Ctype = TensorContainer<FGSTensor>;
  FoldedGXContainer(const _Ctype* gs, int ngs, int nu) : GXContainer<FGSTensor>(gs, ngs, nu)
  {
  }
};

class UnfoldedZXContainer : public ZXContainer<UGSTensor>, public UnfoldedStackContainer
{
public:
  using _Ctype = TensorContainer<UGSTensor>;
  UnfoldedZXContainer(const _Ctype* gss, int ngss, const _Ctype* g, int ng, int ny, int nu) :
      ZXContainer<UGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

class FoldedZXContainer : public ZXContainer<FGSTensor>, public FoldedStackContainer
{
public:
  using _Ctype = TensorContainer<FGSTensor>;
  FoldedZXContainer(const _Ctype* gss, int ngss, const _Ctype* g, int ng, int ny, int nu) :
      ZXContainer<FGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

/* This matrix corresponds to

    [f_{y}]+ [0 [f_y**₊]·[h**_y*] 0]

   This is almost the same as MatrixA, the only difference that the MatrixA is
   constructed from whole h_y*, not only from h**_y*, hence the new
   abstraction. */

class MatrixAA : public PLUMatrix
{
public:
  MatrixAA(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gyss,
           const PartitionY& ypart);
};

/* This class calculates derivatives of g given implicitly by
   f(h(g*(y,u,σ),σ),g(y,u,σ),y,u), where h(y,σ) is given from outside.

   Structurally, the class is very similar to KOrder, but calculations are much
   easier. The two constructors construct an object from sparse derivatives of
   f, and derivatives of h. The caller must ensure that the both derivatives
   are done at the same point.

   The calculation for order k (including k=1) is done by a call
   performStep(k). The derivatives can be retrived by getFoldDers() or
   getUnfoldDers(). */

class KOrderStoch
{
protected:
  IntSequence nvs;
  PartitionY ypart;
  Journal& journal;
  UGSContainer _ug;
  FGSContainer _fg;
  UGSContainer _ugs;
  FGSContainer _fgs;
  UGSContainer _uG;
  FGSContainer _fG;
  const UGSContainer* _uh;
  const FGSContainer* _fh;
  UnfoldedZXContainer _uZstack;
  FoldedZXContainer _fZstack;
  UnfoldedGXContainer _uGstack;
  FoldedGXContainer _fGstack;
  const TensorContainer<FSSparseTensor>& f;
  MatrixAA matA;

public:
  KOrderStoch(const PartitionY& ypart, int nu, const TensorContainer<FSSparseTensor>& fcont,
              const FGSContainer& hh, Journal& jr);
  KOrderStoch(const PartitionY& ypart, int nu, const TensorContainer<FSSparseTensor>& fcont,
              const UGSContainer& hh, Journal& jr);
  template<Storage t>
  void performStep(int order);
  [[nodiscard]] const FGSContainer&
  getFoldDers() const
  {
    return _fg;
  }
  [[nodiscard]] const UGSContainer&
  getUnfoldDers() const
  {
    return _ug;
  }
  [[nodiscard]] const FGSContainer&
  getFoldDersS() const
  {
    return _fgs;
  }

protected:
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoZ(const Symmetry& sym) const;
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoG(const Symmetry& sym) const;

  // Convenience access methods
  template<Storage t>
  typename ctraits<t>::Tg& g();
  template<Storage t>
  const typename ctraits<t>::Tg& g() const;

  template<Storage t>
  typename ctraits<t>::Tgs& gs();
  template<Storage t>
  const typename ctraits<t>::Tgs& gs() const;

  template<Storage t>
  const typename ctraits<t>::Tgss& h() const;

  template<Storage t>
  typename ctraits<t>::TG& G();
  template<Storage t>
  const typename ctraits<t>::TG& G() const;

  template<Storage t>
  typename ctraits<t>::TZXstack& Zstack();
  template<Storage t>
  const typename ctraits<t>::TZXstack& Zstack() const;

  template<Storage t>
  typename ctraits<t>::TGXstack& Gstack();
  template<Storage t>
  const typename ctraits<t>::TGXstack& Gstack() const;
};

/* This calculates a derivative of f(G(y,u,σ),g(y,u,σ),y,u) of a given
   symmetry. */

template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrderStoch::faaDiBrunoZ(const Symmetry& sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faà Di Bruno ZX container for " << sym << endrec;
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(ypart.ny(), TensorDimens(sym, nvs));
  FaaDiBruno bruno(journal);
  bruno.calculate(Zstack<t>(), f, *res);
  return res;
}

/* This calculates a derivative of G(y,u,σ)=h(g*(y,u,σ),σ) of a given
   symmetry. */

template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrderStoch::faaDiBrunoG(const Symmetry& sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faà Di Bruno GX container for " << sym << endrec;
  TensorDimens tdims(sym, nvs);
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(ypart.nyss(), tdims);
  FaaDiBruno bruno(journal);
  bruno.calculate(Gstack<t>(), h<t>(), *res);
  return res;
}

/* This retrieves all g derivatives of a given dimension from implicit
   f(h(g*(y,u,σ),σ),g(y,u,σ),y,u). It supposes that all derivatives of smaller
   dimensions have been retrieved.

   So, we go through all symmetries s, calculate Gₛ conditional on gₛ=0, insert
   the derivative to the G container, then calculate Fₛ conditional on gₛ=0.
   This is a righthand side. The left hand side is matA·gₛ. The gₛ is retrieved
   as gₛ=-matA⁻¹·RHS. Finally we have to update Gₛ by calling
   Gstack<t>().multAndAdd(1, h<t>(), *G_sym_ptr). */

template<Storage t>
void
KOrderStoch::performStep(int order)
{
  int maxd = g<t>().getMaxDim();
  KORD_RAISE_IF(order - 1 != maxd && (order != 1 || maxd != -1),
                "Wrong order for KOrderStoch::performStep");
  for (auto& si : SymmetrySet(order, 4))
    if (si[2] == 0)
      {
        JournalRecordPair pa(journal);
        pa << "Recovering symmetry " << si << endrec;

        auto G_sym = faaDiBrunoG<t>(si);
        auto G_sym_ptr = G_sym.get();
        G<t>().insert(std::move(G_sym));

        auto g_sym = faaDiBrunoZ<t>(si);
        auto g_sym_ptr = g_sym.get();
        g_sym->mult(-1.0);
        matA.multInv(*g_sym);
        g<t>().insert(std::move(g_sym));
        gs<t>().insert(
            std::make_unique<typename ctraits<t>::Ttensor>(ypart.nstat, ypart.nys(), *g_sym_ptr));

        Gstack<t>().multAndAdd(1, h<t>(), *G_sym_ptr);
      }
}

#endif
