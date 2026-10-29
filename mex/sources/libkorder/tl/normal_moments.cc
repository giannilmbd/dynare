/*
 * Copyright © 2004 Ondra Kamenik
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

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>

#include "kron_prod.hh"
#include "normal_moments.hh"
#include "permutation.hh"
#include "tl_static.hh"

std::vector<double>
skewNormalRawMoments(double variance, double skewness, int maxdim)
{
  if (variance < 0)
    TL_RAISE("Variance of a skew normal variable must be non-negative");
  constexpr double pi = std::numbers::pi;
  const double c = std::sqrt(2 / pi);

  double delta = 0;
  if (skewness != 0)
    {
      const double s23 = std::pow(std::abs(skewness), 2.0 / 3.0);
      delta = std::copysign(std::sqrt(pi / 2 * s23 / (s23 + std::pow((4 - pi) / 2, 2.0 / 3.0))),
                            skewness);
    }
  if (std::abs(delta) >= 1)
    TL_RAISE("Skewness of a skew normal variable must be smaller than 0.9953 in absolute value");
  const double omega = std::sqrt(variance / (1 - 2 / pi * delta * delta));
  const double sd1 = std::sqrt(1 - delta * delta);

  // binom[n][k] = C(n,k)
  std::vector<std::vector<double>> binom(maxdim + 1, std::vector<double>(maxdim + 1, 0));
  for (int n = 0; n <= maxdim; n++)
    {
      binom[n][0] = 1;
      for (int k = 1; k <= n; k++)
        binom[n][k] = binom[n - 1][k - 1] + (k < n ? binom[n - 1][k] : 0);
    }

  // E|U₀|ⁱ = 2^{i/2}·Γ((i+1)/2)/√π, and E[(√(1−δ²)·U₁)ⁿ] = (1−δ²)^{n/2}·(n−1)!! for even n
  std::vector<double> abs_mom(maxdim + 1), mom_b(maxdim + 1, 0);
  for (int i = 0; i <= maxdim; i++)
    abs_mom[i] = std::pow(2, i / 2.0) * std::tgamma((i + 1) / 2.0) / std::sqrt(pi);
  mom_b[0] = 1;
  for (int n = 2; n <= maxdim; n += 2)
    mom_b[n] = mom_b[n - 2] * (n - 1) * sd1 * sd1;

  // E[(δ·(|U₀|−√(2/π)))ʲ]
  std::vector<double> mom_a(maxdim + 1, 0);
  for (int j = 0; j <= maxdim; j++)
    {
      double sum = 0;
      for (int i = 0; i <= j; i++)
        sum += binom[j][i] * abs_mom[i] * std::pow(-c, j - i);
      mom_a[j] = std::pow(delta, j) * sum;
    }

  std::vector<double> res(maxdim + 1, 0);
  for (int d = 0; d <= maxdim; d++)
    {
      double sum = 0;
      for (int j = 0; j <= d; j++)
        sum += binom[d][j] * mom_a[j] * mom_b[d - j];
      res[d] = std::pow(omega, d) * sum;
    }
  // The mean is zero by construction
  if (maxdim >= 1)
    res[1] = 0;
  return res;
}

UNormalMoments::UNormalMoments(int maxdim, const TwoDMatrix& v) : TensorContainer<URSingleTensor>(1)
{
  if (maxdim >= 2)
    generateMoments(maxdim, v);
}

UNormalMoments::UNormalMoments(int maxdim, const std::vector<std::vector<double>>& raw_moments) :
    TensorContainer<URSingleTensor>(1)
{
  if (maxdim >= 2)
    generateIndependentMoments(maxdim, raw_moments);
}

/* For independent components, E[u_{i₁}·…·u_{i_d}] is the product over the
   components i of E[uᵢ^{nᵢ}], where nᵢ is the number of occurrences of i in
   the index (i₁,…,i_d). */
void
UNormalMoments::generateIndependentMoments(int maxdim,
                                           const std::vector<std::vector<double>>& raw_moments)
{
  int nv = raw_moments.size();
  for (const auto& m : raw_moments)
    if (static_cast<int>(m.size()) <= maxdim)
      TL_RAISE("Not enough raw moments in UNormalMoments constructor");

  std::vector<int> count(nv);
  for (int d = 2; d <= maxdim; d++)
    {
      auto mom = std::make_unique<URSingleTensor>(nv, d);
      bool nonzero = false;
      for (Tensor::index it = mom->begin(); it != mom->end(); ++it)
        {
          std::fill(count.begin(), count.end(), 0);
          const IntSequence& coor = it.getCoor();
          for (int m = 0; m < d; m++)
            count[coor[m]]++;
          double val = 1;
          for (int i = 0; i < nv && val != 0; i++)
            val *= raw_moments[i][count[i]];
          mom->get(*it, 0) = val;
          nonzero = nonzero || val != 0;
        }
      if (nonzero)
        insert(std::move(mom));
    }
}

/* Here we fill up the container with the tensors for d=2,4,6,… up to the given
   dimension. Each tensor of moments is equal to Fₙ(⊗ⁿv). This has a dimension
   equal to 2n. See the header file for proof and details.

   Here we sequentially construct the Kronecker power ⊗ⁿv and apply Fₙ.
*/
void
UNormalMoments::generateMoments(int maxdim, const TwoDMatrix& v)
{
  TL_RAISE_IF(v.nrows() != v.ncols(),
              "Variance-covariance matrix is not square in UNormalMoments constructor");

  int nv = v.nrows();
  auto mom2 = std::make_unique<URSingleTensor>(nv, 2);
  mom2->getData() = v.getData();
  insert(std::move(mom2));
  auto kronv = std::make_unique<URSingleTensor>(nv, 2);
  kronv->getData() = v.getData();
  for (int d = 4; d <= maxdim; d += 2)
    {
      auto newkronv = std::make_unique<URSingleTensor>(nv, d);
      KronProd::kronMult(ConstVector(v.getData()), ConstVector(kronv->getData()),
                         newkronv->getData());
      kronv = std::move(newkronv);
      auto mom = std::make_unique<URSingleTensor>(nv, d);
      // apply Fₙ to ‘kronv’
      /* Here we go through all equivalences, select only those having 2
         elements in each class, then go through all elements in ‘kronv’ and
         add to permuted location of ‘mom’.

         The permutation must be taken as inverse of the permutation implied by
         the equivalence, since we need a permutation which after application
         to identity of indices yileds indices in the equivalence classes. Note
         how the Equivalence::apply() method works. */
      mom->zeros();
      const EquivalenceSet eset = TLStatic::getEquiv(d);
      for (const auto& cit : eset)
        if (selectEquiv(cit))
          {
            Permutation per(cit);
            per.inverse();
            for (Tensor::index it = kronv->begin(); it != kronv->end(); ++it)
              {
                IntSequence ind(kronv->dimen());
                per.apply(it.getCoor(), ind);
                Tensor::index it2(*mom, ind);
                mom->get(*it2, 0) += kronv->get(*it, 0);
              }
          }
      insert(std::move(mom));
    }
}

/* We return true for an equivalence whose each class has 2 elements. */

bool
UNormalMoments::selectEquiv(const Equivalence& e)
{
  if (2 * e.numClasses() != e.getN())
    return false;
  for (const auto& si : e)
    if (si.length() != 2)
      return false;
  return true;
}

/* Here we go through all the unfolded container, fold each tensor and
   insert it. */
FNormalMoments::FNormalMoments(const UNormalMoments& moms) : TensorContainer<FRSingleTensor>(1)
{
  for (const auto& mom : moms)
    insert(std::make_unique<FRSingleTensor>(*(mom.second)));
}
