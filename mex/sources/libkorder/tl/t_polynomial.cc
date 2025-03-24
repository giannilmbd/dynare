/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2025 Dynare Team
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

#include "t_polynomial.hh"

UTensorPolynomial::UTensorPolynomial(const FTensorPolynomial& fp) :
    TensorPolynomial<UFSTensor, UGSTensor, URSingleTensor>(fp.nrows(), fp.nvars())
{
  for (const auto& it : fp)
    insert(std::make_unique<UFSTensor>(*(it.second)));
}

FTensorPolynomial::FTensorPolynomial(const UTensorPolynomial& up) :
    TensorPolynomial<FFSTensor, FGSTensor, FRSingleTensor>(up.nrows(), up.nvars())
{
  for (const auto& it : up)
    insert(std::make_unique<FFSTensor>(*(it.second)));
}
