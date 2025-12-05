! Copyright © 2025 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the Free Software Foundation, either version 3 of
! the License, or (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with Dynare. If not,
! see <https://www.gnu.org/licenses/>.
!
! Original author: Normann Rion <normann@dynare.org>

module interpolation
    use iso_fortran_env, only: real64, int32
    implicit none (type, external)
contains
    subroutine bracket_linear_weight(x, n, xq, nq, xqi, xqpi)
        !-----------------------------------------------------------------------
        ! Subroutine: bracket_linear_weight
        !
        ! Purpose:
        !   For a given monotonic increasing grid `x` and a vector of query
        !   points `xq`, computes for each query:
        !     - the index `xqi` such that xq ∈ [x(xqi), x(xqi+1)]
        !     - the linear interpolation weight `xqpi` for xqi
        !
        !   This allows for fast and robust piecewise linear interpolation.
        !
        ! Arguments:
        !   x     (in)  [real64, dimension(n)]
        !     Monotonically increasing 1D grid.
        !
        !   n     (in)  [int32]
        !     Number of elements in x (length of x).
        !
        !   xq    (in)  [real64, dimension(nq)]
        !     Query points to locate in the grid.
        !
        !   nq    (in)  [int32]
        !     Number of query points (length of xq).
        !
        !   xqi   (out) [int32, dimension(nq)]
        !     For each xq(i), the index `ilow` such that:
        !       x(ilow) ≤ xq(i) ≤ x(ilow+1)
        !
        !   xqpi  (out) [real64, dimension(nq)]
        !     For each xq(i), the relative weight in the bracketed interval:
        !       xqpi(i) = (x(ilow+1) - xq(i)) / (x(ilow+1) - x(ilow))
        !
        ! Notes:
        !   - Extrapolation is clamped to the first or last interval.
        !   - Grid `x` must be strictly increasing for meaningful output.
        !
        !-----------------------------------------------------------------------

        integer(int32), intent(in) :: n, nq
        real(real64), intent(in)  :: x(n), xq(nq)
        integer(int32), intent(out) :: xqi(nq)
        real(real64), intent(out)   :: xqpi(nq)

        integer(int32) :: iq, ilow, ihigh, imid

        do iq = 1, nq
            if (xq(iq) < x(1)) then
                ilow = 1
            else if (xq(iq) > x(n-1)) then
                ilow = n - 1
            else
                ilow = 1
                ihigh = n
                do while (ihigh - ilow > 1)
                    imid = (ihigh + ilow) / 2
                    if (xq(iq) > x(imid)) then
                        ilow = imid
                    else
                        ihigh = imid
                    end if
                end do
            end if
            xqi(iq) = ilow
            xqpi(iq) = (x(ilow+1) - xq(iq)) / (x(ilow+1) - x(ilow))
        end do
    end subroutine bracket_linear_weight
end module interpolation