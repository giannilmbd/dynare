! Copyright © 2021-2026 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
!
! Pruned simulation of an order-k Markov-switching decision rule (matlab/ms_korder). Same algorithm as
! simulate_pruning in libkordersim/simulation.f08 (layer rule by integer partitions over the unfolded tensors
! [g_vˡ], v = (y*, u, σ)), with the tensors of period t taken from the regime of period t.
!
!  input:
!       order    the order of approximation
!       nstat, npred, nboth, nforw, nexog
!       ystart   starting value (full vector of endogenous, dr order)
!       shocks   matrix of shocks (nexog x number of periods)
!       ysteady  full vector of decision rule's steady state (dr order)
!       drs      struct array with one element per regime, fields g_1, …, g_order: unfolded fully symmetric
!                tensors [g_vˡ]/l! (as dr.pruning of k_order_perturbation), rows in dr order
!       regimes  regime (1 … number of regimes) of each period (1 x number of periods)
!  output:
!       res      simulated results (endo_nbr x (number of periods + 1), dr order; first column ystart)

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_c_binding
   use iso_fortran_env
   use matlab_mex
   use tensors
   use partitions
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in) :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: drs_mx, regimes_mx, ystart_mx, shocks_mx, ysteady_mx, g_q
   integer :: order, nstatic, npred, nboth, nfwrd, exo_nbr, endo_nbr, nys, nvar, nper, nreg
   real(real64), pointer, contiguous :: ysteady(:), ystart(:), sim(:,:), shocks(:,:), regimes(:)
   type(unfolded_tensor), dimension(:,:), allocatable :: g
   type(tensor), dimension(:), allocatable :: psim
   type(partition_triangle), dimension(:,:), allocatable, target :: part
   type(partition_triangle), pointer :: part_set
   real(real64), dimension(:,:), allocatable :: phat
   real(real64), dimension(:), contiguous, pointer :: u
   real(real64) :: prod, sum_part
   character(kind=c_char, len=10) :: fieldname
   integer :: q, l, j, m, r, t, c, s, rs

   if (nrhs /= 11 .or. nlhs /= 1) call mexErrMsgTxt("Must have exactly 11 inputs and 1 output")
   order = int(mxGetScalar(prhs(1)))
   nstatic = int(mxGetScalar(prhs(2)))
   npred = int(mxGetScalar(prhs(3)))
   nboth = int(mxGetScalar(prhs(4)))
   nfwrd = int(mxGetScalar(prhs(5)))
   exo_nbr = int(mxGetScalar(prhs(6)))
   ystart_mx = prhs(7)
   shocks_mx = prhs(8)
   ysteady_mx = prhs(9)
   drs_mx = prhs(10)
   regimes_mx = prhs(11)
   endo_nbr = nstatic+npred+nboth+nfwrd
   nys = npred+nboth
   nvar = nys+exo_nbr

   if (.not. (mxIsDouble(ystart_mx) .and. int(mxGetNumberOfElements(ystart_mx)) == endo_nbr)) &
        call mexErrMsgTxt("7th argument (ystart) should be a real vector with nstat+npred+nboth+nforw elements")
   if (.not. (mxIsDouble(shocks_mx) .and. int(mxGetM(shocks_mx)) == exo_nbr)) &
        call mexErrMsgTxt("8th argument (shocks) should be a real matrix with nexog rows")
   if (.not. (mxIsDouble(ysteady_mx) .and. int(mxGetNumberOfElements(ysteady_mx)) == endo_nbr)) &
        call mexErrMsgTxt("9th argument (ysteady) should be a real vector with nstat+npred+nboth+nforw elements")
   if (.not. mxIsStruct(drs_mx)) call mexErrMsgTxt("10th argument (drs) should be a struct array")
   nper = int(mxGetN(shocks_mx))
   if (.not. (mxIsDouble(regimes_mx) .and. int(mxGetNumberOfElements(regimes_mx)) == nper)) &
        call mexErrMsgTxt("11th argument (regimes) should be a real vector with one element per period")
   ystart => mxGetDoubles(ystart_mx)
   ysteady => mxGetDoubles(ysteady_mx)
   shocks(1:exo_nbr,1:nper) => mxGetDoubles(shocks_mx)
   regimes => mxGetDoubles(regimes_mx)
   nreg = int(mxGetNumberOfElements(drs_mx))
   do t=1,nper
      if (regimes(t) < 1._real64 .or. regimes(t) > real(nreg, real64) .or. regimes(t) /= aint(regimes(t))) &
           call mexErrMsgTxt("regimes should contain integers between 1 and the number of elements of drs")
   end do

   ! Import the decision rule tensors of every regime
   allocate(g(order, nreg))
   do s=1,nreg
      do q=1,order
         write (fieldname, '(a2, i1)') "g_", q
         g_q = mxGetField(drs_mx, int(s, mwIndex), trim(fieldname))
         if (.not. (c_associated(g_q) .and. mxIsDouble(g_q) .and. .not. mxIsComplex(g_q) .and. .not. mxIsSparse(g_q))) &
              call mexErrMsgTxt(trim(fieldname)//" is not allocated in drs")
         if (int(mxGetM(g_q)) /= endo_nbr .or. int(mxGetN(g_q)) /= (nvar+1)**q) &
              call mexErrMsgTxt(trim(fieldname)//" has wrong dimensions in drs")
         g(q,s)%m(1:mxGetM(g_q),1:mxGetN(g_q)) => mxGetDoubles(g_q)
      end do
      call fill_list_unfolded_tensor(g(:,s), nys, nvar+1)
   end do

   plhs(1) = mxCreateDoubleMatrix(int(endo_nbr, mwSize), int(nper+1, mwSize), mxREAL)
   sim(1:endo_nbr,1:(nper+1)) => mxGetDoubles(plhs(1))
   sim(:,1) = ystart

   allocate(psim(order), phat(nys,order), part(order,order))
   do l=1,order
      allocate(psim(l)%m(endo_nbr, 1))
      phat(:,l) = 0._real64
   end do
   phat(:,1) = ystart(nstatic+1:nstatic+nys)-ysteady(nstatic+1:nstatic+nys)
   call fill_partition_triangle(part)
   do t=2,nper+1
      u => shocks(:,t-1)
      rs = int(regimes(t-1))
      do l=1,order
         psim(l)%m(:,1) = 0._real64
      end do
      do l=1,order
         do j=1,size(g(l,rs)%ind)
            if (g(l,rs)%s(j) > 0) then
               do q=l,order
                  sum_part = 0._real64
                  part_set => part(g(l,rs)%s(j),q-(l-g(l,rs)%s(j)))
                  do r=1,size(part_set%partition)
                     prod = real(part_set%count(r),real64)
                     c = 1
                     do m=1,l
                        if (g(l,rs)%ind(j)%coor(m) <= nys) then
                           prod = prod*phat(g(l,rs)%ind(j)%coor(m), part_set%partition(r)%coor(c))
                           c = c+1
                        end if
                        if ((g(l,rs)%ind(j)%coor(m) > nys) .and. (g(l,rs)%ind(j)%coor(m) <= nvar)) &
                             prod = prod*u(g(l,rs)%ind(j)%coor(m)-nys)
                     end do
                     sum_part = sum_part+prod
                  end do
                  psim(q)%m(:,1) = psim(q)%m(:,1)+g(l,rs)%m(:,j)*sum_part
               end do
            else
               prod = 1._real64
               do m=1,l
                  if ((g(l,rs)%ind(j)%coor(m) > nys) .and. (g(l,rs)%ind(j)%coor(m) <= nvar)) &
                       prod = prod*u(g(l,rs)%ind(j)%coor(m)-nys)
               end do
               psim(l)%m(:,1) = psim(l)%m(:,1)+g(l,rs)%m(:,j)*prod
            end if
         end do
      end do
      sim(:,t) = ysteady
      do l=1,order
         phat(:,l) = psim(l)%m(nstatic+1:nstatic+nys,1)
         sim(:,t) = sim(:,t)+psim(l)%m(:,1)
      end do
   end do
end subroutine mexFunction
