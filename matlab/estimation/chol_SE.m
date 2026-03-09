function [R,indef, E, P]=chol_SE(A,pivoting)
% [R,indef, E, P]=chol_SE(A,pivoting)
% Performs a (modified) Cholesky factorization of the form
%
%     P'*A*P  + E = R'*R
%
% As detailed in Schnabel/Eskow (1990), the factorization has 2 phases:
%   Phase 1: While A can still be positive definite, pivot on the maximum diagonal element.
%            Check that the standard Cholesky update would result in a positive diagonal
%            at the current iteration. If so, continue with the normal Cholesky update.
%            Otherwise switch to phase 2.
%            If A is safely positive definite, stage 1 is never left, resulting in
%            the standard Cholesky decomposition.
%
%    Phase 2: Pivot on the minimum of the negatives of the lower Gershgorin bound
%            estimates. To prevent negative diagonals, compute the amount to add to the
%            pivot element and add this. Then, do the Cholesky update and update the estimates of the
%            Gershgorin bounds.
%
% Notes:
%   -   During factorization, L=R' is stored in the lower triangle of the original matrix A,
%       minimizing the memory requirements
%   -   Conforming with the original Schnabel/Eskow (1990) algorithm:
%            - at each iteration the updated Gershgorin bounds are estimated instead of recomputed,
%              reducing the computational requirements from o(n^3) to o (n^2)
%           -  For the last 2 by 2 submatrix, an eigenvalue-based decomposition is used
%   -   While pivoting is not necessary, it improves the size of E, the add-on on to the diagonal. But this comes at
%       the cost of introducing a permutation.
%
%
% INPUTS
%  - A           [n*n]     Matrix to be factorized
%  - pivoting    [scalar]  dummy whether pivoting is used
%
% OUTPUTS
%  - R           [n*n]     originally stored in lower triangular portion of A, including the main diagonal
%
%  - E           [n*1]     Elements added to the diagonal of A
%  - P           [n*1]     record of how the rows and columns of the matrix were permuted while
%                          performing the decomposition
%
% REFERENCES:
%   This implementation is based on
%
%       Robert B. Schnabel and Elizabeth Eskow. 1990. "A New Modified Cholesky
%       Factorization," SIAM Journal of Scientific Statistical Computating,
%       11, 6: 1136-58.
%
%       Elizabeth Eskow and Robert B. Schnabel 1991. "Algorithm 695 - Software for a New Modified Cholesky
%       Factorization," ACM Transactions on Mathematical Software, Vol 17, No 3: 306-312
%
%
% Author: Johannes Pfeifer based on Eskow/Schnabel (1991)

% Copyright © 2015 Johannes Pfeifer
% Copyright © 2015-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

if sum(sum(abs(A-A'))) > 1e-8
    error('A is not symmetric')
elseif sum(sum(abs(A-A'))) > 0
    A=(A+A')/2;
end

if nargin==1
    pivoting=0;
end


n=size(A,1);

tau1=eps^(1/3); %tolerance parameter for determining when to switch phase 2
tau2=eps^(1/3); %tolerance used for determining the maximum condition number of the final 2X2 submatrix.

phase1 = 1;
delta = 0;

P=1:n;
g=zeros(n,1);
E=zeros(n,1);

% Find the maximum magnitude of the diagonal elements. If any diagonal element is negative, then phase1 is false.
gammma=max(diag(A));
if any(diag(A) < 0)
    phase1 = 0;
end

taugam = tau1*gammma;
% If not in phase1, then calculate the initial Gershgorin bounds needed for the start of phase2.
if ~phase1
    g=gersh_nested(A,1,n);
end

% check for n=1
if n==1
    delta = tau2*abs(A(1,1)) - A(1,1);
    if delta > 0
        E(1) = delta;
    end
    if A(1,1) == 0
        E(1) = tau2;
    end
    A(1,1)=sqrt(A(1,1)+E(1));
end

for j = 1:n-1
    % PHASE 1
    if phase1
        if pivoting==1
            % Find index of maximum diagonal element A(i,i) where i>=j
            [~,imaxd] = max(diag(A(j:n,j:n)));
            imaxd=imaxd+j-1;
            % Pivot to the top the row and column with the max diag
            if (imaxd ~= j)
                % Swap row j with row of max diag
                for ii = 1:j-1
                    temp = A(j,ii);
                    A(j,ii) = A(imaxd,ii);
                    A(imaxd,ii) = temp;
                end
                % Swap colj and row maxdiag between j and maxdiag
                for ii = j+1:imaxd-1
                    temp = A(ii,j);
                    A(ii,j) = A(imaxd,ii);
                    A(imaxd,ii) = temp;
                end
                % Swap column j with column of max diag
                for ii = imaxd+1:n
                    temp = A(ii,j);
                    A(ii,j) = A(ii,imaxd);
                    A(ii,imaxd) = temp;
                end
                % Swap diag elements
                temp = A(j,j);
                A(j,j) = A(imaxd,imaxd);
                A(imaxd,imaxd) = temp;
                % Swap elements of the permutation vector
                itemp = P(j);
                P(j) = P(imaxd);
                P(imaxd) = itemp;
            end
        end
        % check to see whether the normal Cholesky update for this
        % iteration would result in a positive diagonal,
        % and if not then switch to phase 2.
        jp1 = j+1;
        if A(j,j)>0
            if min((diag(A(j+1:n,j+1:n)) - A(j+1:n,j).^2/A(j,j))) < taugam %test whether stage 2 is required
                phase1=0;
            end
        else
            phase1 = 0;
        end
        if phase1
            % Do the normal cholesky update if still in phase 1
            A(j,j) = sqrt(A(j,j));
            tempjj = A(j,j);
            for ii = jp1: n
                A(ii,j) = A(ii,j)/tempjj;
            end
            for ii=jp1:n
                temp=A(ii,j);
                for k = jp1:ii
                    A(ii,k) = A(ii,k)-(temp * A(k,j));
                end
            end
            if j == n-1
                A(n,n)=sqrt(A(n,n));
            end
        else
            % Calculate the negatives of the lower Gershgorin bounds
            g=gersh_nested(A,j,n);
        end
    end

    % PHASE 2
    if ~phase1
        if j ~= n-1
            if pivoting
                % Find the minimum negative Gershgorin bound
                [~,iming] = min(g(j:n));
                iming=iming+j-1;
                % Pivot to the top the row and column with the
                % minimum negative Gershgorin bound
                if iming ~= j
                    % Swap row j with row of min Gershgorin bound
                    for ii = 1:j-1
                        temp = A(j,ii);
                        A(j,ii) = A(iming,ii);
                        A(iming,ii) = temp;
                    end

                    % Swap col j with row iming from j to iming
                    for ii = j+1:iming-1
                        temp = A(ii,j);
                        A(ii,j) = A(iming,ii);
                        A(iming,ii) = temp;
                    end
                    % Swap column j with column of min Gershgorin bound
                    for ii = iming+1:n
                        temp = A(ii,j);
                        A(ii,j) = A(ii,iming);
                        A(ii,iming) = temp;
                    end
                    % Swap diagonal elements
                    temp = A(j,j);
                    A(j,j) = A(iming,iming);
                    A(iming,iming) = temp;
                    % Swap elements of the permutation vector
                    itemp = P(j);
                    P(j) = P(iming);
                    P(iming) = itemp;
                    % Swap elements of the negative Gershgorin bounds vector
                    temp = g(j);
                    g(j) = g(iming);
                    g(iming) = temp;
                end
            end
            % Calculate delta and add to the diagonal. delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
            % where normj=sum of |A(i,j)|,for i=1,n, delta_previous is the delta computed at the previous iter and taugam is tau1*gamma.

            normj=sum(abs(A(j+1:n,j)));

            delta = max([0;delta;-A(j,j)+normj;-A(j,j)+taugam]); % get adjustment based on formula on bottom of p. 309 of Eskow/Schnabel (1991)

            E(j) =  delta;
            A(j,j) = A(j,j) + E(j);
            % Update the Gershgorin bound estimates (note: g(i) is the negative of the Gershgorin lower bound.)
            if A(j,j) ~= normj
                temp = (normj/A(j,j)) - 1;
                for ii = j+1:n
                    g(ii) = g(ii) + abs(A(ii,j)) * temp;
                end
            end
            % Do the cholesky update
            A(j,j) = sqrt(A(j,j));
            tempjj = A(j,j);
            for ii = j+1:n
                A(ii,j) = A(ii,j) / tempjj;
            end
            for ii = j+1:n
                temp = A(ii,j);
                for k = j+1: ii
                    A(ii,k) = A(ii,k) - (temp * A(k,j));
                end
            end
        else
            % Find eigenvalues of final 2 by 2 submatrix
            % Find delta such that:
            % 1.  the l2 condition number of the final 2X2 submatrix + delta*I <= tau2
            % 2. delta >= previous delta,
            % 3. min(eigvals) + delta >= tau2 * gamma, where min(eigvals) is the smallest eigenvalue of the final 2X2 submatrix

            A(n-1,n)=A(n,n-1); %set value above diagonal for computation of eigenvalues
            eigvals  = eig(A(n-1:n,n-1:n));
            delta    = max([ 0 ; delta ; -min(eigvals)+tau2*max((max(eigvals)-min(eigvals))/(1-tau1),gammma) ]); %Formula 5.3.2 of Schnabel/Eskow (1990)

            if delta > 0
                A(n-1,n-1) = A(n-1,n-1) + delta;
                A(n,n) = A(n,n) + delta;
                E(n-1) = delta;
                E(n) = delta;
            end

            % Final update
            A(n-1,n-1) = sqrt(A(n-1,n-1));
            A(n,n-1) = A(n,n-1)/A(n-1,n-1);
            A(n,n) = A(n,n) - (A(n,n-1)^2);
            A(n,n) = sqrt(A(n,n));
        end
    end
end

R=(tril(A))';
indef=~phase1;
Pprod=zeros(n,n);
Pprod(sub2ind([n,n],P,1:n))=1;
P=Pprod;

function  g=gersh_nested(A,j,n)

g=zeros(n,1);
for ii = j:n
    if ii == 1
        sum_up_to_i = 0;
    else
        sum_up_to_i = sum(abs(A(ii,j:(ii-1))));
    end
    if ii == n
        sum_after_i = 0;
    else
        sum_after_i = sum(abs(A((ii+1):n,ii)));
    end
    g(ii) = sum_up_to_i + sum_after_i- A(ii,ii);
end

return

%@test:1
% Test 1: Positive definite matrix, no pivoting
% Standard Cholesky should result with E=0 and indef=0
A = [4 2 1; 2 5 3; 1 3 6];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, false);
t(2) = dassert(E, zeros(3,1));
t(3) = dassert(P, eye(3));
t(4) = dassert(R'*R, A, 1e-12);
% R should be upper triangular
t(5) = dassert(R, triu(R));
% Compare with MATLAB's built-in chol
R_matlab = chol(A);
t(6) = dassert(R, R_matlab, 1e-12);
T = all(t);
%@eof:1

%@test:2
% Test 2: Positive definite matrix, with pivoting
% Should still give E=0 and indef=0, but P may differ
A = [4 2 1; 2 5 3; 1 3 6];
[R, indef, E, P] = chol_SE(A, 1);
t(1) = dassert(indef, false);
t(2) = dassert(E, zeros(3,1));
% Verify factorization identity: P'*A*P + diag(E) = R'*R
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-12);
% R should be upper triangular
t(4) = dassert(R, triu(R));
% P should be a valid permutation matrix
t(5) = dassert(P*P', eye(3));
T = all(t);
%@eof:2

%@test:3
% Test 3: Indefinite matrix, no pivoting
% Should produce indef=1 and valid factorization with E>=0
A = [1 1 0; 1 1 1; 0 1 1];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, true);
% E should be non-negative
t(2) = dassert(all(E >= 0), true);
% Verify factorization identity: P'*A*P + diag(E) = R'*R
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-10);
% R should be upper triangular
t(4) = dassert(R, triu(R));
% R'*R should be positive definite (all eigenvalues > 0)
eigvals = eig(R'*R);
t(5) = dassert(all(eigvals > 0), true);
T = all(t);
%@eof:3

%@test:4
% Test 4: Indefinite matrix, with pivoting
A = [1 1 0; 1 1 1; 0 1 1];
[R, indef, E, P] = chol_SE(A, 1);
t(1) = dassert(indef, true);
t(2) = dassert(all(E >= 0), true);
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-10);
t(4) = dassert(R, triu(R));
t(5) = dassert(P*P', eye(3));
T = all(t);
%@eof:4

%@test:5
% Test 5: 1x1 positive scalar
A = [4];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, false);
t(2) = dassert(E, 0);
t(3) = dassert(R, 2, 1e-14);
T = all(t);
%@eof:5

%@test:6
% Test 6: 1x1 negative scalar
A = [-3];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, true);
t(2) = dassert(E(1) > 0, true);
t(3) = dassert(R^2, A + E(1), 1e-14);
T = all(t);
%@eof:6

%@test:7
% Test 7: 1x1 zero matrix
A = [0];
[R, indef, E, P] = chol_SE(A, 0);
tau2 = eps^(1/3);
t(1) = dassert(E(1), tau2);
t(2) = dassert(R, sqrt(tau2), 1e-14);
T = all(t);
%@eof:7

%@test:8
% Test 8: 2x2 indefinite matrix (exercises final 2x2 eigenvalue branch)
A = [1 2; 2 1];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, true);
t(2) = dassert(all(E >= 0), true);
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-12);
t(4) = dassert(R, triu(R));
T = all(t);
%@eof:8

%@test:9
% Test 9: Negative definite matrix
A = -[4 2 1; 2 5 3; 1 3 6];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, true);
t(2) = dassert(all(E >= 0), true);
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-10);
t(4) = dassert(R, triu(R));
eigvals = eig(R'*R);
t(5) = dassert(all(eigvals > 0), true);
T = all(t);
%@eof:9

%@test:10
% Test 10: Negative definite matrix with pivoting
A = -[4 2 1; 2 5 3; 1 3 6];
[R, indef, E, P] = chol_SE(A, 1);
t(1) = dassert(indef, true);
t(2) = dassert(all(E >= 0), true);
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-10);
t(4) = dassert(P*P', eye(3));
T = all(t);
%@eof:10

%@test:11
% Test 11: Non-symmetric matrix should error
A = [1 2; 3 4];
try
    [R, indef, E, P] = chol_SE(A, 0);
    t(1) = false;
catch
    t(1) = true;
end
T = all(t);
%@eof:11

%@test:12
% Test 12: Larger 5x5 indefinite matrix
A = [2 -1 0 0 0; -1 2 -1 0 0; 0 -1 0 -1 0; 0 0 -1 2 -1; 0 0 0 -1 2];
for piv = 0:1
    [R, indef, E, P] = chol_SE(A, piv);
    t(piv+1) = dassert(P'*A*P + diag(E), R'*R, 1e-10);
    t(piv+3) = dassert(R, triu(R));
    t(piv+5) = dassert(P*P', eye(5));
end
T = all(t);
%@eof:12

%@test:13
% Test 13: Nearly singular positive definite matrix
A = [1 0.999; 0.999 1];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(P'*A*P + diag(E), R'*R, 1e-12);
t(2) = dassert(R, triu(R));
T = all(t);
%@eof:13

%@test:14
% Test 14: Default pivoting argument (nargin==1)
A = [4 2; 2 5];
[R, indef, E, P] = chol_SE(A);
t(1) = dassert(indef, false);
t(2) = dassert(P, eye(2));
t(3) = dassert(R'*R, A, 1e-12);
T = all(t);
%@eof:14

%@test:15
% Test 15: Nearly symmetric matrix (small asymmetry gets symmetrized)
A = [4 2+1e-10; 2 5];
[R, indef, E, P] = chol_SE(A, 0);
t(1) = dassert(indef, false);
% Should succeed without error after symmetrization
t(2) = dassert(R, triu(R));
T = all(t);
%@eof:15

%@test:16
% Test 16: 4x4 positive definite matrix with pivoting
% Verify that pivoting produces correct result for well-conditioned PD matrix
A = [10 3 2 1; 3 8 1 2; 2 1 6 3; 1 2 3 9];
[R, indef, E, P] = chol_SE(A, 1);
t(1) = dassert(indef, false);
t(2) = dassert(E, zeros(4,1));
t(3) = dassert(P'*A*P + diag(E), R'*R, 1e-12);
t(4) = dassert(P*P', eye(4));
T = all(t);
%@eof:16
