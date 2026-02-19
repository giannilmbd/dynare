function oo_=PCL_Part_info_moments(M_, oo_, options_, varobs, dr, ivar)
% function oo_=PCL_Part_info_moments(M_, oo_, options_, varobs, dr, ivar)
% sets up parameters and calls part-info Kalman filter
% developed by G Perendia, July 2006 for implementation from notes by Prof. Joe Pearlman to
% suit partial information RE solution in accordance with, and based on, the
% Pearlman, Currie and Levine 1986 solution.
% 22/10/06 - Version 2 for new Riccati with 4 params instead 5

% Copyright © 2006-2025 Dynare Team
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

% Recall that the state space is given by the
% predetermined variables s(t-1), x(t-1)
% and the jump variables x(t).
% The jump variables have dimension NETA

if ~isfolder([M_.dname '/Output'])
    mkdir(M_.dname,'Output');
end

warning_old_state = warning;
warning off

OBS = [];
for i=1:rows(varobs)
    OBS = [OBS find(strcmp(varobs{i}, M_.endo_names)) ];
end
NOBS = length(OBS);

G1=dr.PI_ghx;
impact=dr.PI_ghu;
nmat=dr.PI_nmat;
NX=M_.exo_nbr; % no of exogenous varexo shock variables.
FL_RANK=dr.PI_FL_RANK;
NY=M_.endo_nbr;
LL = sparse(1:NOBS,OBS,ones(NOBS,1),NY,NY);

ss=size(G1,1);

pd=ss-size(nmat,1);
SDX=M_.Sigma_e^0.5; % =SD,not V-COV, of Exog shocks or M_.Sigma_e^0.5 num_exog x num_exog matrix
H=M_.H;
VV=H; % V-COV of observation errors.
MM=impact*SDX; % R*(Q^0.5) in standard KF notation
               % observation vector indices
               % mapping to endogenous variables.

L1=LL*dr.PI_TT1;
L2=LL*dr.PI_TT2;

MM1=MM(1:ss-FL_RANK,:);
U11=MM1*MM1';
% SDX
U22=0;
% determine K1 and K2 observation mapping matrices
% This uses the fact that measurements are given by L1*s(t)+L2*x(t)
% and s(t) is expressed in the dynamics as
% H1*eps(t)+G11*s(t-1)+G12*x(t-1)+G13*x(t).
% Thus the observations o(t) can be written in the form
% o(t)=K1*[eps(t)' s(t-1)' x(t-1)']' + K2*x(t) where
% K1=[L1*H1 L1*G11 L1*G12] K2=L1*G13+L2

G12=G1(NX+1:ss-2*FL_RANK,:);
KK1=L1*G12;
K1=KK1(:,1:ss-FL_RANK);
K2=KK1(:,ss-FL_RANK+1:ss)+L2;

%pre calculate time-invariant factors
A11=G1(1:pd,1:pd);
A22=G1(pd+1:end, pd+1:end);
A12=G1(1:pd, pd+1:end);
A21=G1(pd+1:end,1:pd);
BB=A12/A22;
FF=K2/A22;
QQ=BB*U22*BB' + U11;
UFT=U22*FF';
% kf_param structure:
AA=A11-BB*A21;
CCCC=A11-A12*nmat; % F in new notation
DD=K1-FF*A21; % H in new notation
RR=FF*UFT+VV;
if ~any(RR)
    % if zero add some dummy measurement err. variance-covariances
    % with diagonals 0.000001. This would not be needed if we used
    % the slow solver, or the generalised eigenvalue approach,
    % but these are both slower.
    RR=eye(size(RR,1))*1.0e-6;
end
SS=BB*UFT;
% initialise pshat
AQDS=AA*QQ*DD'+SS;
DQDR=DD*QQ*DD'+RR;
AQDQ=AQDS/DQDR;
ff=AA-AQDQ*DD;
hh=AA*QQ*AA'-AQDQ*AQDS';%*(DD*QQ*AA'+SS');
rr=DD*QQ*DD'+RR;
ZSIG0=disc_riccati_fast(ff,DD,rr,hh);
PP=ZSIG0 +QQ;

DPDR=DD*PP*DD'+RR;
PDIDPDRD=PP*DD'/DPDR*DD;

MSIG=lyapunov_solver(CCCC,CCCC,PDIDPDRD*PP,options_);

COV_P=[ PP, PP; PP, PP+MSIG]; % P0
COV_P=(COV_P+COV_P')/2; %assure symmetry

dr.PI_GG=[CCCC (AA-CCCC)*(eye(ss-FL_RANK)-PDIDPDRD); zeros(ss-FL_RANK) AA*(eye(ss-FL_RANK)-PDIDPDRD)];

GAM= [ AA*(eye(ss-FL_RANK)-PDIDPDRD) zeros(ss-FL_RANK); (AA-CCCC)*(eye(ss-FL_RANK)-PDIDPDRD),  CCCC];

VV = [  dr.PI_TT1 dr.PI_TT2];
nn=size(VV,1);
COV_OMEGA= COV_P(end-nn+1:end,end-nn+1:end);
COV_YR0= VV*COV_OMEGA*VV';
diagCovYR0=diag(COV_YR0);
diagCovYR0(abs(diagCovYR0)<1e-10)=0;
labels=get_labels_transformed_vars(M_.endo_names,ivar,options_,false);

if ~options_.nomoments
    z = [ sqrt(diagCovYR0(ivar)) diagCovYR0(ivar) ];
    title = 'THEORETICAL MOMENTS';
    headers = {'VARIABLE'; 'STD. DEV.'; 'VARIANCE'};
    dyntable(options_, title, headers, labels, z, cellofchararraymaxlength(labels)+2, 16, 10);
end
if ~options_.nocorr
    diagSqrtCovYR0 = sqrt(diagCovYR0);
    non_zero_indices=find(diagSqrtCovYR0>1e-10);
    temp=1./diagSqrtCovYR0(non_zero_indices);
    DELTA = zeros(length(diagCovYR0),1);
    DELTA(non_zero_indices)=temp;
    DELTA=diag(DELTA);
    COR_Y = DELTA*COV_YR0*DELTA;
    i1=intersect(non_zero_indices,ivar);
    if options_.contemporaneous_correlation
        oo_.contemporaneous_correlation = COR_Y(ivar,ivar);
    end
    if ~options_.noprint
        skipline()
        title = 'MATRIX OF CORRELATIONS';
        labels=get_labels_transformed_vars(M_.endo_names,i1,options_,false);
        headers = vertcat('Variables', labels);
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_, title, headers, labels, COR_Y(i1,i1), lh, 8, 4);
        if options_.TeX
            labels=get_labels_transformed_vars(M_.endo_names_tex,i1,options_,true);
            headers = vertcat('Variables', labels);
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'th_corr_matrix', headers, labels, COR_Y(i1,i1), lh, 8, 4);
        end
    end
else
    COR_Y=[];
end

ar = options_.ar;
if ar > 0
    AutoCOR_YRk= zeros(nn, ar);
    for k = 1:ar
        COV_P = GAM*COV_P;
        COV_OMEGA = COV_P( end-nn+1:end, end-nn+1:end);
        COV_YRk = VV*COV_OMEGA*VV';
        AutoCOR_YRkMAT = DELTA*COV_YRk*DELTA;
        oo_.autocorr{k} = AutoCOR_YRkMAT(ivar,ivar);
        AutoCOR_YRk(:,k) = diag(COV_YRk)./diagCovYR0;
    end
    title = 'COEFFICIENTS OF AUTOCORRELATION';
    labels=get_labels_transformed_vars(M_.endo_names,i1,options_,false);
    headers = vertcat('Order ', cellstr(int2str((1:options_.ar)')));
    lh = cellofchararraymaxlength(labels)+2;
    dyntable(options_, title, headers, labels, AutoCOR_YRk(i1,:), lh, 8, 4);
    if options_.TeX
        labels=get_labels_transformed_vars(M_.endo_names_tex,i1,options_,true);
        headers = vertcat('Order ', cellstr(int2str((1:options_.ar)')));
        lh = cellofchararraymaxlength(labels)+2;
        dyn_latex_table(M_, options_, title, 'th_autocorr_matrix', headers, labels, AutoCOR_YRk(i1,:), lh, 8, 4);
    end
else
    AutoCOR_YRk = [];
end
save ([M_.dname filesep 'Output' filesep M_.fname '_PCL_moments'], 'COV_YR0','AutoCOR_YRk', 'COR_Y');
warning(warning_old_state);
