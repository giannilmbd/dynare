var r, pi, n, w, y, mc, x, rr, yobs, piobs,robs;
      
var y_star, n_star, w_star;
var  y_f,pi_f,y_e,pi_e, pi_exp,y_exp;
%Exogenous process (common to both specifications) 
var a, nu, u;
%Shocks (common to both specifications) 
varexo eps_a,eps_nu,eps_u;
%model parameters (common to both specifications) 
parameters h, gama, phi, bet, var_omega, kappa, iota_r, iota_y, iota_pi,iota_piexp;
%shock persistence (common to both specifications) 
parameters rho_a,rho_nu,rho_u;
%some steady state variables (common to both specifications) 
parameters  a_ss, n_ss, mc_ss;
% learning parameters
parameters theta, rho,sigmaScale;

parameters alpha_pi_e, alpha_pi_f, alpha_y_f, alpha_y_e;

%----------------------------------------------------------------
% Parametrization
%----------------------------------------------------------------
%model parameters
h=0.5; %%%%%%%%%%to be matched with moments  habit (0,1) 
gama=1; %%%%%%%%%%to be matched with moments  inverse of the inter. elasticity of subst. 
phi=1;%%%%%%%%%%to be matched with moments %inverse of the Frisch 
bet=0.99; %%%%%%%%%%set
var_omega=0.75; %%%%%%%%%%to be matched with moments   Calvo's probability 
kappa=0.5;%%%%%%%%%%to be matched with moments   Past indexation (0,1) 
%labor 
n_ss=0.2;
%elasticity of substituion across varieties sigma
sigma=6;
%marginal cost
mc_ss=(sigma-1)/sigma;
%profit_share
 profit_share=1/sigma;
 wage_share=1-profit_share;
%use the demand for labor
a_ss=(wage_share)/(mc_ss*n_ss);
% % Taylor rule parameters
iota_r=0.2; %set %%%%%%%%%% to be estimated
iota_y=0.5/4; %set%%%%%%%%%% to be estimated
iota_pi=1.5;  %set%%%%%%%%%% to be estimated
iota_piexp=0; %set%%%%%%%%%% to be estimated
rho_a=0.8; %%%%%%%%%%to be matched with moments persistence tech shock 
rho_nu=0.8;   %%%%%%%%%%to be matched with moments persistence demand shock 
rho_u=0.5; %%%%%%%%%%to be matched with moments existence monetary policy shock (usually lower than the other two)
%Behavioral Parameters
theta=2; %%%%%%%%%% to be estimated
rho=0.5; %%%%%%%%%% to be estimated
sigmaScale=0.005; %%%%%%%%%% Scale Parameter of the logistic function
alpha_pi_e=0.5; 
alpha_pi_f=0.5; 
alpha_y_f=0.5; 
alpha_y_e=0.5; 

%----------------------------------------------------------------
% The Model 
%----------------------------------------------------------------
model; 
//Composite parameters 
#messy1=((1-var_omega*bet)/(1+bet*kappa))*(1-var_omega)/var_omega; % I coeff. NKPC
#messy2=bet/(1+bet*kappa); % II coeff. NKPC
#messy3=kappa/(1+bet*kappa); % III coeff. NKPC
#epsilon_n=(a_ss*n_ss); %elasticity wrt n
%sticky prices
[name='Rep Euler Equation']
y*(1+h)=y(+1)+h*y(-1)-((1-h)/gama)*(r-pi(+1))+nu;
%eq2
[name='Labor supply']
n=(w-gama*y/(1-h)+gama*h*y(-1)/(1-h))/phi;
%eq3
[name='Production function']
y=epsilon_n*(a+n);
%eq4
[name='Labor Demand']
w=mc+(y-n)+a;
[name='Rep NKPC']
pi=messy1*mc+messy2*pi(+1)+messy3*pi(-1);
%eq9
[name='Taylor Rule']
r=(1-iota_r)*(iota_pi*pi+iota_y*(y-y_star)+iota_piexp*pi_exp)+iota_r*r(-1)+u;
%eq10
[name='Output Gap']
x=y-y_star;
%%%%Flex prices
%eq11
[name='Labor Supply - Flex Prices']
n_star=(w_star-gama*y_star/(1-h)+gama*h*y_star(-1)/(1-h))/phi;
%eq12
[name='Production Function - Flex Prices']
y_star=epsilon_n*(a+n_star);
%eq13
[name='Labor Demand - Flex Prices']
w_star=y_star-n_star+a;
% % % % % % % % Equation swap
%eq17
 [name='Rational Expectation Y']
 y_f=y(+1);
% %eq18
 [name='Rational Expectation PI']
 pi_f=pi(+1);
% %  Fundamentalist
%eq19
[name='Extrapolative Expectation Y']
y_e=y(-1);
%eq20
[name='Extrapolative Expectation PI']
pi_e=pi(-1);
%eq21
[name='Mkt Expectation Y']
y_exp=alpha_y_f*y_f+alpha_y_e*y_e;
%eq22
[name='Mkt Expectation PI']
pi_exp=alpha_pi_f*pi_f+alpha_pi_e*pi_e;
%eq31
[name='Real Interest Rate']
rr=r-pi(+1);
% Observation Equations
yobs=y;
piobs=pi;
robs=r;
%Exogenous process
%eq33-36
[name='TFP Shock']
a=rho_a*a(-1)+eps_a;
[name='Demand Schock']
nu=rho_nu*nu(-1)+eps_nu;
[name='Monetary Policy Shock']
u=rho_u*u(-1)+eps_u;
end;

steady_state_model; 
mc=0;
n=0; 
w=0;
pi=0; 
piobs=0;                                                                              
r=0;
robs=0;                                                                                  
y=0;
yobs=0;
x=0; 
rr=0;
n_star=0; 
y_star=0; 
w_star=0; 
a=0; 
nu=0; 
u=0; 
y_f=0;
pi_f=0;
y_e=0;
pi_e=0; 
pi_exp=0;
y_exp=0;
end;


steady;
check;

              
              
       M_.Sigma_e=eye(M_.exo_nbr)*0;
     shocks;
      var eps_a; stderr 0.01;
        end;
     
       
      varobs yobs piobs robs  ;
 
       
  
stoch_simul(order = 1,irf=40, periods=0, partial_information,contemporaneous_correlation,TeX)y,x,pi,r,rr;
collect_latex_files;

[status, cmdout]=system(['pdflatex -halt-on-error -interaction=nonstopmode ' M_.fname '_TeX_binder.tex']);
if status
    cmdout
    error('TeX-File did not compile.')
end
