/*
 * This file provides replication files for 
 * Smets, Frank and Wouters, Rafael (2007): "Shocks and Frictions in US Business Cycles: A Bayesian
 * DSGE Approach", American Economic Review, 97(3), 586-606, that are compatible with Dynare 4.5 onwards
 *
 * To replicate the full results, you have to get back to the original replication files available at
 * https://www.aeaweb.org/articles.php?doi=10.1257/aer.97.3.586 and include the respective estimation commands and mode-files.
 *
 * Notes:
 *  - The consumption Euler equation in the paper, equation (2), premultiplies the risk premium process \varepsilon_t^b,
 *      denoted by b in this code, by the coefficient c_3. In the code this prefactor is omitted by setting the 
 *      coefficient to 1. As a consequence, b in this code actually is b:=c_3*\varepsilon_t^b. As a consequence, in 
 *      the arbitrage equation for the value of capital in the paper, equation (4), the term 1*\varepsilon_t^b
 *      is replaced by 1/c_3*b, which is equal to \varepsilon_t^b given the above redefinition. This rescaling also explains why the 
 *      standard deviation of the risk premium shock in the AR(1)-process for b has a different standard deviation than reported
 *      in the paper. However, the results are unaffected by this scaling factor (except for the fact that the posterior distribution
 *      reported in the paper cannot be directly translated to the present mod-file due to parameter correlation in the posterior.  
 *  - As pointed out in Del Negro/Schorfheide (2012): "Notes on New-Keynesian Models"
 *      in the code implementation of equation (8) for both the flex price and the sticky price/wage economy, 
 *      there is a (1+cbetabar*cgamma) missing in the i_2 in front of q_t (denoted qs in the code). 
 *      Equation (8) in the paper reads:  
 *          (1-(1-delta)/gamma)*(1+beta*gamma^(1-sigma))*gamma^2*varphi
 *      which translates to the code snippet:
 *          (1-(1-ctou)/cgamma)*(1+cbetabar*cgamma)*cgamma^2*csadjcost
 *      But the code implements
 *          (1-(1-ctou)/cgamma)*cgamma^2*csadjcost
 *      which corresponds to an equation reading
 *          (1-(1-delta)/gamma)*gamma^2*varphi
 *  - Chib/Ramamurthy (2010): "Tailored randomized block MCMC methods with application to DSGE models", Journal of Econometrics, 155, pp. 19-38
 *      have pointed out that the mode reported in the original Smets/Wouters (2007) paper is not actually the mode. \bar \pi (constepinf) is estimated lower
 *      while \bar \l (constelab) is higher.
 *  - Note that at the prior mean, [cmap,crhopinf] and [cmaw,crhow] are pairwise collinear. Thus, running identification at the prior
 *      mean will return a warning. But this is only a local issue. These parameters are only indistinguishable at the prior mean, but not 
 *      at different points.
 *  - In the prior Table 1A in the paper, the 
 *          - habit parameter $\lambda$ is erroneously labeled h
 *          - the fixed cost parameter $\phi_p$ is labeled $\Phi$ 
 *  - Table 1B claims that $\rho_{ga}$ follows a beta prior with B(0.5,0.2^2), but the code shows that it actually
 *      follows a normal distribution with N(0.5,0.25^2)
 *
 * This file was originally written by Frank Smets and Rafeal Wouters and has been updated by
 * Johannes Pfeifer. 
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright © 2007-2013 Frank Smets and Raf Wouters
 * Copyright © 2013-2025 Johannes Pfeifer
 * Copyright © 2025-2026 Dynare Team
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */
// @#define observe_levels = 1
// @#define set_options_occbin = set_options_occbin + ["likelihood.first_period_binding_regime_allowed=101"]
// @#if dgp_estim
// 	@#define data_file = "usdata_pkf"//usmodel_data"//usdata_pkf"
// @#endif
// -----------------Dynare ----------------------------------------------//   
var labobs      ${lHOURS}$      (long_name='log hours worked') 
    robs        ${FEDFUNDS}$    (long_name='Federal funds rate') 
    pinfobs     ${dlP}$         (long_name='Inflation') 
    dy          ${dlGDP}$       (long_name='Output growth rate') 
    dc          ${dlCONS}$      (long_name='Consumption growth rate') 
    dinve       ${dlINV}$       (long_name='Investment growth rate') 
    dw          ${dlWAG}$       (long_name='Wage growth rate') 
    ewma        ${\eta^{w,aux}}$ (long_name='Auxiliary wage markup moving average variable')  
    epinfma     ${\eta^{p,aux}}$ (long_name='Auxiliary price markup moving average variable')
        ygap        ${ygap}$        (long_name='ygap') 
    mc          ${\mu_p}$       (long_name='gross price markup') 
    zcap        ${z}$           (long_name='Capital utilization rate') 
    rk          ${r^{k}}$       (long_name='rental rate of capital') 
    k           ${k^{s}}$       (long_name='Capital services') 
    pk          ${q}$           (long_name='real value of existing capital stock') 
    c           ${c}$           (long_name='Consumption')
    inve        ${i}$           (long_name='Investment')
    y           ${y}$           (long_name='Output')
    lab         ${l}$           (long_name='hours worked')
    pinf        ${\pi}$         (long_name='Inflation')
    w           ${w}$           (long_name='real wage')
    r           ${r}$           (long_name='nominal interest rate')
    rnot        ${rnot}$           (long_name='notional interest rate')
    rlag     	${rlag}$           (long_name='Inertia interest rate')
    a           ${\varepsilon_a}$       (long_name='productivity process')
    b           ${c_2*\varepsilon_t^b}$ (long_name='Scaled risk premium shock')
    g           ${\varepsilon^g}$       (long_name='Exogenous spending')
    qs          ${\varepsilon^i}$       (long_name='Investment-specific technology')
    ms          ${\varepsilon^r}$       (long_name='Monetary policy shock process') 
    spinf       ${\varepsilon^p}$       (long_name='Price markup shock process')
    sw          ${\varepsilon^w}$       (long_name='Wage markup shock process')
    kp          ${k}$           (long_name='Capital stock') 
    ;    
 
varexo ea       ${\eta^a}$      (long_name='productivity shock')
    eb          ${\eta^b}$      (long_name='risk premium shock')
    eg          ${\eta^g}$      (long_name='Spending shock')
    eqs         ${\eta^i}$      (long_name='Investment-specific technology shock')
    em          ${\eta^m}$      (long_name='Monetary policy shock')
    epinf       ${\eta^{p}}$    (long_name='Price markup shock')  
    ew          ${\eta^{w}}$    (long_name='Wage markup shock')  
        ;  
 
parameters curvw ${\varepsilon_w}$  (long_name='Curvature Kimball aggregator wages')  
    cgy         ${\rho_{ga}}$       (long_name='Feedback technology on exogenous spending')  
    curvp       ${\varepsilon_p}$   (long_name='Curvature Kimball aggregator prices')  
    constelab   ${\bar l}$          (long_name='steady state hours')  
    constepinf  ${\bar \pi}$        (long_name='steady state inflation rate')  
    constebeta  ${100(\beta^{-1}-1)}$ (long_name='time preference rate in percent')  
    cmaw        ${\mu_w}$           (long_name='coefficient on MA term wage markup')  
    cmap        ${\mu_p}$           (long_name='coefficient on MA term price markup')  
    calfa       ${\alpha}$          (long_name='capital share')  
    czcap       ${\psi}$            (long_name='capacity utilization cost')  
    csadjcost   ${\varphi}$         (long_name='investment adjustment cost')  
    ctou        ${\delta}$          (long_name='depreciation rate')  
    csigma      ${\sigma_c}$        (long_name='risk aversion')  
    chabb       ${\lambda}$         (long_name='external habit degree')  
//     ccs         ${d_4}$             (long_name='Unused parameter')  
//     cinvs       ${d_3}$             (long_name='Unused parameter')  
    cfc         ${\phi_p}$          (long_name='fixed cost share')  
    cindw       ${\iota_w}$         (long_name='Indexation to past wages')  
    cprobw      ${\xi_w}$           (long_name='Calvo parameter wages')   
    cindp       ${\iota_p}$         (long_name='Indexation to past prices')  
    cprobp      ${\xi_p}$           (long_name='Calvo parameter prices')   
    csigl       ${\sigma_l}$        (long_name='Frisch elasticity')   
    clandaw     ${\phi_w}$          (long_name='Gross markup wages')   
//     crdpi       ${r_{\Delta \pi}}$  (long_name='Unused parameter')  
    crpi        ${r_{\pi}}$         (long_name='Taylor rule inflation feedback') 
    crdy        ${r_{\Delta y}}$    (long_name='Taylor rule output growth feedback') 
    cry         ${r_{y}}$           (long_name='Taylor rule output level feedback') 
    crr         ${\rho}$            (long_name='interest rate persistence')  
    crhoa       ${\rho_a}$          (long_name='persistence productivity shock')  
    crhoas      ${d_2}$             (long_name='Unused parameter')  
    crhob       ${\rho_b}$          (long_name='persistence risk premium shock')  
    crhog       ${\rho_g}$          (long_name='persistence spending shock')  
    crhols      ${d_1}$             (long_name='Unused parameter')  
    crhoqs      ${\rho_i}$          (long_name='persistence risk premium shock')  
    crhoms      ${\rho_r}$          (long_name='persistence monetary policy shock')  
    crhopinf    ${\rho_p}$          (long_name='persistence price markup shock')  
    crhow       ${\rho_w}$          (long_name='persistence wage markup shock')  
    ctrend      ${\bar \gamma}$     (long_name='net growth rate in percent')  
    cg          ${\frac{\bar g}{\bar y}}$     (long_name='steady state exogenous spending share')
    rlb         ${rlb}$             (long_name='effective lower bound')
    zlb         ${zlb}$             (long_name='zlb dummy')
	inomobszlb	${iobszlb}$         (long_name='switch if inom is observed at the zlb')
    
cpie    //gross inflation rate
cgamma          // gross growth rate
cbeta     //discount factor 
clandap                  // fixed cost share/gross price markup
cbetabar   // growth-adjusted discount factor in Euler equation
cr  // steady state gross real interest rate
crk // R^k_{*}: steady state rental rate
cw       // steady state real wage
cikbar      // (1-k_1) in equation LOM capital, equation (8)
cik    // i_k: investment-capital ratio
clk   // labor to capital ratio
cky         // k_y: steady state output ratio
ciy                       // investment-output ratio
ccy                 // consumption-output ratio
crkky                    // z_y=R_{*}^k*k_y
cwhlc // W^{h}_{*}*L_{*}/C_{*} used in c_2 in equation (2)
cwly                   // unused parameter
conster               // steady state federal funds rate ($\bar r$)
// state init values
ainit      
binit      
ginit      
qsinit     
msinit     
epinfmainit
ewmainit   
spinfinit  
swinit     
yinit      
ygapinit      
rlaginit   
rnotinit   
kpinit     
pinfinit   
winit      
inveinit   
cinit    
// gdp level
ybar
;
// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;
zlb = 0;
        inomobszlb=1;
rlb = 0.05;//1.05;
// estimated parameters initialisation
calfa=.24;
cbeta=.9995;
csigma=1.5;
cfc=1.5;
cgy=0.51;
csadjcost= 6.0144;
chabb=    0.6361;    
cprobw=   0.8087;
csigl=    1.9423;
cprobp=   0.6;
cindw=    0.3243;
cindp=    0.47;
czcap=    0.2696;
crpi=     1.488;
crr=      0.8762;
cry=      0.0593;
crdy=     0.2347;
crhoa=    0.9977;
crhob=    0.5799;
crhog=    0.9957;
crhols=   0.9928;
crhoqs=   0.7165;
crhoas=1; 
crhoms=0;
crhopinf=0;
crhow=0;
cmap = 0;
cmaw  = 0;
constelab=0;
// Trends from data
load('usdata_pkf')
constepinf=0.7;
constebeta=0.7420;
//ctrend=0.3982;
set_param_value('ctrend',mean(dy(74:end)))
ainit = 0;      
binit = 0;      
ginit = 0;      
qsinit = 0;     
msinit = 0;     
epinfmainit = 0;
ewmainit = 0;   
spinfinit = 0;  
swinit = 0;     
yinit = 0;      
ygapinit = 0;      
rlaginit = 0;   
rnotinit = 0;   
kpinit = 0;     
pinfinit = 0;   
winit = 0;      
inveinit = 0;   
cinit = 0;     
ybar = 19;
shocks;
var ea;
stderr 0.4618;
var eb;
stderr 1.8513;
var eg;
stderr 0.6090;
var eqs;
stderr 0.6017;
var em;
stderr 0.2397;
var epinf;
stderr 0.1455;
var ew;
stderr 0.2089;
end;
monetary_shock_ar = 0;
mup_muw_shock_arma = 0;
het_shock_stddev = 1;
//@#if observe_levels && !dgp_estim
//model;
//@#else
model(linear); 
//@#endif
//deal with parameter dependencies; taken from usmodel_stst.mod 
// #cpie=1+constepinf/100;         %gross inflation rate
// #cgamma=1+ctrend/100 ;          %gross growth rate
// #cbeta=1/(1+constebeta/100);    %discount factor
// 
// #clandap=cfc;                   %fixed cost share/gross price markup
// #cbetabar=cbeta*cgamma^(-csigma);   %growth-adjusted discount factor in Euler equation
// #cr=cpie/(cbeta*cgamma^(-csigma));  %steady state gross real interest rate
// #crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou); %R^k_{*}: steady state rental rate
// #cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));      %steady state real wage
// //cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
// #cikbar=(1-(1-ctou)/cgamma);        %(1-k_1) in equation LOM capital, equation (8)
// #cik=(1-(1-ctou)/cgamma)*cgamma;    %i_k: investment-capital ratio
// #clk=((1-calfa)/calfa)*(crk/cw);    %labor to capital ratio
// #cky=cfc*(clk)^(calfa-1);           %k_y: steady state output ratio
// #ciy=cik*cky;                       %investment-output ratio
// #ccy=1-cg-cik*cky;                  %consumption-output ratio
// #crkky=crk*cky;                     %z_y=R_{*}^k*k_y
// #cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy; %W^{h}_{*}*L_{*}/C_{*} used in c_2 in equation (2)
// #cwly=1-crk*cky;                    %unused parameter
// #conster=(cr-1)*100;                %steady state federal funds rate ($\bar r$)
// flexible economy
// sticky price - wage economy
[name='FOC labor with mpl expressed as function of rk and w, SW Equation (9)']
mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
[name='FOC capacity utilization, SW Equation (7)']
zcap =  (1/(czcap/(1-czcap)))* rk ;
[name='Firm FOC capital, SW Equation (11)']
rk =  w+lab-k ;
[name='Definition capital services, SW Equation (6)']
k =  kp(-1)+zcap ;
[name='Investment Euler Equation, SW Equation (3)']
inve = (1/(1+cbetabar*cgamma))* (inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
[name='Arbitrage equation value of capital, SW Equation (4)']
pk = -r+pinf(1)-0*b 
+ (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b 
+ (crk/(crk+(1-ctou)))*rk(1) 
+ ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
[name='Consumption Euler Equation, SW Equation (2)']
c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) 
+ (1/(1+chabb/cgamma))*c(+1) 
+((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) 
- (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
[name='Aggregate Resource Constraint, SW Equation (1)']
y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
[name='Aggregate Production Function, SW Equation (5)']
y = cfc*( calfa*k+(1-calfa)*lab +a );
[name='New Keynesian Phillips Curve, SW Equation (10)']
pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
+((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
[name='Wage Phillips Curve, SW Equation (13), with (12) plugged for mu_w']
w =  (1/(1+cbetabar*cgamma))*w(-1)
+(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
+(cindw/(1+cbetabar*cgamma))*pinf(-1)
-(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
+(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
+(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
(csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
+ 1*sw ;
[name='Taylor rule, SW Equation (14)']
rnot =  crpi*(1-crr)*pinf
  //Note: y = cfc*(calfa*(zcap)+(1-calfa)*lab-stst(lab));
//    +cry*(1-crr)*(cfc*(calfa*(zcap)+(1-calfa)*lab))     
//    +crdy*((cfc*(calfa*(zcap)+(1-calfa)*lab))-(cfc*(calfa*(zcap(-1))+(1-calfa)*lab(-1))))
   +cry*(1-crr)*ygap     
   +crdy*(ygap-ygap(-1))
   +
   crr*rlag(-1) ;
[name = 'Observed interest rate',relax='zlb']
    r = rnot + ms;  
[name = 'Observed interest rate',bind='zlb']     
    r = rlb-conster + inomobszlb*ms; //-> and remove inom from obs
   
[name = 'Auxiliary notional interest rate',bind='zlb']     
rlag = rnot;
[name = 'Auxiliary notional interest rate',relax='zlb']     
rlag = r;
    [name='Output gap check']
        ygap = cfc*(calfa*(zcap)+(1-calfa)*lab);     
[name='Law of motion for productivity']              
//a = crhoa*a(-1)  + ea*exp(lea);  // for log normal prior
a = crhoa*a(-1)  + ea;
[name='Law of motion for risk premium']              
b = crhob*b(-1) + eb;
[name='Law of motion for spending process']              
g = crhog*(g(-1)) + eg + cgy*ea;
[name='Law of motion for investment specific technology shock process']              
qs = crhoqs*qs(-1) + eqs;
[name='Law of motion for monetary policy shock process']          
ms = em;
[name='Law of motion for price markup shock process']              
spinf = epinf;
epinfma=epinf;
[name='Law of motion for wage markup shock process']              
sw = ew;
ewma=ew; 
[name='Law of motion for capital, SW Equation (8) (see header notes)']              
kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;
// measurement equations
[name='Observation equation output']              
dy=y-y(-1)+ctrend;
[name='Observation equation consumption']              
dc=c-c(-1)+ctrend;
[name='Observation equation investment']              
dinve=inve-inve(-1)+ctrend;
[name='Observation equation real wage']              
dw=w-w(-1)+ctrend;
[name='Observation equation inflation']              
pinfobs = 1*(pinf) + constepinf;
[name='Observation equation interest rate']              
robs =    1*(r) + conster;
[name='Observation equation hours worked']              
labobs = lab + constelab;
end; 
    occbin_constraints;
    name 'zlb'; bind rnot+conster<=rlb; //relax rnot+conster> rlb;
    end;   
steady_state_model;
dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs = constepinf;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
//robs = ((cpie/(cbeta*cgamma^(-csigma)))-1)*100;
//robs = (cpie/cbetabar-1)*100;
//cr=cpie/cbetabar
//conster=(cr-1)*100 -> conster=(cpie/cbetabar-1)*100 = robs
labobs = constelab;
cpie=1+constepinf/100;         %gross inflation rate
cgamma=1+ctrend/100 ;          %gross growth rate
cbeta=1/(1+constebeta/100);    %discount factor
clandap=cfc;                   %fixed cost share/gross price markup
cbetabar=cbeta*cgamma^(-csigma);   %growth-adjusted discount factor in Euler equation
cr=cpie/(cbeta*cgamma^(-csigma));  %steady state gross real interest rate
crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou); %R^k_{*}: steady state rental rate
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));      %steady state real wage
// //cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
cikbar=(1-(1-ctou)/cgamma);        %(1-k_1) in equation LOM capital, equation (8)
cik=(1-(1-ctou)/cgamma)*cgamma;    %i_k: investment-capital ratio
clk=((1-calfa)/calfa)*(crk/cw);    %labor to capital ratio
cky=cfc*(clk)^(calfa-1);           %k_y: steady state output ratio
ciy=cik*cky;                       %investment-output ratio
ccy=1-cg-cik*cky;                  %consumption-output ratio
crkky=crk*cky;                     %z_y=R_{*}^k*k_y
cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy; %W^{h}_{*}*L_{*}/C_{*} used in c_2 in equation (2)
cwly=1-crk*cky;                    %unused parameter
conster=(cr-1)*100;                %steady state federal funds rate ($\bar r$)
end;

varobs dy dc dinve labobs pinfobs dw robs;

options_.occbin.smoother.max_number_of_iterations=20;
options_.occbin.smoother.first_period_occbin_update=1;
options_.occbin.likelihood.first_period_occbin_update=1;
options_.occbin.likelihood.check_ahead_periods=60;
options_.occbin.likelihood.max_check_ahead_periods=60;
options_.occbin.simul.check_ahead_periods=60;
options_.occbin.simul.max_check_ahead_periods=60;
options_.occbin.smoother.check_ahead_periods=60;
options_.occbin.smoother.max_check_ahead_periods=60;

estimated_params;
    stderr ea,,0.01,3,INV_GAMMA_PDF,0.1,2;
    stderr eb,,0.025,5,INV_GAMMA_PDF,0.1,2;
    stderr eg,,0.01,3,INV_GAMMA_PDF,0.1,2;
    stderr eqs,,0.01,3,INV_GAMMA_PDF,0.1,2;
    stderr em,,0.01,3,INV_GAMMA_PDF,0.1,2;
    stderr epinf,,0.01,3,INV_GAMMA_PDF,0.1,2;
    stderr ew,,0.01,3,INV_GAMMA_PDF,0.1,2;
    crhoa,,.01,.9999,BETA_PDF,0.5,0.20;
        crhob,,.01,.9999,BETA_PDF,0.5,0.20;     
    crhog,,.01,.9999,BETA_PDF,0.5,0.20;
    crhoqs,,.01,.9999,BETA_PDF,0.5,0.20;
    csadjcost,,2,15,NORMAL_PDF,4,1.5;
    csigma,,0.25,3,NORMAL_PDF,1.50,0.375;
    chabb,,0.001,0.99,BETA_PDF,0.7,0.1;
    cprobw,,0.3,0.95,BETA_PDF,0.5,0.1;
    csigl,,0.25,10,NORMAL_PDF,2,0.75;
    cprobp,,0.5,0.95,BETA_PDF,0.5,0.10;
    cindw,,0.01,0.99,BETA_PDF,0.5,0.15;
    cindp,,0.01,0.99,BETA_PDF,0.5,0.15;
    czcap,,0.01,1,BETA_PDF,0.5,0.15;
    cfc,,1.0,3,NORMAL_PDF,1.25,0.125;
    crpi,,1.0,3,NORMAL_PDF,1.5,0.25;
    crr,,0.5,0.975,BETA_PDF,0.75,0.10;
    cry,,0.001,0.5,NORMAL_PDF,0.125,0.05;
    crdy,,0.001,0.5,NORMAL_PDF,0.125,0.05;
    constepinf,,0.1,2.0,GAMMA_PDF,0.625,0.1;//20;
    constebeta,,0.01,2.0,GAMMA_PDF,0.25,0.1;//0.20;
    constelab,,-10.0,10.0,NORMAL_PDF,0.0,2.0;
    ctrend,,0.1,0.8,NORMAL_PDF,0.4,0.10;
    cgy,,0.01,2.0,NORMAL_PDF,0.5,0.25;
    calfa,,0.01,1.0,NORMAL_PDF,0.3,0.05;    
end;

options_.debug = 0;  
options_.occbin.likelihood.periodic_solution    = true;
options_.occbin.simul.periodic_solution         = true;
options_.occbin.smoother.periodic_solution      = true;
    
options_.noprint = false;
options_.nograph = false;
options_.plot_priors = true;


estimation(datafile='usdata_pkf', 
	//lik_init=2,//3, 
	//lik_init=3, //kalman_algo=4,//kalman_algo=3,
	nobs=210,
	first_obs=74,
	use_univariate_filters_if_singularity_is_detected=0,
	//diffuse_kalman_tol=1.e-4,
	mh_replic=0, 
	mode_compute=0, 
    mode_file = 'swpkfest_mh_mode',
	posterior_sampling_method='slice', 
	mh_nblocks=4,
	mh_drop=0.5, 
	//bayesian_irf,
	presample = 20,
	nodisplay, 
	graph_format=(eps,fig), 
    filter_covariance, smoothed_state_uncertainty,
	smoother, smoother_redux,consider_all_endogenous);
	//options_.parameter_set = 'calibration';
// steady;
// check;
// stoch_simul(order=1, nograph, nomoments);
//load swpkf/Output/swpkf_mh_mode xparam1
//options_.datafile = 'usdata_pkf';
//cli.evaluate.smoother(xparam1);
     //   identification;    
load swpkfest_mh_mode.mat fval xparam1;

if abs(cli.evaluate.posterior_kernel(xparam1)-fval)>1.e-10
    error('Computed PKF postrior kernel at the mode is not recovered!')
else
    disp('Computed PKF posterior kernel at the mode is recovered!')
    cli.evaluate.posterior_kernel(xparam1)
end

//triggers PPF
disp('PPF posterior kernel at the mode')
options_.occbin.filter.particle.status=true;
cli.evaluate.posterior_kernel(xparam1)

//triggers diagnostic plots with PPF
disp('PPF diagnostic plots')
options_.occbin.filter.particle.diagnostics.graph_periods=[172 200]; 
options_.occbin.filter.particle.diagnostics.nograph=false; 
options_.occbin.filter.particle.diagnostics.status=true;
cli.evaluate.posterior_kernel(xparam1)
options_.occbin.filter.particle.diagnostics.nograph=true; 
options_.occbin.filter.particle.diagnostics.status=false;

// posterior importance sampling using ppf starting from pkf draws
options_.occbin.posterior_importance_sampling.status = true;
options_.occbin.posterior_importance_sampling.sub_draws = 10; // set to 600 to use all draws
options_.occbin.posterior_importance_sampling.filter = 'ppf'; % can be pkf or ppf
options_.occbin.posterior_importance_sampling.orig_dname = 'swpkfest';
options_.occbin.posterior_importance_sampling.orig_fname = 'swpkfest';
options_.occbin.posterior_importance_sampling.orig_filter = 'pkf'; % can be linear or pkf
oo_ = occbin.posterior_importance_sampling(M_, estim_params_, oo_, options_, bayestopt_);