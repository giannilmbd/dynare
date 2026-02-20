/*
 * This file replicates the model studied in:
 * Lawrence J. Christiano, Roberto Motto and Massimo Rostagno (2007):
 * "Notes on Ramsey-Optimal Monetary Policy", Section 2
 * The paper is available at http://faculty.wcas.northwestern.edu/~lchrist/d16/d1606/ramsey.pdf
 *
 * Notes:
 *
 *  - This mod-file tests OSR with a planner_objective
 *
 *  - Due to divine coincidence, the first best policy involves fully stabilizing inflation
 *      and thereby the output gap. As a consequence, the optimal inflation feedback coefficient
 *      in a Taylor rule would be infinity. The OSR command therefore estimates it to be at the
 *      upper bound defined via osr_params_bounds.
 *
 *  - The mod-file also allows to conduct estimation under Ramsey policy by setting the
 *      Estimation_under_Ramsey switch to 1.
 *
 * This implementation was written by Johannes Pfeifer.
 *
 * If you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare
 * implementation of the model.
 */

/*
 * Copyright © 2019 Dynare Team
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <https://www.gnu.org/licenses/>.
 */

//**********Define which monetary policy setup to use ***********

@#ifndef Ramsey
    @#define Ramsey=0
@#endif

//**********Define whether to use distorted steady state***********

var C       $C$                 (long_name='Consumption')
    pi      $\pi$               (long_name='Gross inflation')
    h       $h$                 (long_name='hours worked')
    Z       $Z$                 (long_name='TFP')
    R       $R$                 (long_name='Net nominal interest rate')
    log_C   ${\ln C}$           (long_name='Log Consumption')
    log_h   ${\ln h}$           (long_name='Log hours worked')
    pi_ann  ${\pi^{ann}}$       (long_name='Annualized net inflation')
    R_ann   ${R^{ann}}$         (long_name='Annualized net nominal interest rate')
    r_real  ${r^{ann,real}}$    (long_name='Annualized net real interest rate')
    y_nat   ${y^{nat}}$         (long_name='Natural (flex price) output')
    y_gap   ${r^{gap}}$         (long_name='Output gap')
;

varexo epsilon ${\varepsilon}$     (long_name='TFP shock')
    ;

parameters beta     ${\beta}$       (long_name='discount factor')
        theta       ${\theta}$      (long_name='substitution elasticity')
        tau         ${\tau}$        (long_name='labor subsidy')
        chi         ${\chi}$        (long_name='labor disutility')
        phi         ${\phi}$        (long_name='price adjustment costs')
        rho         ${\rho}$        (long_name='TFP autocorrelation')
        @# if !defined(Ramsey) || Ramsey==0
            pi_star     ${\pi^*}$   (long_name='steady state inflation')
            alpha       ${\alpha}$  (long_name='inflation feedback Taylor rule')
        @# endif
        ;

beta=0.99;
theta=5;
phi=100;
rho=0.9;
@# if !defined(Ramsey) || Ramsey==0
    alpha=1.5;
    pi_star=1;
@# endif
tau=0;
chi=1;

model;
    [name='Euler equation']
    1/(1+R)=beta*C/(C(+1)*pi(+1));
    [name='Firm FOC']
    (tau-1/(theta-1))*(1-theta)+theta*(chi*h*C/(exp(Z))-1)=phi*(pi-1)*pi-beta*phi*(pi(+1)-1)*pi(+1);
    [name='Resource constraint']
    C*(1+phi/2*(pi-1)^2)=exp(Z)*h;
    [name='TFP process']
    Z=rho*Z(-1)+epsilon;
    @#if !defined(Ramsey) || Ramsey==0
        [name='Taylor rule']
        R=pi_star/beta-1+alpha*(pi-pi_star);
    @#endif
    [name='Definition log consumption']
    log_C=log(C);
    [name='Definition log hours worked']
    log_h=log(h);
    [name='Definition annualized inflation rate']
    pi_ann=4*log(pi);
    [name='Definition annualized nominal interest rate']
    R_ann=4*R;
    [name='Definition annualized real interest rate']
    r_real=4*log((1+R)/pi(+1));
    [name='Definition natural output']
    y_nat=exp(Z)*sqrt((theta-1)/theta*(1+tau)/chi);
    [name='output gap']
    y_gap=log_C-log(y_nat);
end;

steady_state_model;
    Z=0;
    @# if !defined(Ramsey) || Ramsey==0
        R=pi_star/beta-1; %only set this if not conditional steady state file for Ramsey
    @# endif
    pi=(R+1)*beta;
    C=sqrt((1+1/theta*((1-beta)*(pi-1)*pi-(tau-1/(theta-1))*(1-theta)))/(chi*(1+phi/2*(pi-1)^2)));
    h=C*(1+phi/2*(pi-1)^2);
    log_C=log(C);
    log_h=log(h);
    pi_ann=4*log(pi);
    R_ann=4*R;
    r_real=4*log((1+R)/pi);
    y_nat=sqrt((theta-1)/theta*(1+tau)/chi);
    y_gap=log_C-log(y_nat);
end;

@# if defined(Ramsey) && Ramsey==1
    //define initial value of instrument for Ramsey
    initval;
        R=1/beta-1;
    end;
@# endif

shocks;
    var epsilon = 0.01^2;
end;

@# if !defined(Ramsey) || Ramsey==0

    //starting value for OSR parameter
    alpha = 1.5;

    //define planner objective, which corresponds to utility function of agents
    planner_objective log(C)-chi/2*h^2;
    stoch_simul(order=2,irf=20,periods=0) pi_ann log_h R_ann log_C Z r_real;

        //define OSR parameters to be optimized
        osr_params alpha;

        //starting value for OSR parameter
        alpha = 1.5;

        //define bounds for OSR during optimization
        osr_params_bounds;
        alpha, 1, 200;
    end;

    //compute OSR and provide output
    osr(opt_algo=9,planner_discount=beta) pi_ann log_h R_ann log_C Z r_real;
    evaluate_planner_objective;

@# else
    //use Ramsey optimal policy

    //define planner objective, which corresponds to utility function of agents
    planner_objective log(C)-chi/2*h^2;

    //set up Ramsey optimal policy problem with interest rate R as the instrument,...
    // defining the discount factor in the planner objective to be the one of private agents
    ramsey_model(instruments=(R),planner_discount=beta,planner_discount_latex_name=$\beta$);

    //conduct stochastic simulations of the Ramsey problem
    stoch_simul(order=1,irf=20,periods=0) pi_ann log_h R_ann log_C Z r_real;
    evaluate_planner_objective;
@# endif

//Compare results to Ramsey results
if max(abs(-51.15717757-oo_.planner_objective_value.unconditional))>2e-6
   error('Planner objective value does not match.')
end
if max(abs(200-oo_.osr.optim_params.alpha))>1e-6
   error('Optimal parameter value does not match.')
end
