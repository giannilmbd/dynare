@#include "ls2003_model.inc"

estimated_params;
psi1 , gamma_pdf,1.5,0.5;
psi2 , gamma_pdf,0.25,0.125;
psi3 , gamma_pdf,0.25,0.125;
rho_R ,beta_pdf,0.5,0.2;
alpha ,beta_pdf,0.3,0.1;
rr ,gamma_pdf,2.5,1;
k , gamma_pdf,0.5,0.25;
tau ,gamma_pdf,0.5,0.2;
rho_q ,beta_pdf,0.4,0.2;
rho_A ,beta_pdf,0.5,0.2;
rho_ys ,beta_pdf,0.8,0.1;
rho_pies,beta_pdf,0.7,0.15;
stderr e_R,inv_gamma_pdf,(1.2533/3),(0.6551/10);
stderr e_q,inv_gamma_pdf,(2.5066/3),(1.3103/10);
stderr e_A,inv_gamma_pdf,(1.2533/3),(0.6551/10);
stderr e_ys,inv_gamma_pdf,(1.2533/3),(0.6551/10);
stderr e_pies,inv_gamma_pdf,(1.88/3),(0.9827/10);
end;

options_.prior_mc=5000;

% tests the prior command
prior plot table moments;
prior simulate moments(distribution) irfs(distribution);
prior optimize;
