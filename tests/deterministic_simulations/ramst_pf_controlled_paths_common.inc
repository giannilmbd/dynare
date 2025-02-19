var c k y;
varexo x e_y;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
y = 0.98*y(-1)+e_y;
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
y = 0;
end;

steady;

@#if shocks
perfect_foresight_controlled_paths;
  exogenize c;
  periods 2 4:5;
  values 1.6 1.7;
  endogenize x;

  exogenize k;
  periods 7:9;
  values 13;
  endogenize x;

  exogenize y;
  periods 2003Y:2005Y;
  values 0.5;
  endogenize e_y;
end;

shocks;
  var x;
  periods 1;
  values 1.2;
end;
@#endif
