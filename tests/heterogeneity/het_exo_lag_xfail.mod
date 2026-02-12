// EXPECTEDFAIL Test: Lagged heterogeneous exogenous variables are not supported
//
// This test verifies that the preprocessor correctly rejects lagged
// heterogeneous exogenous variables in model(heterogeneity=...) blocks.
//
// Example: e(-1) where e is varexo(heterogeneity=households)
//
// Expected result: Preprocessor error with message about lagged het exo

heterogeneity_dimension households;

var(heterogeneity=households) c k;
varexo(heterogeneity=households) e;
var r;
parameters beta;

beta = 0.99;

model(heterogeneity=households);
    c - e(-1);  // ERROR: lagged het exo
    k - 1;
end;

model;
    r - 0.04;
end;
