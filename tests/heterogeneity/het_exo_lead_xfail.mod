// EXPECTEDFAIL Test: Lead on heterogeneous exogenous variable is not supported
//
// This test verifies that the preprocessor correctly rejects leads
// on heterogeneous exogenous variables in model(heterogeneity=...) blocks.
//
// Example: e(+1) where e is varexo(heterogeneity=households)
//
// Expected result: Preprocessor error with message about het exo lead

heterogeneity_dimension households;

var(heterogeneity=households) c k;
varexo(heterogeneity=households) e;
var r;
parameters beta;

beta = 0.99;

model(heterogeneity=households);
    c - e(+1);  // ERROR: het exo lead
    k - 1;
end;

model;
    r - 0.04;
end;
