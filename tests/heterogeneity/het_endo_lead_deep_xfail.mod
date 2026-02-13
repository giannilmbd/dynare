// EXPECTEDFAIL Test: Heterogeneous endogenous with lead > 1 is not supported
//
// This test verifies that the preprocessor correctly rejects heterogeneous
// endogenous variables with leads greater than +1 in model(heterogeneity=...) blocks.
//
// Example: c(+2) where c is var(heterogeneity=households)
//
// Expected result: Preprocessor error with message about lead +2 not supported

heterogeneity_dimension households;

var(heterogeneity=households) c k;
varexo(heterogeneity=households) e;
var r;
parameters beta;

beta = 0.99;

model(heterogeneity=households);
    k(-1) - c(+2);  // ERROR: lead +2 not supported
    k - e;
end;

model;
    r - 0.04;
end;
