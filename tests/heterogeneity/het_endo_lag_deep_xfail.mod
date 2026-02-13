// EXPECTEDFAIL Test: Heterogeneous endogenous with lag < -1 is not supported
//
// This test verifies that the preprocessor correctly rejects heterogeneous
// endogenous variables with lags deeper than -1 in model(heterogeneity=...) blocks.
//
// Example: k(-2) where k is var(heterogeneity=households)
//
// Expected result: Preprocessor error with message about lag -2 not supported

heterogeneity_dimension households;

var(heterogeneity=households) c k;
varexo(heterogeneity=households) e;
var r;
parameters beta;

beta = 0.99;

model(heterogeneity=households);
    c - k(-2);  // ERROR: lag -2 not supported
    k - e;
end;

model;
    r - 0.04;
end;
