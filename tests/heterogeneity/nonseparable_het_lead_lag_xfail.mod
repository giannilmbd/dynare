// EXPECTEDFAIL Test: Non-separable het lead/lag expressions are not supported
//
// This test verifies that the preprocessor correctly rejects non-separable
// expressions that combine leads with lagged states in model(heterogeneity=...) blocks.
//
// Example: log(k(-1) + c(+1)) - the log() makes it non-separable
//
// Expected result: Preprocessor error with message about non-separable expression

heterogeneity_dimension households;

var(heterogeneity=households) c k;
varexo(heterogeneity=households) e;
var r;
parameters beta;

beta = 0.99;

model(heterogeneity=households);
    // Non-separable term combining lead and lagged state
    log(k(-1) + c(+1)) - c;
    k - e;
end;

model;
    r - 0.04;
end;
