/* Tests perfect_foresight_with_expectation_errors_{setup,solver}
   using the shocks(learnt_in=…), mshocks(learnt_in=…) and endval(learnt_in=…) syntax
   with dates (instead of integer indices) for periods */

@#define bare_first_info_period = false
@#define dates = true

@#include "pfwee_learnt_in.inc"
