/* Tests perfect_foresight_with_expectation_errors_{setup,solver}
   using the shocks(learnt_in=…), mshocks(learnt_in=…) and endval(learnt_in=…) syntax,
   except for the first informational period where the learnt_in keyword is omitted */

@#define bare_first_info_period = true
@#define dates = false

@#include "pfwee_learnt_in.inc"
