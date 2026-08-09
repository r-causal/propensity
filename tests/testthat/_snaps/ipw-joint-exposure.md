# a joint exposure refuses a focal estimand

    Code
      expr
    Condition <propensity_ipw_joint_estimand_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint exposure for the "ate" estimand only.
      x The weights were built for the "att" estimand.
      i Every cell mean on the joint surface standardizes to one population, and a tilted estimand standardizes each of them to a population the simple effects and the interaction are not defined over.
      i Refit the weights for "ate", or drop the declaration with `factor(x)` to report each cell against the reference cell under the estimand you have.

# a joint exposure refuses .by

    Code
      expr
    Condition <propensity_ipw_joint_by_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a joint exposure.
      x A joint exposure already reports an interaction between two treatments, and reporting it within the levels of a modifier is a three-way question this surface does not answer.
      i Drop `.by` to report the joint surface, or drop the declaration with `factor(x)` to report each cell against the reference cell within each subgroup.

