# the conflict warning names the record it dropped

    Code
      out <- c(normal, heavier)
    Condition <propensity_metadata_conflict_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Dropping the joint_wt_meta attribute from the result.
      i The two sets of weights record it differently, so neither value describes the result.
      i Every attribute the two agree on is carried through.

# wt_joint() refuses a component that records no exposure type

    Code
      expr
    Condition <propensity_wt_joint_exposure_type_error>
      Error in `wt_joint()`:
      ! `wt_joint()` requires each component to record the exposure type it weights.
      x `w_e` records none.
      i A weight built by hand, or by a version of the package that did not record it, carries an estimand and a stabilization status and nothing about the exposure, so the requirement that a continuous component be stabilized could not be applied to it.
      i Rebuild it with a weight function such as `wt_ate()`, or name both types in `exposure_type`, for example `exposure_type = c("binary", "continuous")`.

# wt_joint() refuses a component that does not target the ate

    Code
      expr
    Condition <propensity_wt_joint_estimand_error>
      Error in `wt_joint()`:
      ! `wt_joint()` requires both components to target the "ate" estimand.
      x `w_a` targets "att" and `w_e` targets "ate".
      i The product targets the joint "ate", and a tilted component targets a population the other one does not, so the product targets no single population.
      i Rebuild both components with `wt_ate()`.

# wt_joint() refuses an unstabilized continuous component

    Code
      expr
    Condition <propensity_wt_joint_stabilize_error>
      Error in `wt_joint()`:
      ! `wt_joint()` requires a continuous component to be stabilized.
      x `w_e` weights a continuous exposure and is not stabilized.
      i The unstabilized density ratio has a heavy right tail on its own, and multiplying it by a second weight inherits that tail, leaving the product with no usable variance.
      i Rebuild it with `stabilize = TRUE`. A binary or categorical component needs no stabilization.

# wt_joint() refuses components of different lengths

    Code
      expr
    Condition <propensity_wt_joint_length_error>
      Error in `wt_joint()`:
      ! `wt_joint()` requires its two components to be the same length.
      x `w_a` is length 400; `w_e` is length 10.
      i The product is elementwise, so each unit needs one weight from each component.
      i Build both components from the same observations.

# wt_joint() refuses a component that is not a propensity score weight

    Code
      expr
    Condition <propensity_wt_joint_class_error>
      Error in `wt_joint()`:
      ! `wt_joint()` requires both components to be <psw> vectors.
      x `w_a` is not.
      i A bare numeric records neither an estimand nor a stabilization status, so the checks the product depends on cannot be applied to it and the result could record nothing about where it came from.
      i Build each component with a weight function such as `wt_ate()`.

# wt_joint() refuses an exposure_type that does not name both components

    Code
      expr
    Condition <propensity_wt_joint_exposure_type_error>
      Error in `wt_joint()`:
      ! `exposure_type` must name the exposure type of each component.
      i Supply one per component, in the order the components were given, for example `exposure_type = c("binary", "continuous")`, or leave `exposure_type` unset to read each component's own record.
      Caused by error in `wt_joint()`:
      ! `exposure_type` must be one of "binary", "categorical", or "continuous", not "ordinal".

# the unsupported-model refusal names the family of an additive fit

    Code
      expr
    Condition <propensity_wt_joint_models_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` must be given models it can read a treatment density from.
      x The model named `k` is <gam> fit with `poisson()`.
      i Supported models: a binomial `glm()` for a binary treatment, a `nnet::multinom()` for a categorical one, and an `lm()`, a gaussian `glm()`, a gaussian `mgcv::gam()`, or a `MASS::rlm()` for a continuous one.

# joint_wt_models() requires a discrete second model to condition on the first treatment

    Code
      expr
    Condition <propensity_wt_joint_factorization_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` requires the second model to condition on the first treatment.
      x `e` does not read "a" on its right-hand side.
      i A joint weight factorizes as f(a | L) f(e | a, L). The product of two marginal models, f(a | L) f(e | L), is a different quantity, and it is not the joint weight wherever "e" depends on "a".
      i Nothing downstream can tell the two apart: the product is an ordinary vector of positive numbers either way.
      i Add "a" to the formula of `e`, and model that dependence flexibly rather than as a single additive term.

# joint_wt_models() refuses models that condition on each other

    Code
      expr
    Condition <propensity_wt_joint_circular_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` requires a sequential factorization, and each of these models conditions on the other's treatment.
      x `a` reads "e", and `e` reads "a".
      i A sequential factorization has a first factor that is marginal in the second treatment, so such a pair are not the two factors of any factorization, in either order.
      i Refit the first treatment's model without "e".

# the factorization refusal offers the swapped order when the pair is one reversed

    Code
      expr
    Condition <propensity_wt_joint_factorization_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` requires the second model to condition on the first treatment.
      x `e` does not read "a" on its right-hand side.
      i A joint weight factorizes as f(a | L) f(e | a, L). The product of two marginal models, f(a | L) f(e | L), is a different quantity, and it is not the joint weight wherever "e" depends on "a".
      i Nothing downstream can tell the two apart: the product is an ordinary vector of positive numbers either way.
      i Add "a" to the formula of `e`, and model that dependence flexibly rather than as a single additive term.
      i `a` already reads "e", so the pair in the other order, `joint_wt_models(e = ..., a = ...)`, is the factorization f(e | L) f(a | e, L). Supply them that way if that is the order the treatments are assigned in.

# joint_wt_models() refuses arguments that do not name exactly two treatments

    Code
      expr
    Condition <propensity_wt_joint_models_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` must be given exactly two models, each named for the treatment it fits.
      x 1 of them is unnamed.
      i The names are the treatment names, and their order is the factorization's order: the first treatment is the one the second model conditions on.

# joint_wt_models() refuses a model whose response is not its treatment

    Code
      expr
    Condition <propensity_wt_joint_response_error>
      Error in `joint_wt_models()`:
      ! `joint_wt_models()` requires each model's response to be the treatment it is named for.
      x The model named `a` has response `e`.
      i The factorization check looks for the first treatment by name in the second model's formula, so a model recorded under the wrong name would be checked for the wrong variable.
      i Rename the argument to the treatment the model fits, or supply the model that fits "a".

