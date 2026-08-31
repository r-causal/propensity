# a marginal structural model reading a covariate is still refused

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a marginal structural model whose exposure terms read the exposure alone.
      x `e:x1` is not a term of "e" alone.
      i A term reading anything else contributes a coefficient that depends on what it reads, so there is no one effect for a row to report.
      i A term reading the exposure alone is admitted however it is written, so a curve such as `e + I(e^2)` reports one row per coefficient.
      i Read the full coefficient vector from the returned fit object for a model this surface cannot report.

---

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a marginal structural model whose exposure terms read the exposure alone.
      x `a:e` is not a term of "e" alone.
      i A term reading anything else contributes a coefficient that depends on what it reads, so there is no one effect for a row to report.
      i A term reading the exposure alone is admitted however it is written, so a curve such as `e + I(e^2)` reports one row per coefficient.
      i Read the full coefficient vector from the returned fit object for a model this surface cannot report.

# ipw() refuses dose weights built from a kernel density

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` cannot build a sandwich variance for weights built with a "kernel" density.
      x The bandwidth of a kernel estimate is chosen from the residuals of the propensity score model, so the weights are not a differentiable function of that model's parameters.
      i This route builds standard errors from the stacked system alone. Build the dose weights from a model and a density that system can differentiate, or bootstrap the whole joint fit yourself: resample the rows, refit both treatment models, rebuild the weights with `wt_joint()`, and refit the outcome model on each resample.

# ipw() refuses a dose psi it cannot write, under the component's name

    Code
      expr
    Condition <propensity_ipw_robust_psi_error>
      Error in `ipw()`:
      ! `ipw()` cannot write the equation this <rlm/lm> propensity score model of a continuous exposure is the root of.
      x `e` was fit with a psi function this path cannot recognize.
      i Refit `e` with `MASS::psi.huber()`, `MASS::psi.bisquare()`, or `MASS::psi.hampel()`, whose constants `ipw()` reads off the fit.
      i This route builds standard errors from the stacked system alone. Build the dose weights from a model and a density that system can differentiate, or bootstrap the whole joint fit yourself: resample the rows, refit both treatment models, rebuild the weights with `wt_joint()`, and refit the outcome model on each resample.

# the weights mismatch names the ratio the dose records

    Code
      expr
    Condition <propensity_ipw_weights_mismatch_error>
      Error in `ipw()`:
      ! The "ate" weights recomputed from `wt_mod` differ from the weights supplied to `outcome_mod` (compared at relative tolerance 1e-6).
      i The estimand the weights were built for may differ from the one `ipw()` resolved.
      i `ipw()` rebuilt these weights as a "t(df = 4)" density with a "marginal" numerator.
      i Weights built with an observation-level `.sigma`, such as `influence(model)$sigma`, are one cause: `ipw()` models the conditional density with a single pooled residual root mean square, which is what `wt_ate()` uses when no `.sigma` is given.
      i Weights trimmed, truncated, or normalized after `wt_mod` was fit differ from the ones rebuilt here, which come from that model alone.
      i `.data` values that differ from the data the models were fit to move the recomputed weights on their own and leave the supplied weights exactly right.
      i Refit `outcome_mod` with weights from the two treatment models `wt_mod` holds, and this estimand, if the weights are the cause.

# ipw() still refuses a term reading a treatment and a covariate

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a joint marginal structural model whose treatment terms read the treatments alone.
      x `a:x1`, `e:x1`, and `a:e:x1` in `outcome_mod` are not terms of "a" and "e" alone.
      i A coefficient of a term reading a covariate is a change in an effect per unit of that covariate, so no row could name the effect it stands for.
      i A covariate entering on its own is admitted and contributes no row, so adjust `outcome_mod` for it that way, or cross two discrete treatments to report a surface that standardizes over the covariates.

---

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a joint marginal structural model whose treatment terms read the treatments alone.
      x `e:x1` in `outcome_mod` is not a term of "a" and "e" alone.
      i A coefficient of a term reading a covariate is a change in an effect per unit of that covariate, so no row could name the effect it stands for.
      i A covariate entering on its own is admitted and contributes no row, so adjust `outcome_mod` for it that way, or cross two discrete treatments to report a surface that standardizes over the covariates.

# ipw() refuses a treatment column coded some other way

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `a` and `a:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "a" enters as 0 for "no" and 1 for "yes", "e" enters as itself, and their interaction is the product of the two.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i A model with no intercept, written `- 1` or `+ 0`, expands a factor treatment to an indicator for every level, so its first column is the reference-level indicator rather than the 0/1 indicator the rows describe. Keep the intercept, or code "a" as a 0/1 numeric.
      i Refit `outcome_mod` with "a" as a 0/1 numeric, or as an unordered factor under treatment contrasts.

---

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `a` and `a:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "a" enters as 0 for "no" and 1 for "yes", "e" enters as itself, and their interaction is the product of the two.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i A model with no intercept, written `- 1` or `+ 0`, expands a factor treatment to an indicator for every level, so its first column is the reference-level indicator rather than the 0/1 indicator the rows describe. Keep the intercept, or code "a" as a 0/1 numeric.
      i Refit `outcome_mod` with "a" as a 0/1 numeric, or as an unordered factor under treatment contrasts.

# the weights mismatch names a fixed stabilization score

    Code
      expr
    Condition <propensity_ipw_weights_mismatch_error>
      Error in `ipw()`:
      ! The "ate" weights recomputed from `wt_mod` differ from the weights supplied to `outcome_mod` (compared at relative tolerance 1e-6).
      i The estimand the weights were built for may differ from the one `ipw()` resolved.
      i `ipw()` rebuilt these weights as a "normal" density with a "marginal" numerator.
      i Weights built with an observation-level `.sigma`, such as `influence(model)$sigma`, are one cause: `ipw()` models the conditional density with a single pooled residual root mean square, which is what `wt_ate()` uses when no `.sigma` is given.
      i A component built with a fixed `stabilization_score` is one cause: the product records that a numerator was a score without recording the vector it was, so `ipw()` rebuilt the numerator of `e` from the exposure's own marginal distribution instead.
      i Weights trimmed, truncated, or normalized after `wt_mod` was fit differ from the ones rebuilt here, which come from that model alone.
      i `.data` values that differ from the data the models were fit to move the recomputed weights on their own and leave the supplied weights exactly right.
      i Refit `outcome_mod` with weights from the two treatment models `wt_mod` holds, and this estimand, if the weights are the cause.

# a bare-term model with no intercept is refused, not errored

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `a` in `outcome_mod` contributes a column coded some other way.
      i The reported rows name the coefficients of a model in which "a" enters as 0 for "no" and 1 for "yes", "e" enters as itself, and their interaction is the product of the two.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i A model with no intercept, written `- 1` or `+ 0`, expands a factor treatment to an indicator for every level, so its first column is the reference-level indicator rather than the 0/1 indicator the rows describe. Keep the intercept, or code "a" as a 0/1 numeric.
      i Refit `outcome_mod` with "a" as a 0/1 numeric, or as an unordered factor under treatment contrasts.

# ipw() refuses .by on a joint continuous fit

    Code
      expr
    Condition <propensity_ipw_joint_by_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a joint exposure.
      x A joint exposure already reports an interaction between two treatments, and reporting it within the levels of a modifier is a three-way question this surface does not answer.
      i Drop `.by` to report the joint surface, or weight the crossing of the two treatments as one plain categorical exposure to report each cell against the reference cell within each subgroup.

# ipw() refuses linearization on a joint continuous fit

    Code
      expr
    Condition <propensity_ipw_joint_models_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support "linearization" standard errors for a joint treatment model.
      x The cell means and every contrast built from them are parameters of the stacked estimating equations, and the linearization path solves no such system.
      i Use `se_method = "mestimation"` for a joint treatment model.

# the joint dose weights refusal comes before the estimates do

    Code
      expr
    Condition <propensity_ipw_ps_weights_error>
      Error in `ipw()`:
      ! `ipw()` does not support a treatment model fit with case weights.
      x `e` was fit with non-unit `weights`, so its coefficients are not the root of the unweighted score stacked for it.
      i Refit `e` without `weights`.

