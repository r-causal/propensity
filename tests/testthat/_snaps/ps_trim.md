# adaptive method: ignores lower/upper, warns appropriately

    Code
      expr
    Condition <propensity_warning>
      Warning in `ps_trim()`:
      For `method = 'adaptive'`, `lower` and `upper` are ignored.

# pref method: requires exposure, fails with all 0 or all 1

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! For `method = 'pref'`, must supply `exposure`.

---

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_trim()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.treated` and `.untreated.`

---

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_trim()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.treated` and `.untreated.`

# cr method: uses min(ps_treat) / max(ps_untrt), warns if cutoffs given

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! For `method = 'cr'`, must supply `exposure`.

---

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_trim()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.treated` and `.untreated.`

---

    Code
      expr
    Condition <propensity_warning>
      Warning in `ps_trim()`:
      For `method = 'cr'`, `lower` and `upper` are ignored.
    Output
      <ps_trim; trimmed 4 of [50]>
              1         2         3         4         5         6         7         8 
      0.2854244 0.4690101 0.2562434 0.3559890 0.4981327 0.3030427 0.3129802 0.3006030 
              9        10        11        12        13        14        15        16 
      0.3274483 0.3615652 0.4551177 0.2882753 0.2683247 0.3375986 0.2689295 0.3390795 
             17        18        19        20        21        22        23        24 
      0.3034643        NA 0.3700628 0.3295265 0.4263702 0.4299815 0.4764418 0.4095624 
             25        26        27        28        29        30        31        32 
      0.4192767 0.3268432 0.4720514 0.4791835 0.2989776 0.2845567 0.3763277 0.4447085 
             33        34        35        36        37        38        39        40 
             NA 0.4541869 0.4774435 0.4308300 0.2733151        NA        NA 0.3387900 
             41        42        43        44        45        46        47        48 
      0.4837499 0.2882282 0.3442460 0.5145541 0.3141775 0.3971577 0.2820792 0.3138951 
             49        50 
      0.2939755 0.3446783 

# ps_refit() refits on keep_idx, warns if everything trimmed, etc.

    Code
      expr
    Condition <propensity_length_error>
      Error in `ps_refit()`:
      ! `.data` must have the same number of rows as observations in `trimmed_ps`.
      x `.data` has 10 rows.
      x `trimmed_ps` has 20 observations.

---

    Code
      expr
    Condition <propensity_no_data_error>
      Error in `ps_refit()`:
      ! No retained rows to refit on (all were trimmed).

# Combining two ps_trim with different parameters triggers warning

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trim.ps_trim()`:
      Converting ps_trim to numeric: different trimming parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# Combining ps_trim with double => double

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.ps_trim.double()`:
      Converting ps_trim to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# ps_trim errors when exposure is missing for methods that require it

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! For `method = 'pref'`, must supply `exposure`.

---

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! For `method = 'cr'`, must supply `exposure`.

# ps_trim warns when combining objects with different parameters

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trim.ps_trim()`:
      Converting ps_trim to numeric: different trimming parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

