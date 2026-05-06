# Changelog

## propensity (development version)

- Fixed `broom::tidy(glm_fit, conf.int = TRUE)` failing on GLMs weighted
  by `psw` vectors. `confint.glm()` builds profile-likelihood intervals
  via `profile.glm()`, which refits through
  [`glm.fit()`](https://rdrr.io/r/stats/glm.html); the refit indexes
  `weights[good]` with a matrix subscript, which `[.vctrs_vctr`
  rejected. Added a `[.psw` method that falls back to base R linear
  indexing for matrix/array subscripts and delegates everything else to
  `[.vctrs_vctr`.

- Comparison operators on `psw` (`==`, `!=`, `<`, `>`, `<=`, `>=`) now
  short-circuit `vec_equal()` / `vec_compare()` and return a logical
  vector silently. Previously each comparison fired a
  `propensity_class_downgrade` warning via `vec_ptype2.psw.double()`,
  producing 100+ warnings during a single `tidy(glm, conf.int = TRUE)`
  call. Combine and cast paths still warn.

- Added a `NEWS.md` file to track changes to the package.
