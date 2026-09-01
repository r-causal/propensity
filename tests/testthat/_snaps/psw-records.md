# the printed weights read as they are pinned

    Code
      print(continuous_records_psw())
    Output
      <psw{estimand = ate; stabilized}[6]>
      [1] 0.13806589 0.10848184 0.15284739 0.12707211 0.19968723 0.07330774
      density:   normal
      numerator: marginal
      sigma:     pooled
    Code
      print(continuous_records_psw(.density = dens_t(4)))
    Output
      <psw{estimand = ate; stabilized}[6]>
      [1] 0.13555778 0.10683962 0.15808856 0.12701939 0.19685396 0.07671633
      density:   t(df = 4)
      numerator: marginal
      sigma:     pooled
    Code
      print(continuous_records_psw(stabilize = records_numerator_model()))
    Output
      <psw{estimand = ate; stabilized}[6]>
      [1] 0.3754406 0.0993653 0.1745891 0.1223126 0.1632425 0.1326700
      density:   normal
      numerator: model
      sigma:     pooled
      stabilize: exposure ~ covariate
    Code
      print(wt_ate(records_binary_ps, records_binary_exposure))
    Output
      <psw{estimand = ate}[6]>
      [1] 1.250000 3.333333 1.666667 2.000000 2.500000 1.428571

