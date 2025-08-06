# categorical exposure validation works

    Code
      expr
    Condition <propensity_categorical_levels_error>
      Error:
      ! Categorical exposure must have more than 2 levels.
      i Found 2 levels.
      i Use binary exposure methods for 2-level exposures.

---

    Code
      expr
    Message <cliMessage>
      i Treating `.exposure` as categorical
    Condition <propensity_focal_category_error>
      Error:
      ! Focal category must be one of the exposure levels.
      i Focal category: "D"
      i Available levels: "A", "B", and "C"

# propensity score matrix validation works

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

---

    Code
      expr
    Condition <propensity_matrix_dims_error>
      Error:
      ! Number of rows in propensity score matrix must match number of observations.
      i Matrix rows: 3
      i Observations: 5

---

    Code
      expr
    Condition <propensity_matrix_dims_error>
      Error:
      ! Number of columns in propensity score matrix must match number of exposure categories.
      i Matrix columns: 2
      i Categories: 3

---

    Code
      expr
    Condition <propensity_matrix_sum_error>
      Error:
      ! Propensity score matrix rows must sum to 1.
      i Problem rows: 1

---

    Code
      expr
    Condition <propensity_range_error>
      Error:
      ! All propensity scores must be between 0 and 1.

# ATT weights work for categorical exposures

    Code
      expr
    Condition <propensity_focal_required_error>
      Error in `calculate_categorical_weights()`:
      ! Focal category must be specified for ATT with categorical exposures.

# categorical weights error on mismatched column names

    Code
      expr
    Condition <propensity_matrix_names_error>
      Error:
      ! Column names of propensity score matrix must match exposure levels.
      i Expected levels: "A", "B", and "C"
      i Found columns: "X", "Y", and "Z"

# categorical weights warn on unnamed columns

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = ate}[50]>
       [1]  2.143341  1.954254 12.093275  5.416589  3.323649 87.044712  3.460027
       [8]  3.150962  3.489142  2.506162  3.998280  3.160266  1.360694  1.460999
      [15]  2.815738  2.776858  2.019128 19.152358  4.206868  3.142638  1.867195
      [22]  3.665408  2.080233  2.499039  3.391793  2.561385  2.015574  1.729226
      [29]  3.406204  2.097834  3.542601  3.712682  2.527807  2.050151  2.334594
      [36]  1.778397  8.078721  4.868670 16.059930 39.913352  2.411579 10.288830
      [43] 15.330933  3.938484  2.763845  2.035122  2.197903  1.881399 10.170810
      [50]  2.358510

---

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = att}[50]>
       [1] 1.00000000 0.87352399 3.66175927 0.57595796 1.00000000 1.00000000
       [7] 1.26801436 1.42913785 1.29054472 1.00000000 2.32855347 0.62473452
      [13] 0.13230026 0.35851118 0.66350203 1.39859132 0.97990967 1.00000000
      [19] 1.00000000 1.33171040 0.32581785 1.00000000 0.70412427 1.00000000
      [25] 1.23587052 1.00000000 1.00000000 0.59803469 1.00000000 0.67693241
      [31] 1.53356706 1.00000000 1.05844676 0.06069529 1.00000000 0.58568250
      [37] 0.77538083 3.12025028 1.00000000 1.00000000 0.97621084 4.29259953
      [43] 1.00000000 1.77147262 1.00000000 1.00000000 1.00000000 0.83511422
      [49] 3.34284185 0.05248966

---

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = atu}[50]>
       [1]  1.143341  1.080730  8.431515  4.840631  2.323649 86.044712  2.192012
       [8]  1.721824  2.198597  1.506162  1.669727  2.535532  1.228394  1.102488
      [15]  2.152236  1.378267  1.039218 18.152358  3.206868  1.810928  1.541377
      [22]  2.665408  1.376109  1.499039  2.155922  1.561385  1.015574  1.131191
      [29]  2.406204  1.420902  2.009033  2.712682  1.469360  1.989456  1.334594
      [36]  1.192714  7.303340  1.748420 15.059930 38.913352  1.435369  5.996230
      [43] 14.330933  2.167012  1.763845  1.035122  1.197903  1.046284  6.827968
      [50]  2.306021

---

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = atm}[50]>
       [1] 0.24469633 0.08073041 1.00000000 0.57595796 0.43131206 1.00000000
       [7] 1.00000000 0.72182428 1.00000000 0.54789132 0.66972665 0.62473452
      [13] 0.13230026 0.10248791 0.66350203 0.37826672 0.03921846 1.00000000
      [19] 1.00000000 0.81092775 0.32581785 1.00000000 0.37610855 0.27411758
      [25] 1.00000000 0.27654766 0.01649960 0.13119138 1.00000000 0.42090165
      [31] 1.00000000 1.00000000 0.46936028 0.06069529 0.43403078 0.19271423
      [37] 0.77538083 0.74841982 1.00000000 1.00000000 0.43536855 1.00000000
      [43] 1.00000000 1.00000000 0.83487696 0.07564075 0.23413933 0.04628430
      [49] 1.00000000 0.05248966

---

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = ato}[50]>
       [1] 0.16130372 0.06881510 0.71040143 0.33371024 0.25994592 0.95393404
       [7] 0.38058238 0.32413847 0.38326285 0.25848300 0.34216136 0.30751062
      [13] 0.07729787 0.07381947 0.29629345 0.22942923 0.03633892 0.81911707
      [19] 0.42653956 0.33511284 0.16902337 0.39984160 0.19688874 0.18300109
      [25] 0.37393572 0.18538011 0.01597229 0.09713835 0.36972899 0.20605368
      [31] 0.37834047 0.40286757 0.24537829 0.05409382 0.22653137 0.12663934
      [37] 0.40844086 0.37641595 0.77889229 0.90342595 0.23141327 0.69778315
      [43] 0.70706441 0.41298568 0.30541395 0.06551956 0.15851480 0.04201145
      [49] 0.67993304 0.04803753

---

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.
    Output
      <psw{estimand = entropy}[50]>
       [1]  2.074514  1.630655 10.485926  4.300818  3.147657 63.485806  3.784388
       [8]  3.341375  3.813926  2.673043  3.841354  3.271719  1.023953  1.155122
      [15]  3.023828  2.734616  1.565673 16.488798  4.412698  3.386981  1.863534
      [22]  3.995362  2.138524  2.395140  3.713377  2.442551  1.481370  1.520968
      [29]  3.699735  2.182655  3.816076  4.036321  2.639062  1.652381  2.435938
      [36]  1.654485  5.470641  4.372544 13.693339 30.894723  2.508410  9.692657
      [43] 11.144542  4.205670  3.028906  1.681041  2.106362  1.481779  9.284428
      [50]  1.829680

