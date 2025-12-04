# train() stops on NA in S or covariates; treatment variation check

    `S` contains NA.

---

    Covariates contain NA. Please impute/drop before training.

---

    Cannot train model: sample indicator column `S` has no variation.

# estimate() input validation branches

    `testing_data` must be a data frame.

---

    Column 'S' not found in testing_data.

---

    `S` contains NA in testing data.

---

    Covariates contain NA in testing data.

---

    `Y` or `Tr` contains NA among S==1 in testing data.

---

    missing value where TRUE/FALSE needed

---

    Invalid `pi`: P(S=1) must be between 0 and 1 (exclusive).

---

    Invalid `pi`: P(S=1) must be between 0 and 1 (exclusive).

---

    `pi_m` and `e_m` must be model objects (e.g. from glm) for prediction.

---

    `pi_m` and `e_m` must be model objects (e.g. from glm) for prediction.

# stratified_kfold basic properties and K>n warning path

    `S` should be a vector or factor, not a data frame or matrix.

# estimate_dml input guards hit key branches

    `data` must be a data frame.

---

    Column 'S' not found in data.

---

    `S` contains NA.

---

    Covariates contain NA. Please impute/drop before training.

---

    `Y` or `Tr` contains NA among S==1 rows.

---

    `crossfit` must be an integer >= 2.

---

    `crossfit` must be an integer >= 2.

---

    Sample indicator `S` has no variation (all 0 or all 1). Cannot perform cross-fitting.

---

    Sample indicator `S` has no variation (all 0 or all 1). Cannot perform cross-fitting.

