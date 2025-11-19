library(devtools)
library(tidyverse)
#library(pbapply)

remove.packages("ROOT")
.rs.restartR()
devtools::document()
devtools::check()
res <- devtools::check()
res$errors
res$warnings
res$notes

devtools::install()

library(ROOT)
devtools::build_manual()

pkgload::load_all()
styler::style_pkg()
devtools::install()






#########################################################################################
#########################################################################################
#                                                                                       #
#                                                                                       #
#                           Example of a synthetic dataset                              #
#                                                                                       #
#                                                                                       #
#########################################################################################
#########################################################################################

###Simulate a synthetic target and RCT dataset
make_root_datasets <- function(
    n_rct      = 1500,
    n_target   = 4000,
    p          = 6,      # number of continuous covariates
    rho        = 0.35,   # common correlation between covariates
    mean_shift = 0.45,   # shift of target means vs. RCT
    seed       = 123
) {
  stopifnot(p >= 1)
  set.seed(seed)

  # (1) helpers
  Sigma <- matrix(rho, p, p); diag(Sigma) <- 1

  # fast log N(x | mu, Sigma) without extra deps
  logdmvnorm <- function(X, mu, Sigma) {
    X <- as.matrix(X)
    p <- ncol(X)
    R <- chol(Sigma)
    z <- t(backsolve(R, t(X) - mu, transpose = TRUE))
    quad <- rowSums(z * z)
    const <- -0.5 * p * log(2 * pi) - sum(log(diag(R)))
    const - 0.5 * quad
  }

  # (2) generate covariates (ALL continuous)
  mu_rct <- rep(0, p)
  mu_tgt <- mu_rct + mean_shift

  X_rct <- MASS::mvrnorm(n_rct, mu = mu_rct, Sigma = Sigma)
  X_tgt <- MASS::mvrnorm(n_target, mu = mu_tgt, Sigma = Sigma)
  colnames(X_rct) <- colnames(X_tgt) <- paste0("X", seq_len(p))

  # (3) treatment assignment in the RCT
  beta_t <- seq(0.8, 0.2, length.out = p)        # weights
  lin_t  <- -0.2 + as.vector(X_rct %*% beta_t)
  pr_t   <- 1 / (1 + exp(-lin_t))
  Tr     <- rbinom(n_rct, 1, pr_t)

  # (4) outcome in the RCT (heterogeneous treatment effect)
  beta_y <- seq(0.6, 0.1, length.out = p)
  theta  <- seq(0.5, 0.1, length.out = p) * 0.5
  tau0   <- 0.4
  b0     <- -0.3
  lin_y  <- b0 + as.vector(X_rct %*% beta_y) + Tr * (tau0 + as.vector(X_rct %*% theta))
  pr_y   <- 1 / (1 + exp(-lin_y))
  Y      <- rbinom(n_rct, 1, pr_y)

  # (5) rarity score lX for target rows (under RCT X distribution)
  logp_target <- logdmvnorm(X_tgt, mu = mu_rct, Sigma = Sigma)
  lX          <- as.numeric(scale(logp_target))  # mean 0, sd 1

  # (6) assemble outputs
  DataRCT    <- data.frame(X_rct, Tr = Tr, Y = Y, check.names = FALSE)
  DataTarget <- data.frame(X_tgt, lX = lX, check.names = FALSE)
  cov_names  <- colnames(X_rct)

  # ROOT-ready combined frame: S=1 for RCT, S=0 for Target; Y/Tr NA on Target
  D_root <- rbind(
    data.frame(X_rct, Yobs = Y, Tr = Tr, S = 1L, check.names = FALSE),
    data.frame(X_tgt, Yobs = NA_real_, Tr = NA_integer_, S = 0L, check.names = FALSE)
  )
  rownames(D_root) <- NULL

  list(
    DataRCT              = DataRCT,
    DataTarget           = DataTarget,
    covariate_DataRCT    = cov_names,
    covariate_DataTarget = cov_names,
    treatment_DataRCT    = "Tr",
    outcome_DataRCT      = "Y",
    D_root               = D_root
  )
}


###Generate data
gen <- make_root_datasets(
  n_rct = 1500, n_target = 4000, p = 6, rho = 0.35, mean_shift = 0.45, seed = 123
)

### Run ROOT directly using an arbitrary data
outputROOT <- ROOT(
  data     = gen$D_root,
  outcome  = "Yobs",
  treatment= "Tr",
  sample   = "S",
  leaf_proba = 0.25,
  seed       = 3,
  num_trees  = 50,
  vote_threshold = 2/3,
  explore_proba  = 0.05,
  feature_est    = "Ridge",
  top_k_trees    = FALSE,
  verbose  = TRUE,
  cutoff         = "baseline"
)
summary(outputROOT)
plot(outputROOT)


### Run characterizing_underrep function
### Output underrepresented population characterization tree
outputUnderrep <- characterizing_underrep(
  DataRCT               = gen$DataRCT,
  covariateColName_RCT  = gen$covariate_DataRCT,
  trtColName_RCT        = gen$treatment_DataRCT,
  outcomeColName_RCT    = gen$outcome_DataRCT,
  DataTarget            = gen$DataTarget,
  covariateColName_TargetData  = gen$covariate_DataTarget,
  num_trees             = 50,
  seed                  = 3,
  verbose               = TRUE
 )
summary(outputUnderrep)
plot(outputUnderrep)



# --------------- Custom loss functions ----------------

# Factory returning an objective_fn(D) -> scalar
objective_variance_plus_drop <- function(lambda = 0.05) {
  force(lambda)
  function(D) {
    wsum <- sum(D$w)
    if (wsum <= 0) return(Inf)                           # no kept units → invalid
    # default variance proxy
    var_term <- sqrt(sum(D$vsq * D$w) / (wsum^2))
    # penalize dropping too many (1 - mean(w) is drop rate)
    pen <- lambda * (1 - mean(D$w))
    var_term + pen
  }
}
obj <- objective_variance_plus_drop(lambda = 0.05)

### Run ROOT directly using an arbitrary data with a user-specified loss function
outputROOT2 <- ROOT(
  data     = gen$D_root,
  outcome  = "Yobs",
  treatment= "Tr",
  sample   = "S",
  leaf_proba = 0.25,
  seed       = 3,
  num_trees  = 50,
  vote_threshold = 2/3,
  explore_proba  = 0.05,
  feature_est    = "Ridge",
  top_k_trees    = FALSE,
  verbose  = TRUE,
  cutoff         = "baseline",
  global_objective_fn = obj,
)
summary(outputROOT2)
plot(outputROOT2)


### Run characterizing_underrep function with a user-specified loss function
outputUnderrep2 <- characterizing_underrep(
  DataRCT               = gen$DataRCT,
  covariateColName_RCT  = gen$covariate_DataRCT,
  trtColName_RCT        = gen$treatment_DataRCT,
  outcomeColName_RCT    = gen$outcome_DataRCT,
  DataTarget            = gen$DataTarget,
  covariateColName_TargetData  = gen$covariate_DataTarget,
  num_trees             = 50,
  seed                  = 3,
  verbose               = TRUE,
  global_objective_fn = obj
)
summary(outputUnderrep2)
plot(outputUnderrep2)




#########################################################################################
#########################################################################################
#                                                                                       #
#                                                                                       #
#                  Example of a synthetic trial data for diabetes                       #
#                                                                                       #
#                                                                                       #
#########################################################################################
#########################################################################################
# ---- Illustrative example: T2D blood sugar RCT vs target ----
# Assumptions (from your prompt):
# - Covariates: race (Black vs non-Black), birth-sex (male vs female),
#               diet management (yes/no), age (<45 vs >=45)
# - Target prevalences: Black female = 20%; diet yes = 10%; age>=45 = 15%
# - Trial imbalance: Black female only 5% in the TRIAL sample
# - Outcome model (ε ~ N(0,1)):
#   Y = 1 + I(Treatment==1) * (1 + 1.6*DietYes + 1.3*Age45)
#         - 0.5*Black + 0.3*Male + DietYes + Age45 + ε

simulate_diabetes_illustrative <- function(n_target = 10000,
                                           n_trial  = 2000,
                                           seed = 1,
                                           trial_fraction = NULL) {
  set.seed(seed)

  # --- 1) Target covariates (four binaries) ---
  p_BF  <- 0.20; p_BM <- 0.10; p_NBF <- 0.35; p_NBM <- 0.40
  joint_levels <- c("Black_Female","Black_Male","NonBlack_Female","NonBlack_Male")
  joint_draw   <- sample(joint_levels, n_target, replace = TRUE,
                         prob = c(p_BF, p_BM, p_NBF, p_NBM))
  Race_Black <- as.integer(grepl("^Black", joint_draw))
  Sex_Male   <- as.integer(grepl("Male$", joint_draw))
  DietYes    <- rbinom(n_target, 1, 0.10)
  Age45      <- rbinom(n_target, 1, 0.15)

  target <- data.frame(Race_Black, Sex_Male, DietYes, Age45)

  # --- 2) Build TRIAL with 5% Black females ---
  if (is.null(trial_fraction)) trial_fraction <- n_trial / n_target
  n_trial <- round(trial_fraction * n_target)

  is_BF <- target$Race_Black == 1 & target$Sex_Male == 0
  n_BF_trial    <- round(0.05 * n_trial)
  n_notBF_trial <- n_trial - n_BF_trial

  idx_BF     <- which(is_BF)
  idx_notBF  <- which(!is_BF)
  trial_idx  <- c(sample(idx_BF, n_BF_trial, replace = FALSE),
                  sample(idx_notBF, n_notBF_trial, replace = FALSE))

  S <- integer(n_target); S[trial_idx] <- 1L  # S=1 trial, S=0 target
  target$S <- S

  # --- 3) Randomize treatment IN TRIAL ONLY: Tr ∈ {0,1} with 1=Treatment 1 ---
  Tr <- rep(NA_integer_, n_target)
  Tr[S == 1L] <- rbinom(sum(S == 1L), 1, 0.5)  # 1 => Treatment 1, 0 => Treatment 2

  # --- 4) Outcome DGP (ε ~ N(0,1)) ---
  # Y = 1 + Tr*(1 + 1.6*DietYes + 1.3*Age45)
  #       - 0.5*Race_Black + 0.3*Sex_Male + DietYes + Age45 + ε
  eps       <- rnorm(n_target, 0, 1)
  base_term <- 1 - 0.5*target$Race_Black + 0.3*target$Sex_Male + target$DietYes + target$Age45
  treat_bump <- 1 + 1.6*target$DietYes + 1.3*target$Age45

  # Potential outcomes (for clarity)
  Y1_pot <- base_term + treat_bump + eps      # if Tr=1 (Treatment 1)
  Y0_pot <- base_term + 0          + eps      # if Tr=0 (Treatment 2)

  # Observed Y only for trial rows where Tr is assigned
  Yobs <- rep(NA_real_, n_target)
  Yobs[S == 1L] <- ifelse(Tr[S == 1L] == 1L, Y1_pot[S == 1L], Y0_pot[S == 1L])

  dat_all <- cbind(
    target[, c("Race_Black","Sex_Male","DietYes","Age45","S")],
    Tr = Tr,                # 1 = Treatment 1, 0 = Treatment 2
    Y  = Yobs,
    Y1_pot = Y1_pot,        # optional: for diagnostics
    Y0_pot = Y0_pot
  )

  # --- 5) π(X) = P(S=1 | X) for Figure 4 ---
  sel_fit <- glm(S ~ Race_Black + Sex_Male + DietYes + Age45 + Race_Black:Sex_Male,
                 data = dat_all, family = binomial())
  dat_all$pi_hat <- as.numeric(predict(sel_fit, type = "response"))

  # --- 6) Quick checks ---
  prop_target_BF <- mean(dat_all$Race_Black==1 & dat_all$Sex_Male==0 & dat_all$S==0)
  prop_trial_BF  <- mean(dat_all$Race_Black==1 & dat_all$Sex_Male==0 & dat_all$S==1)

  list(
    data = dat_all,
    selection_model = sel_fit,
    checks = c(
      target_BF = prop_target_BF,  # ~0.20
      trial_BF  = prop_trial_BF,   # ~0.05
      diet_overall  = mean(dat_all$DietYes==1),
      age45_overall = mean(dat_all$Age45==1)
    )
  )
}

# ---- Example run ----
data <- simulate_diabetes_illustrative(n_target = 10000, n_trial = 2000, seed = 12)
head(data)
newdata = data$data[,c(1:7)]
trial  <- subset(newdata, S == 1)
target <- subset(newdata, S == 0)


### Run ROOT directly using DGP of a trial on Type 2 diabetes
outputROOT.DGP <- ROOT(
  data     = newdata,
  outcome  = "Y",
  treatment= "Tr",
  sample   = "S",
  leaf_proba = 0.25,
  seed       = 3,
  num_trees  = 50,
  vote_threshold = 2/3,
  explore_proba  = 0.05,
  feature_est    = "Ridge",
  top_k_trees    = FALSE,
  verbose  = TRUE,
  cutoff         = "baseline"
)
summary(outputROOT.DGP)
plot(outputROOT.DGP)

### Run characterizing_underrep function
outputUnderrep.DGP <- characterizing_underrep(
  DataRCT               = trial,
  covariateColName_RCT  = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
  trtColName_RCT        = "Tr",
  outcomeColName_RCT    = "Y",
  DataTarget            = target,
  covariateColName_TargetData  = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
  num_trees             = 50,
  seed                  = 3,
  verbose               = TRUE
 )
summary(outputUnderrep.DGP)
plot(outputUnderrep.DGP)



