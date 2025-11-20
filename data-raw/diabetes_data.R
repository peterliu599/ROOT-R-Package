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
diabetes_data <- simulate_diabetes_illustrative(n_target = 10000, n_trial = 2000, seed = 12)
diabetes_data <- diabetes_data$data[, c(1:7)]
# usethis::use_data(diabetes_data)
