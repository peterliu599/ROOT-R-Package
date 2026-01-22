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

### A simple generalizability example
### Using the diabetes_data.rda dataset in the data folder

ROOT.output <- ROOT(diabetes_data, generalizability_path = TRUE, seed=123)
summary(ROOT.output)
plot(ROOT.output)

char.output <- characterizing_underrep(diabetes_data,generalizability_path = TRUE, seed = 123)
summary(char.output)
plot(char.output)

### A simple optimization example
### ROOT minimizes variance globally
set.seed(123)
n <- 1000

X1 <- sample(0:1, n, replace = TRUE)
X2 <- sample(0:1, n, replace = TRUE)

# XOR pattern: low variance on diagonal (X1==X2), high variance off-diagonal
Y <- ifelse(X1 == X2,
            rnorm(n, mean = 5, sd = 1),
            rnorm(n, mean = 5, sd = 10))

data <- data.frame(X1 = X1, X2 = X2, v = Y)

variance_objective <- function(D) {
  w <- D$w
  if (sum(w, na.rm = TRUE) < 2) return(Inf)

  y_kept <- D$v[w == 1]
  sd(y_kept)
}

root.result <- ROOT(
  data = data,
  global_objective_fn = variance_objective,
  generalizability_path = FALSE
)
summary(root.result)
plot(root.result)
