## code to prepare `FFR_example` dataset goes here

library(rbenvo)
library(dplyr)
set.seed(341438)
n <- 100
parametric_f <- function(x) pweibull(q = x,shape = 5,scale = .5,lower.tail = FALSE)
has_exp <- rbinom(n,1,prob=.95)
num_dists <- rpois(n = n,lambda = 10)
ldists <- lapply(num_dists,function(x) runif(x))
exp <- purrr::map_dbl(ldists,function(y) sum(parametric_f(y)))
exp <- exp*has_exp
Z <- rbinom(n,1,.5)

y <- 26 + Z *-1.5 + exp + rnorm(n,sd=1)

sdf <- tibble(subject_id = 1:n,
              BMI = y,
              sex = factor(Z,labels=c("Male","Female")))

ddf <- purrr::map2_dfr(1:n,ldists,function(x,y) tibble(subject_id = x,
                                                       Distance = y))

FFR_example <- rbenvo::benvo(subject_data = sdf,sub_bef_data = list(FFR=ddf),by="subject_id")


usethis::use_data(FFR_example, overwrite = TRUE)
