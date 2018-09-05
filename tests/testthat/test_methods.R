data("homog_distance_data")
data("homog_subject_data")


SW <- suppressWarnings

capture.output(
    ITER <- 10,
    CHAINS <- 2,
    REFRESH <- 0,
    stap_glm1 <- SW(stap_glm(formula = y ~ sex + sap(Fast_Food),
                    subject_data = homog_subject_data,
                    distance_data = homog_distance_data,
                    family = gaussian(link = 'identity'),
                    id_key = 'subj_id',
                    prior = normal(location = 0, scale = 5,autoscale = F),
                    prior_intercept = normal(location = 25, scale = 5, autoscale = F),
                    prior_stap = normal(location = 0, scale = 3, autoscale = F),
                    prior_theta = log_normal(location = 1, scale = 1),
                    prior_aux = cauchy(location = 0,scale = 5),
                    max_distance = max(homog_distance_data$Distance),
                    chains = CHAINS, iter = ITER,
                    refresh = REFRESH)),
    glm1 <- glm(y ~ sex,data=homog_subject_data,family=gaussian), file = 'NULL')
    

context("methods for stanreg objects")


# extractors --------------------------------------------------------------
test_that("stapreg extractor methods work properly", {
    expect_equal(resid(stap_glm1), stap_glm1$residuals)
    expect_equal(coef(stap_glm1), stap_glm1$coefficients)
    expect_equal(vcov(stap_glm1), stap_glm1$covmat)
    expect_equal(fitted(stap_glm1), stap_glm1$fitted.values)
    expect_equal(se(stap_glm1), stap_glm1$ses)
    
    expect_error(confint(stap_glm1), regexp = "use posterior_interval")
})


# posterior_interval ------------------------------------------------------
test_that("posterior_interval returns correct structure",{
    expect_silent(ci <- posterior_interval(stap_glm1, prob = 0.5))
    expect_silent(ci2 <- posterior_interval(stap_glm1, prob = 0.95, regex_pars = "sex"))
    expect_identical(rownames(ci),c("(Intercept)","sexF","Fast_Food","Fast_Food_spatial_scale","sigma"))
    expect_identical(rownames(ci2),c("sexF"))
})



# log_lik -----------------------------------------------------------------
test_that("log_lik method works",{
    expect_silent(log_lik(stap_glm1))
    # Compute log-lik matrix using different method than log_lik.stanreg
    # and compare
    samp <- as.matrix(stap_glm1)
    y <- get_y(stap_glm1)
    z <- get_z(stap_glm1)
    sigma <- samp[, 5]
    f1 <- y ~ sex + sap(Fast_Food)
    stap_data <- rstap:::extract_stap_data(f1)
    crs_data <- rstap:::extract_crs_data(stap_data,
                                 subject_data = homog_subject_data,
                                 distance_data = homog_distance_data,
                                 id_key = 'subj_id',max_distance = 3)
    X <- get_x(stap_glm1)
    stap_exp <- t(apply(apply(X,c(2,3), function(y) y * samp[,3]),c(1,2),sum))
    eta <- tcrossprod(z, samp[, 1:2]) + stap_exp
    llmat <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta))
    for (i in 1:nrow(llmat)) 
        llmat[i, ] <- dnorm(y, mean = eta[, i], sd = sigma[i], log = TRUE)
    expect_equal(log_lik(stap_glm1), llmat, check.attributes = FALSE)
    
    ## Need to figure out how to handle "newdata" both for predictions and log_lik functions
})


# ngrps, nobs -------------------------------------------------------------


test_that("ngrps is right", {
    expect_error(ngrps(stap_glm1), "stan_glmer and stan_lmer models only")
})


test_that("nobs is right", {
    expect_equal(nobs(stap_glm1), nrow(homog_subject_data))
})
# vcov --------------------------------------------------------------------

test_that("vcov returns correct structure", {
    expect_equal(rownames(vcov(stap_glm1)), c("(Intercept)","sexF","Fast_Food","Fast_Food_spatial_scale"))
})

# sigma -------------------------------------------------------------------

test_that("sigma method works",{
    expect_double <- function(x) expect_type(x, "double")
    expect_double(sig <- sigma(stap_glm1))
    expect_false(identical(sig, 1))
})


# VarCorr -----------------------------------------------------------------
test_that("VarCorr returns correct structure", {
    expect_error(VarCorr(stap_glm1), "stan_glmer and stan_lmer models only")
})

# ranef,fixef,coef --------------------------------------------------------
test_that("ranef returns correct structure", {
    expect_error(ranef(stap_glm1), "stan_glmer and stan_lmer models only")
})

# as.matrix,as.data.frame,as.array ----------------------------------------

test_that("as.matrix, as.data.frame, as.array methods work for MCMC", {
    last_dimnames <- rstap:::last_dimnames
    # glm
    mat <- as.matrix(stap_glm1)
    df <- as.data.frame(stap_glm1)
    arr <- as.array(stap_glm1)
    expect_identical(df, as.data.frame(mat))
    expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
    expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, 5L))
    expect_equal(dim(arr), c(floor(ITER/2), CHAINS, 5L))
    expect_identical(last_dimnames(mat), c("(Intercept)", "sexF", "Fast_Food", "Fast_Food_spatial_scale",
                                           "sigma"))
    expect_identical(last_dimnames(arr), last_dimnames(mat))
    
    # selecting only 1 parameter
    mat <- as.matrix(stap_glm1, pars = "sexF")
    df <- as.data.frame(stap_glm1, pars = "sexF")
    arr <- as.array(stap_glm1, pars = "sexF")
    expect_identical(df, as.data.frame(mat))
    expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
    expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, 1L))
    expect_equal(dim(arr), c(floor(ITER/2), CHAINS, 1L))
    expect_identical(last_dimnames(mat), "sexF")
    expect_identical(last_dimnames(arr), last_dimnames(mat))
    
    # pars & regex_pars
    nr <- posterior_sample_size(stap_glm1)
    mat <- as.matrix(stap_glm1, pars = "sexF", regex_pars = "Fast_Food")
    df <- as.data.frame(stap_glm1, pars = "sexF", regex_pars = "Fast_Food")
    arr <- as.array(stap_glm1, pars = "sexF", regex_pars = "Fast_Food")
    expect_identical(df, as.data.frame(mat))
    expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
    expect_equal(dim(mat), c(nr, 3L))
    expect_equal(dim(arr), c(nr/2, 2, 3L))
    expect_identical(last_dimnames(mat), c("sexF", "Fast_Food","Fast_Food_spatial_scale"))
    expect_identical(last_dimnames(mat), last_dimnames(arr))
})


# terms, formula, model.frame, model.matrix methods  --------------

context("model.frame methods")
test_that("model.frame works properly", {
    expect_identical(model.frame(stap_glm1), model.frame(glm1))
})

context("terms methods")
test_that("terms works properly", {
    expect_identical(terms(stap_glm1), terms(glm1))
})

context("formula methods")
test_that("formula works properly", {
    expect_equal(formula(stap_glm1),as.formula(y~sex + sap(Fast_Food)))
})

## no update method yet
    

# print and summary -------------------------------------------------------

context("print and summary methods")

test_that("print and summary methods work ok for stap_glm", {
    expect_output(print(stap_glm1,digits = 2), "stap_glm")
    expect_output(print(stap_glm1,digits = 2), 'sigma')
    expect_silent(s <- summary(stap_glm1))
    expect_silent(d <- as.data.frame(s))
    expect_s3_class(s,"summary.stapreg")
    expect_identical(colnames(s),colnames(d))
    expect_identical(rownames(s),rownames(d))
})


# prior_summary -----------------------------------------------------------
test_that("prior_summary doesn't error", {
    expect_output(print(prior_summary(stap_glm1, digits = 2)),
                  "Priors for model 'stap_glm1'")
})

test_that("prior_summary returns correctly named list", {
    expect_named(prior_summary(stap_glm1),
                 c("prior", "prior_stap", "prior_theta", "prior_intercept", "prior_aux"))
})


# predictive_error,predictive_interval ------------------------------------




