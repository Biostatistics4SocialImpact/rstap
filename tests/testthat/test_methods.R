data("homog_distance_data")
data("homog_subject_data")


SW <- suppressWarnings

capture.output(
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
                    chains = 2, iter = 10,
                    refresh = 0)), file = 'NUl')

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
    expect_identical(rownames(ci),c("sexF"))
})



# log_lik -----------------------------------------------------------------
test_that("log_lik method works",{
    
})

# ngrps, nobs -------------------------------------------------------------




# vcov --------------------------------------------------------------------




# sigma -------------------------------------------------------------------




# VarCorr -----------------------------------------------------------------


# ranef,fixef,coef --------------------------------------------------------



# as.matrix,as.data.frame,as.array ----------------------------------------


# terms, formula, model.frame, model.matrix, update methods  --------------




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


# predictive_error,predictive_interval ------------------------------------


# stanreg lists -----------------------------------------------------------


