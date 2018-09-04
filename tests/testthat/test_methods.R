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
    expect_identical(rownames(ci2),c("sexF"))
})



# log_lik -----------------------------------------------------------------
test_that("log_lik method works",{
    expect_silent(log_lik(stap_glm1))
    # Compute log-lik matrix using different method than log_lik.stanreg
    # and compare
    samp <- as.matrix(stap_glm1)
    y <- get_y(stap_glm1)
    y_new <- y[1:10] + rnorm(10)
    z <- get_z(stap_glm1)
    z_new <- cbind(1, z[1:10, 2] + rnorm(10))
    sigma <- samp[, 5]
    f1 <- y ~ sex + sap(Fast_Food)
    stap_data <- rstap:::extract_stap_data(f1)
    crs_data <- rstap:::extract_crs_data(stap_data,
                                 subject_data = homog_subject_data,
                                 distance_data = homog_distance_data,
                                 id_key = 'subj_id',max_distance = 3)
    crs_data_new <- rstap:::extract_crs_data(stap_data,
                                 subject_data = homog_subject_data[1:10,],
                                 distance_data = homog_distance_data[which(homog_distance_data$subj_id%in%1:10),],
                                 id_key = 'subj_id',max_distance = 3)
    
    X <- get_x(stap_glm1)
    stap_exp <- apply(apply(X,c(2,3), function(y) y * samp[,3]),c(1,2),sum)
    X_new <- rstap:::.calculate_stap_X(crs_data_new$d_mat, u_s = crs_data_new$u_s,
                              stap_data = stap_data,
                              scales = samp[,4,drop=F])
    X_new <- X_new[1,,1] + rnorm(10)
    eta <- tcrossprod(z, samp[, 1:2]) + t(stap_exp)
    eta_new <- tcrossprod(z_new, samp[, 1:3])
    llmat <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta))
    llmat_new <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta_new))
    for (i in 1:nrow(llmat)) {
        llmat[i, ] <- dnorm(y, mean = eta[, i], sd = sigma[i], log = TRUE)
        llmat_new[i, ] <- dnorm(y_new, mean = eta_new[, i], sd = sigma[i],
                                log = TRUE)
    }
    expect_equal(log_lik(stap_glm1), llmat, check.attributes = FALSE)
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


