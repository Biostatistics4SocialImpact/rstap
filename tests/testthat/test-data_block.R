context("test internal coding functions")
f1 <- y ~ sex + stap(Fast_Food)
f2 <- y ~ Age + stap(Fast_Food) + sap(Coffee_Shops)
f3 <- y ~ Age + stap(Fast_Food,cerf,exp) + sex + tap(Coffee_Shops)
m3 <- rbind(c(2,3),c(0,1))
a1 <- c("Fast_Food"=2)
a2 <- c(a1,"Coffee_Shops"=0)
a3 <- c("Fast_Food"=2,"Coffee_Shops"=1)

test_that("correctly assigns weights",{
    expect_equal(get_weight_code(all.names(f1),'Fast_Food',c(2)),
                 matrix(c(2,1),nrow=1))
    expect_equal(get_weight_code(all.names(f2),c('Fast_Food',"Coffee_Shops"),c(2,0)),
                 matrix(c(c(2,2),c(1,0)),nrow=2))
    expect_equal(get_weight_code(all.names(f3),c("Fast_Food","Coffee_Shops"),c(2,1)),
                 m3)
})

test_that("correctly assigns stap coding",{
    expect_equal(get_stap_code(all.names(f1),'Fast_Food'),a1)
    expect_equal(get_stap_code(all.names(f2),c("Fast_Food","Coffee_Shops")),a2)
    expect_equal(get_stap_code(all.names(f3),c("Fast_Food","Coffee_Shops")),a3)
})

test_that("weight_switch works",{
    expect_equal(weight_switch(0),"none")
    expect_equal(weight_switch(1),"erf")
    expect_equal(weight_switch(2),"cerf")
    expect_equal(weight_switch(3),"exp")
    expect_equal(weight_switch(4),"cexp")
})

test_that("get_weight_name works",{
    expect_equal(get_weight_name(c(0,1)),list('spatial' = 'none',
                                              'temporal' = 'erf'))
    expect_equal(get_weight_name(c(2,0)),list("spatial" = "cerf",
                                              "temporal" = "none"))
    expect_equal(get_weight_name(c(2,3)), list("spatial" = "cerf",
                                               "temporal" = "exp"))
})

context("test stap data extraction functions")

f1 <- y ~ sex + Age + stap(Coffee_Shops)
f2 <- y ~ tap_log(Con_Stores,erf) + sex + Age + stap(Coffee_Shops)


f1 <- BMI ~ Age +  sap(Coffee_Shops)
f2 <- BMI ~ Age + tap(Coffee_Shops)
distance_data <- data.frame(subj_id = 1:10,
                            BEF = rep("Coffee_Shops",10),
                            dist = rexp(10))
subj_data <- data.frame(subj_id = 1:10,
                        BMI = rnorm(10,mean = 25, sd = 2),
                        Age = rnorm(10,mean = 35, sd = 10))
admat <- matrix(distance_data$dist,nrow=1)
rownames(admat) <- "Coffee_Shops"
a1 <-  list(d_mat = admat,
            t_mat = NA,
            u_t = NA,
            u_s = as.array(cbind(1:10,1:10),dim=c(10,2,1)))
a2 <- list(d_mat = NA,
           t_mat = admat,
           u_t = as.array(cbind(1:10,1:10),dim=c(10,2,1)),
           u_s = NA)
stap_data_1 <- extract_stap_data(f1)
stap_data_2 <- extract_stap_data(f2)

test_that("extract_crs_data correctly errors when no distance or time data are given",{
    expect_error(extract_crs_data(formula = y ~ X,
                                         id_key = 'subj_id',
                                         max_distance = 3))
})

test_that("extract_crs_data correctly extracts data",{
    expect_equal(extract_crs_data(stap_data_1,subject_data = subj_data,
                                  distance_data = distance_data,
                                  time_data = NULL,
                                  id_key = 'subj_id',
                                  max_distance = max(distance_data$dist)),
                 a1)
    expect_equal(extract_crs_data(stap_data_2,subj_data, time_data = distance_data,
                                  id_key = 'subj_id',
                                  max_distance = max(distance_data$dist)),
                 a2)
})
