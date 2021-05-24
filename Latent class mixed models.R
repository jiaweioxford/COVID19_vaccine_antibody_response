library(tidyverse)
library(lcmm)

###################latent class, Oxford-AZ
data_az <- as.data.frame(data_az)

# add in time non-linear:
time <- data.frame(splines::ns(data_az$time,knots = c(0,12,24,37),Boundary.knots = c(-11,47)))
names(time) <- c("TIME1","TIME2","TIME3","TIME4","TIME5")
data_az <- cbind(data_az,time)

# add in age non-linear: 
age <- data.frame(splines::ns(data_az$age,knots = c(30,44,57,71),Boundary.knots = c(19,82)))
names(age) <- c("AGE1","AGE2","AGE3","AGE4","AGE5")
data_az <- cbind(data_az,age)

#baseline
model_az <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                 random = ~ 1, subject='ID',
                 data = data_az)

#two classes
model_az2 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~ AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male, 
                  random = ~ 1, subject='ID', ng=2, B=model_az,
                  data = data_az)

summary(model_az2)
postprob(model_az2)

#three classes
model_az3 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~ AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male, 
                  random = ~ 1, subject='ID', ng=3, B=model_az,
                  data = data_az)

summary(model_az3)
postprob(model_az3)

#four classes
model_az4 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male,maxiter = 500,
                  random = ~ 1, subject='ID', ng=4, B=model_az,
                  data = data_az)

summary(model_az4)
postprob(model_az4)

###################latent class-Pfizer
data_pf <- as.data.frame(data_pf)

# add in time non-linear:
time <- data.frame(splines::ns(data_pf$time,knots = c(2,18,34,50), Boundary.knots = c(-10,62)))
names(time) <- c("TIME1","TIME2","TIME3","TIME4","TIME5")
data_pf <- cbind(data_pf,time)

# add in age non-linear: 
age <- data.frame(splines::ns(data_pf$age,knots = c(30,44,57,71),Boundary.knots = c(19,82)))
names(age) <- c("AGE1","AGE2","AGE3","AGE4","AGE5")
data_pf <- cbind(data_pf,age)

#baseline
model_pf <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                 random = ~ 1, subject='ID',
                 data = data_pf)

#two classes
model_pf2 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male,
                  random = ~ 1, subject='ID', ng=2, B=model_pf,
                  data = data_pf)

summary(model_pf2)
postprob(model_pf2)

#three classes
model_pf3 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male,
                  random = ~ 1, subject='ID', ng=3, B=model_pf,
                  data = data_pf)

summary(model_pf3)
postprob(model_pf3)

#four class
model_pf4 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5, 
                  mixture = ~ TIME1 + TIME2 + TIME3 + TIME4 + TIME5,
                  classmb = ~AGE1 + AGE2 + AGE3 + AGE4 + AGE5 + 
                    ever_lthc+ever_patientfacing_hcw+Male, maxiter = 500,
                  random = ~ 1, subject='ID', ng=4, B=model_pf,
                  data = data_pf)

summary(model_pf4)
postprob(model_pf4)