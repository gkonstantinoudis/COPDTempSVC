



# Created 05.08.2021


# Simulate data for the models COPD


######################################################################################

library(tidyr)
library(dplyr)

n <- round(279579/12)
controls <- sample(c(3,4), size = n, replace = TRUE)

data.frame(O = 1, controls = controls, ID = 1:length(controls)) -> dat
dat
head(dat)

dat$control1 <- 0
dat$control2 <- 0
dat$control3 <- 0
dat$control4 <- 0

dat$control4[dat$controls %in% 4] <- NA
dat$case <- 1
dat$O <- NULL

set.seed(11)

data_long <- gather(dat, case_control, O, control1:case, factor_key=TRUE)
data_long <- data_long[order(data_long$ID),]
data_long$controls <- NULL
data_long <- data_long[complete.cases(data_long$O),]
data_long$temperature <- abs(rnorm(n = nrow(data_long), mean = 20, sd = 4))
data_long$hol <- rbinom(n = nrow(data_long), size = 1, prob = 3/90)
data_long$O3 <- abs(rnorm(n = nrow(data_long), mean = 70, sd = 20))
data_long$PM25 <- abs(rnorm(n = nrow(data_long), mean = 9, sd = 5))
data_long$RH <- sin(rnorm(n = nrow(data_long), mean = 1.1, sd = 0.2))^2
data_long$ladid <- sample(1:326, size = n, replace = TRUE)[data_long$ID]

# now assign sex and age
age.group <- c("0-64", "65-74", "75+")
sex <- c("male", "female")

dat.sex.age <- data.frame(ID = unique(data_long$ID))
dat.sex.age$age.group <- sample(x = age.group, size = nrow(dat.sex.age), replace = TRUE)
dat.sex.age$sex <- sample(x = sex, size = nrow(dat.sex.age), replace = TRUE)

# and we will assume that 30% of the hospitalizations are recurrent ones.
# sample the IDs that will be the recurrent ones
round(0.3*max(dat.sex.age$ID))# I will do 6990 so its easier to assign
rec.hosp <- sample(size = 6989, x = 1:nrow(dat.sex.age))

dat.sex.age$ID.rec <- dat.sex.age$ID
dat.sex.age$ID.rec[dat.sex.age$ID.rec %in% rec.hosp[1:(6990/2)]] <- rec.hosp[(6990/2+1):6990]

data_long <- left_join(data_long, dat.sex.age, by = c("ID" = "ID"))
data_long$ID.rec <- as.numeric(factor(data_long$ID.rec))
data_long$ID.rec[is.na(data_long$ID.rec)] <- 1
saveRDS(data_long, file = "E:/Postdoc Imperial/misc/COPDTempSVC/dat")

######################################################################################
######################################################################################
######################################################################################


