library(tidyverse)
library(mgcv)
library(splines)


#############################################
#logistic regression for prior positive
results_logistic_p <- list()
model_logistic_p <- list()

for (i in c("pf_one_positive", "az_one_positive")) {
  
  data2 <- data %>% filter(Group==i) %>% 
    select(participant_id, result_tdi,time,age) %>% 
    unique() %>% 
    filter(time<=quantile(.$time[.$time>0],0.9)) %>% 
    group_by(participant_id)%>%
    mutate(ID=cur_group_id()) %>% 
    mutate(ID=as.factor(ID)) %>% 
    ungroup()
  
  xrange <- range(data2$time)
  yrange <- range(data2$age)
  xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                    seq(from = yrange[1], to = yrange[2],1))  
  colnames(xz) <- c("time", "age")
  plot_data <- data.frame(xz) 
  
  m_logistic <- bam(result_tdi~te(time,age,bs="bs",k=c(round(length(unique(data2$time))/3), 
                                                       round(length(unique(data2$age))/3)))+s(ID, bs="re"),method = "fREML", 
                    data = data2, family = binomial, discrete = T, nthreads = 12)
  
  model_logistic_p[[i]] <- m_logistic
  
  #predicing: 
  pred <- predict(m_logistic, newdata=plot_data, se.fit=TRUE, discrete=F,
                  exclude="s(ID)",newdata.guaranteed=TRUE)
  plot_data$prediction <- pred$fit
  plot_data$ll <- (pred$fit-1.96*pred$se.fit)
  plot_data$ul <- (pred$fit+1.96*pred$se.fit)
  plot_data[,c("prediction", "ll", "ul")] <- logit_to_prob(plot_data[,c("prediction", "ll", "ul")])
  plot_data$model <- paste(i)
  
  results_logistic_p[[i]] <- plot_data 
  
}

combine_logistic_p <- do.call(rbind,results_logistic_p)
rownames(combine_logistic_p) <- NULL
combine_logistic_p$model <- as.factor(combine_logistic_p$model)


###########################################
#linear regression for prior positive
results_linear_p <- list()
model_linear_p <- list()

for (i in c("pf_one_positive", "az_one_positive")) {
  
  data2 <- data %>% filter(Group==i) %>% 
    select(participant_id, assay,time,age) %>% 
    filter(time<=quantile(.$time[.$time>0],0.9)) %>% 
    group_by(participant_id)%>%
    mutate(ID=cur_group_id()) %>% 
    mutate(ID=as.factor(ID)) %>% 
    ungroup()
  
  xrange <- range(data2$time)
  yrange <- range(data2$age)
  
  xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                    seq(from = yrange[1], to = yrange[2],1))  
  colnames(xz) <- c("time", "age")
  plot_data <- data.frame(xz) 
  
  m_linear <- bam(log10(assay)~te(time,age,bs="bs",k=c(round(length(unique(data2$time))/3), 
                                                       round(length(unique(data2$age))/3)))+s(ID, bs="re"),method = "fREML", 
                  data = data2, discrete = T, nthreads = 12)
  
  model_linear_p[[i]] <- m_linear
  
  #predicing: 
  pred <- predict(m_linear, newdata=plot_data, se.fit=TRUE, discrete=F,
                  exclude="s(ID)",newdata.guaranteed=TRUE)
  plot_data$prediction <- pred$fit
  plot_data$se <- pred$se.fit
  plot_data$ll <- (pred$fit-1.96*pred$se.fit)
  plot_data$ul <- (pred$fit+1.96*pred$se.fit)
  
  plot_data$prediction2 <- 10^plot_data$prediction
  plot_data$ll2 <- 10^plot_data$ll
  plot_data$ul2 <- 10^plot_data$ul
  
  plot_data$model <- paste(i)
  
  results_linear_p[[i]] <- plot_data 
  
}

combine_linear_p <- do.call(rbind,results_linear_p)
rownames(combine_linear_p) <- NULL
combine_linear_p$model <- as.factor(combine_linear_p$model)


####################################
#logistic regression for prior negative
results_logistic_n <- list()
model_logistic_n <- list()

for (i in c("pf_one_negative", "az_one_negative", "pf_two_negative")) {
  
  data2 <- data %>% filter(Group==i) %>% 
    mutate(time=ifelse(time<0,0,time)) %>% 
    select(participant_id, result_tdi,time,age) %>% 
    unique() %>% 
    filter(time<=quantile(.$time[.$time!=0],0.9)) %>% 
    group_by(participant_id)%>%
    mutate(ID=cur_group_id()) %>% 
    mutate(ID=as.factor(ID)) %>% 
    ungorup()
  
  xrange <- range(data2$time)
  yrange <- range(data2$age)
  
  xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                    seq(from = yrange[1], to = yrange[2],1))  
  colnames(xz) <- c("time", "age")
  plot_data <- data.frame(xz) 
  
  m_logistic <- bam(result_tdi~te(time,age,bs="bs",k=c(round(length(unique(data2$time))/3), 
                                                       round(length(unique(data2$age))/3)))+
                      s(ID, bs="re"),method = "fREML", data = data2, family = binomial,
                    discrete = T, nthreads = 12)
  
  #predicting: 
  pred <- predict(m_logistic, newdata=plot_data, se.fit=TRUE, discrete=F,
                  exclude="s(ID)",newdata.guaranteed=TRUE)
  plot_data$prediction <- pred$fit
  plot_data$ll <- (pred$fit-1.96*pred$se.fit)
  plot_data$ul <- (pred$fit+1.96*pred$se.fit)
  plot_data[,c("prediction", "ll", "ul")] <- logit_to_prob(plot_data[,c("prediction", "ll", "ul")])
  plot_data$model <- paste(i)
  
  results_logistic_n[[i]] <- plot_data 
  
}

combine_logistic_n <- do.call(rbind,results_logistic_n)
rownames(combine_logistic_n) <- NULL
combine_logistic_n$model <- as.factor(combine_logistic_n$model)


################################
#linear regression for prior negative
results_linear_n <- list()
model_linear_n <- list()

for (i in c("pf_one_negative", "az_one_negative", "pf_two_negative")) {
  
  data2 <- data %>% filter(Group==i) %>% 
    group_by(participant_id)%>%
    mutate(ID=cur_group_id()) %>% 
    mutate(ID=as.factor(ID)) %>% 
    arrange(time, .by_group=T) %>% 
    ungroup() %>% 
    mutate(truncate=ifelse(time< -14,1,0))
  
  data2_trunc <- data2 %>% 
    filter(truncate==1) %>% 
    group_by(ID) %>% 
    slice_max(time) %>% 
    ungroup() %>% 
    mutate(time=-14)
  
  data2 <- data2 %>% 
    filter(truncate==0) %>% 
    rbind(data2_trunc) %>% 
    filter(time<=quantile(.$time[.$time>0],0.9)) 
  
  xrange <- range(data2$time)
  yrange <- range(data2$age)
  
  xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                    seq(from = yrange[1], to = yrange[2],1))  
  colnames(xz) <- c("time", "age")
  plot_data <- data.frame(xz) 
  
  m_linear <- bam(log10(assay)~te(time,age,bs="bs",k=c(round(length(unique(data2$time))/3), 
                                                       round(length(unique(data2$age))/3)))+
                    s(ID, bs="re"),method = "fREML", data = data2, discrete = T, nthreads = 12)
  
  #predicting: 
  pred <- predict(m_linear, newdata=plot_data, se.fit=TRUE, discrete=F,
                  exclude="s(ID)",newdata.guaranteed=TRUE)
  plot_data$prediction <- pred$fit
  plot_data$se <- pred$se.fit
  plot_data$ll <- (pred$fit-1.96*pred$se.fit)
  plot_data$ul <- (pred$fit+1.96*pred$se.fit)
  
  plot_data$prediction2 <- 10^plot_data$prediction
  plot_data$ll2 <- 10^plot_data$ll
  plot_data$ul2 <- 10^plot_data$ul
  
  plot_data$model <- paste(i)
  
  results_linear_n[[i]] <- plot_data 
  
}

combine_linear_n <- do.call(rbind,results_linear_n)
rownames(combine_linear_n) <- NULL
combine_linear_n$model <- as.factor(combine_linear_n$model)

