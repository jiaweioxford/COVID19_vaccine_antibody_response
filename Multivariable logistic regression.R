library(tidyverse)
library(ggeffects)

####### multivariable model
fit_all <- glm(result_tdi~ns(age,knots=c(30,44,57,71),Boundary.knots = c(19,82))*Male+type+ethnicity+IMD+hhsize+ever_patientfacing_hcw+ever_personfacing_socialcare+
                 ever_care_home_worker+ever_lthc, data=data_lr, family = binomial)
summary(fit_all)
pvalue=summary(fit_all)$coefficients[,4]

OR=exp(cbind(coef(fit_all), confint(fit_all))) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(CI=paste0(`2.5 %`,"-",`97.5 %`)) %>% 
  cbind(p=pvalue) %>% 
  select(-`2.5 %`,-`97.5 %`) 

### marginal prediction
p <- ggpredict(fit_all, c("age [all]","Male")) 
plot1 <- plot(p)+labs(title="", x="Age",y="Predicted probability of anti-spike IgG positivity")+
  scale_x_continuous(limits = c(16,85), breaks = c(16,25, 35, 45, 55,65, 75,85))+
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_npg(name="Sex", labels=c("Female","Male"))

p <- ggpredict(fit_all, c("age [all]","type")) 
plot2 <- plot(p)+labs(title="", x="Age",y="Predicted probability of anti-spike IgG positivity")+
  scale_x_continuous(limits = c(16,85), breaks = c(16,25, 35,45, 55, 65, 75, 85))+
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_fill_npg(labels=c("One dose Pfizer","One dose AZ","Two doses Pfizer"))+
  scale_color_npg(labels=c("One dose Pfizer","One dose AZ","Two doses Pfizer"))

####### univariable model
var_list <- c("ns(age,knots=c(30,44,57,71),Boundary.knots = c(19,82))","type","Male","ethnicity","hhsize","IMD",
              "ever_patientfacing_hcw", "ever_personfacing_socialcare","ever_care_home_worker","ever_lthc")

results <- list()
for (var in var_list){
  fit_uni <- glm(paste0("result_tdi~",var), data=data_lr, family = binomial)
  pvalue=summary(fit_uni)$coefficients[,4]
  OR=exp(cbind(coef(fit_uni), confint(fit_uni))) %>% 
    as.data.frame() %>% 
    round(2) %>% 
    rownames_to_column() %>% 
    mutate(CI=paste0(`2.5 %`,"-",`97.5 %`)) %>% 
    cbind(p=pvalue) %>% 
    select(-`2.5 %`,-`97.5 %`) %>% 
    slice(-1)
  
  results[[var]] <- OR
} 

combine_LR <- do.call(rbind,results)
rownames(combine_LR) <- NULL
