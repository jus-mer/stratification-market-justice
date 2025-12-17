m15 <- clmm(just_pension ~ 1 + ola_num +  ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             ideo + sex + age + 
             egp*merit_effort_cwc*age + (1 + ola_num + merit_effort_cwc + age| idencuesta), 
          link = "logit",
Hess = TRUE,
         data = df_study1)
m16 <- clmm(just_pension ~ 1 + ola_num +  ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc +
             merit_effort_mean + merit_talent_mean + 
             ideo + sex + age + 
             egp*merit_talent_cwc*age + (1 + ola_num + merit_talent_cwc + age| idencuesta), 
          link = "logit",
Hess = TRUE,
         data = df_study1)
m17 <- clmm(just_pension ~ 1 + ola_num +  ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             ideo + sex + age + 
             egp*merit_effort_mean*age + (1 + ola_num + merit_effort_mean + age| idencuesta), 
          link = "logit",
Hess = TRUE,
         data = df_study1)
m18 <- clmm(just_pension ~ 1 + ola_num +  ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             ideo + sex + age + 
             egp*merit_talent_mean*age + (1 + merit_effort_mean + ola_num | idencuesta), 
          link = "logit",
Hess = TRUE,
         data = df_study1)