# 0. Identification ---------------------------------------------------

# Title: Data analysis for paper Changes in the Justification of Pension Inequality in Chile (2016–2023) and its Relationship to Social Class and Beliefs in Meritocracy
# Responsible: Andreas Laffert

# Executive Summary: This script contains the code to data analysis for the paper
# Date: May 1, 2025

# 1. Packages  -----------------------------------------------------

if (! require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse, 
               sjmisc, 
               sjPlot, 
               lme4, 
               here, 
               performance,
               influence.ME, 
               srvyr,
               texreg, 
               ggdist,
               misty,
               kableExtra,
               ggalluvial, 
               shadowtext,
               MetBrewer,
               patchwork,
               sjlabelled)


options(scipen=999)
rm(list = ls())

# 2. Data -----------------------------------------------------------------

load(file = here("input/data/proc/df_study1_long_t7.RData"))

# Generate analytical sample

df_study1 <- df_study1_long_t7 %>%
  select(-muestra) %>% 
  na.omit() %>% 
  mutate(ola = case_when(ola == 1 ~ 1,
                         ola == 2 ~ 2, 
                         ola == 3 ~ 3,
                         ola == 4 ~ 4,
                         ola == 6 ~ 5,
                         ola == 7 ~ 6)) %>% 
  mutate(ola = as.factor(ola),
         ola_num = as.numeric(ola),
         ola_2=as.numeric(ola)^2,
         wave = case_when(ola == 1 ~ "2016",
                          ola == 2 ~ "2017",
                          ola == 3 ~ "2018",
                          ola == 4 ~ "2019",
                          ola == 5 ~ "2022",
                          ola == 6 ~ "2023"),
         wave = factor(wave, levels = c("2016",
                                        "2017",
                                        "2018",
                                        "2019",
                                        "2022",
                                        "2023")))

df_study1 <- df_study1 %>%
  group_by(idencuesta) %>%             # Agrupar por el identificador del participante
  mutate(n_participaciones = n()) %>%  # Contar el número de filas (participaciones) por participante
  ungroup()

df_study1 <- df_study1 %>% filter(n_participaciones>1)

# Corregir etiquetas

df_study1$just_pension <- sjlabelled::set_label(df_study1$just_pension, 
                                                label = "Pension distributive justice")

df_study1$egp <- sjlabelled::set_label(df_study1$egp, 
                                       label = "Social class")

df_study1$merit_effort <- sjlabelled::set_label(df_study1$merit_effort, 
                                                label = "People are rewarded for their efforts")

df_study1$merit_talent <- sjlabelled::set_label(df_study1$merit_talent, 
                                                label = "People are rewarded for their intelligence")

df_study1$wave <- sjlabelled::set_label(df_study1$wave, 
                                        label = "Wave")

# 3. Analysis -----------------------------------------------------------

# 3.1 Descriptives --------------------------------------------------------

## Alluvial
datos.pension <- df_study1 %>% 
  mutate(just_pension = factor(just_pension, 
                               levels = c("Strongly agree",
                                          "Agree",
                                          "Neither agree nor disagree",
                                          "Disagree",
                                          "Strongly disagree"))) %>% 
  group_by(idencuesta, wave) %>% 
  count(just_pension) %>% 
  group_by(wave) %>% 
  mutate(porcentaje=n/sum(n)) %>% 
  ungroup() %>% 
  na.omit()



etiquetas.pension <- df_study1 %>%
  mutate(just_pension = factor(just_pension, 
                               levels = c("Strongly agree",
                                          "Agree",
                                          "Neither agree nor disagree",
                                          "Disagree",
                                          "Strongly disagree"))) %>% 
  group_by(wave, just_pension) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(wave) %>%
  mutate(porcentaje = count / sum(count)) %>% 
  na.omit() %>% 
  mutate(idencuesta = 1)


g1 <- datos.pension %>% 
  ggplot(aes(x = wave, fill = just_pension, stratum = just_pension,
             alluvium = idencuesta, y = porcentaje)) +
  ggalluvial::geom_flow(alpha = .4) + 
  ggalluvial::geom_stratum(linetype = 0) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values =  c("#0571B0","#92C5DE","#b3b3b3ff","#F4A582","#CA0020")) +
  geom_shadowtext(data = etiquetas.pension,
                  aes(label = ifelse(porcentaje > 0 , scales::percent(porcentaje, accuracy = .1),"")),
                  position = position_stack(vjust = .5),
                  show.legend = FALSE,
                  size = 3,
                  color = rep('white'),
                  bg.colour='grey30')+
  labs(y = NULL,
       x = NULL,
       fill = NULL,
       title = "It is fair that people with higher incomes have better pensions than people\nwith lower incomes?",
       caption = "Source: own elaboration with data from ELSOC 2016-2023 (N obs = 5,755)") +
  theme_ggdist() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        text = element_text(size = 16)) 

g1

## Summary table

t2 <- df_study1 %>% 
  filter(ola == last(ola)) %>% 
  select(just_pension, egp, merit_effort, merit_talent) 

summarytools::dfSummary(t2)

## Mean Z change

df_pond <- df_study1 %>% 
  as_survey_design(.data = .,
                   ids = segmento, 
                   strata = estrato, 
                   weights = ponderador_long_total)


df_pond <- df_pond %>% 
  select(wave, egp, just_pension,  merit_effort, merit_talent) %>% 
  mutate_at(.vars = 3:5, .funs = ~ as.numeric(.)) %>% 
  mutate(
    just_pension_z    = as.numeric(scale(just_pension)),
    merit_effort_z    = as.numeric(scale(merit_effort)),
    merit_talent_z    = as.numeric(scale(merit_talent))
  )


g2 <- df_pond %>% 
  select(wave, egp, just_pension_z,  merit_effort_z, merit_talent_z) %>% 
  group_by(wave, egp) %>% 
  summarise_all(~survey_mean(., vartype = "ci")) %>% 
  pivot_longer(cols = -c(wave, egp),
               names_to = "temp",
               values_to = "valor") %>%
  mutate(
    ci = case_when(
      str_ends(temp, "_low") ~ "ic_low",
      str_ends(temp, "_upp") ~ "ic_upp",
      TRUE                   ~ "mean"
    ),
    variable = str_remove(temp, "_low|_upp")
  ) %>%
  select(wave, egp, variable, ci, valor) %>%
  pivot_wider(
    names_from  = ci,
    values_from = valor
  ) %>% 
  mutate(variable = case_when(variable == "just_pension_z" ~ "Justification of pension inequality",
                              variable == "merit_effort_z"     ~ "Merit: Effort",
                              variable == "merit_talent_z"     ~ "Merit: Talent"),
         variable = factor(variable, levels = c("Justification of pension inequality",
                                                "Merit: Effort",
                                                "Merit: Talent")
         )) %>% 
  ggplot(aes(x = wave, y = mean, group = egp)) +
  geom_point(aes(shape=egp, color=egp), size = 3.5) +
  geom_line(aes(color = egp), linewidth = 0.8) +
  geom_errorbar(aes(ymin = ic_low, ymax = ic_upp, color = egp),
                width = 0.1) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 5) +
  facet_wrap(~variable, nrow = 3) +
  labs(y = "Mean (Z-score)",
       x = "Wave",
       color = NULL,
       shape = NULL,
       caption = "Source: own elaboration with pooled data from ELSOC 2016-2023 (N obs = 5,755; N groups = 1,027)") +
  theme_ggdist() +
  theme(legend.position = "top",
        text = element_text(size = 12)) 


g2

# 3.2 Longitudinal multilevel models ---------------------------------------------------

df_study1$just_pension <- as_numeric(df_study1$just_pension)

df_study1$merit_effort <- as_numeric(df_study1$merit_effort)
df_study1$merit_talent <- as_numeric(df_study1$merit_talent)

df_study1$merit_effort_dic <- if_else(df_study1$merit_effort <= 3, 0, 1)
df_study1$merit_talent_dic <- if_else(df_study1$merit_talent <= 3, 0, 1)

df_study1$age_2 <- (df_study1$age)^2 

df_study1 <- df_study1 %>% 
  group_by(idencuesta) %>% 
  mutate(merit_effort_mean = mean(merit_effort, na.rm = T),
         merit_effort_cwc = merit_effort - merit_effort_mean,
         merit_talent_mean = mean(merit_talent, na.rm = T),
         merit_talent_cwc = merit_talent - merit_talent_mean,
         merit_effort_dic_mean = mean(merit_effort_dic, na.rm = T),
         merit_effort_dic_cwc = merit_effort_dic - merit_effort_dic_mean,
         merit_talent_dic_mean = mean(merit_talent_dic, na.rm = T),
         merit_talent_dic_cwc = merit_talent_dic - merit_talent_dic_mean,
  ) %>% 
  ungroup()

df_study1$egp <- factor(df_study1$egp, levels = c("Working class (V+VI+VII)",
                                                  "Intermediate class (III+IV)",
                                                  "Service class (I+II)"))

df_study1$ola <- df_study1$wave

## ICC 

m0 <- lmer(just_pension ~ 1 + (1 | idencuesta), 
            data = df_study1)

performance::icc(m0, by_group = T) # 0.23 es between, 0.77 within

## Time effects

m1.1 <- lmer(just_pension ~ 1 + ola + (1 | idencuesta),
              data = df_study1)

m1.2 <- lmer(just_pension ~ 1 + ola_num + (1 | idencuesta),
              data = df_study1)

m1.3 <- lmer(just_pension ~ 1 + ola_num + ola_2 + (1| idencuesta),
              data = df_study1)

m1.4 <- lmer(just_pension ~ 1 + ola_num + ola_2 + (1 + ola_num | idencuesta),
              data = df_study1)

ccoef <- list(
  "(Intercept)" = "Intercept", 
  "ola2017" = "Wave 2017",
  "ola2018" = "Wave 2018",
  "ola2019" = "Wave 2019",
  "ola2022" = "Wave 2022",
  "ola2023" = "Wave 2023",
  ola_num = "Wave",
  ola_2 = "Wave^2")

texreg::screenreg(list(m1.1,m1.2,m1.3,m1.4),
                  caption.above = T,
                  caption = NULL,
                  stars = c(0.05, 0.01, 0.001),
                  custom.coef.map = ccoef,
                  digits = 3,
                  groups = list("Wave (Ref. = 2016)" = 2:6),
                  custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                  leading.zero = T,
                  use.packages = F,
                  booktabs = F,
                  scalebox = 0.80,
                  include.loglik = FALSE,
                  include.aic = FALSE,
                  center = T)

anova(m1.3, m1.4) # quedarse con tiempo continua, cuadratica y con pendiente aleatoria

## WE and BE main effects

m3 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + (1 + ola_num | idencuesta),
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

m4 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
           merit_effort_cwc + merit_talent_cwc + (1 + ola_num | idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

m5 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
           merit_effort_cwc + merit_talent_cwc + 
           merit_effort_mean + merit_talent_mean + (1 + ola_num| idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

m6 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
           merit_effort_cwc + merit_talent_cwc + 
           merit_effort_mean + merit_talent_mean + 
           educ + ideo + sex + age + (1 + ola_num| idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2",
  "egpIntermediate class (III+IV)" = "Intermediate class (III+IV)",
  "egpService class (I+II)" = "Service class (I+II)",
  merit_effort_cwc = "Merit: Effort (WE)",
  merit_talent_cwc = "Merit: Talent (WE)",
  merit_effort_mean = "Merit: Effort (BE)",
  merit_talent_mean = "Merit: Talent (BE)",
  "educUniversitary" = "Universitary education (Ref.= Less than Universitary)",
  ideoCenter = "Center",
  ideoRight = "Right",
  "ideoDoes not identify" = "Does not identify",
  sexFemale = "Female (Ref.= Male)",
  age = "Age"
)

texreg::screenreg(list(m3,m4,m5,m6),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                groups = list("Social Class (Ref.= Working class (V+VI+VII))" = 4:5,
                              "Political identification (Ref.= Left)" =      11:13),
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T)


## WE and BE Interactions without controls

m7 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             egp*merit_effort_cwc + (1 + ola_num + merit_effort_cwc| idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

m8 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             egp*merit_talent_cwc + (1 + ola_num + merit_talent_cwc| idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))


m9 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
             merit_effort_cwc + merit_talent_cwc + 
             merit_effort_mean + merit_talent_mean + 
             egp*merit_effort_mean + (1 + ola_num + merit_effort_mean| idencuesta), 
           data = df_study1,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 1e5)))

m10 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              egp*merit_talent_mean + (1 + ola_num + merit_talent_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2",
  "egpIntermediate class (III+IV)" = "Intermediate class (III+IV)",
  "egpService class (I+II)" = "Service class (I+II)",
  merit_effort_cwc = "Merit: Effort (WE)",
  merit_talent_cwc = "Merit: Talent (WE)",
  merit_effort_mean = "Merit: Effort (BE)",
  merit_talent_mean = "Merit: Talent (BE)",
  "egpIntermediate class (III+IV):merit_effort_cwc" =  "Intermediate class (III+IV) x Merit: Effort (WE)",
  "egpService class (I+II):merit_effort_cwc" = "Service class (I+II) x Merit: Effort (WE)",
  "egpIntermediate class (III+IV):merit_talent_cwc" = "Intermediate class (III+IV) x Merit: Talent (WE)",
  "egpService class (I+II):merit_talent_cwc" = "Service class (I+II) x Merit: Talent (WE)",
  "egpIntermediate class (III+IV):merit_effort_mean" =  "Intermediate class (III+IV) x Merit: Effort (BE)",
  "egpService class (I+II):merit_effort_mean" = "Service class (I+II) x Merit: Effort (BE)",
  "egpIntermediate class (III+IV):merit_talent_mean" = "Intermediate class (III+IV) x Merit: Talent (BE)",
  "egpService class (I+II):merit_talent_mean" = "Service class (I+II) x Merit: Talent (BE)")

texreg::screenreg(list(m7,m8,m9,m10),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                groups = list("Social Class (Ref.= Working class (V+VI+VII))" = 4:5,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (WE)" =      10:13,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (BE)" =      14:17),
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("No",4))))

## Interactions with controls (direct effect)

## WE and BE Interactions with controls

m11 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              ideo + sex + age + 
              egp*merit_effort_cwc + (1 + ola_num + merit_effort_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m12 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              ideo + sex + age + 
              egp*merit_talent_cwc + (1 + ola_num + merit_talent_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


m13 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              ideo + sex + age + 
              egp*merit_effort_mean + (1 + ola_num + merit_effort_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m14 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              ideo + sex + age + 
              egp*merit_talent_mean + (1 + ola_num + merit_talent_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2", 
  "egpIntermediate class (III+IV)" = "Intermediate class (III+IV)",
  "egpService class (I+II)" = "Service class (I+II)",
  merit_effort_cwc = "Merit: Effort (WE)",
  merit_talent_cwc = "Merit: Talent (WE)",
  merit_effort_mean = "Merit: Effort (BE)",
  merit_talent_mean = "Merit: Talent (BE)",
  "egpIntermediate class (III+IV):merit_effort_cwc" =  "Intermediate class (III+IV) x Merit: Effort (WE)",
  "egpService class (I+II):merit_effort_cwc" = "Service class (I+II) x Merit: Effort (WE)",
  "egpIntermediate class (III+IV):merit_talent_cwc" = "Intermediate class (III+IV) x Merit: Talent (WE)",
  "egpService class (I+II):merit_talent_cwc" = "Service class (I+II) x Merit: Talent (WE)",
  "egpIntermediate class (III+IV):merit_effort_mean" =  "Intermediate class (III+IV) x Merit: Effort (BE)",
  "egpService class (I+II):merit_effort_mean" = "Service class (I+II) x Merit: Effort (BE)",
  "egpIntermediate class (III+IV):merit_talent_mean" = "Intermediate class (III+IV) x Merit: Talent (BE)",
  "egpService class (I+II):merit_talent_mean" = "Service class (I+II) x Merit: Talent (BE)")

texreg::screenreg(list(m11,m12,m13,m14),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                groups = list("Social Class (Ref.= Working class (V+VI+VII))" = 4:5,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (WE)" =      10:13,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (BE)" =      14:17),
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("Yes",4))))


## Interactions with meritocratic variables as dichotomized variables without controls


# WE and BE Interactions merit dummy

m15 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              egp*merit_effort_dic_cwc + (1 + ola_num + merit_effort_dic_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m16 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              egp*merit_talent_dic_cwc + (1 + ola_num + merit_talent_dic_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


m17 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              egp*merit_effort_dic_mean + (1 + ola_num + merit_effort_dic_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m18 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              egp*merit_talent_dic_mean + (1 + ola_num + merit_talent_dic_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2", 
  "egpIntermediate class (III+IV)" = "Intermediate class (III+IV)",
  "egpService class (I+II)" = "Service class (I+II)",
  merit_effort_dic_cwc = "Merit: Effort (WE)",
  merit_talent_dic_cwc = "Merit: Talent (WE)",
  merit_effort_dic_mean = "Merit: Effort (BE)",
  merit_talent_dic_mean = "Merit: Talent (BE)",
  "egpIntermediate class (III+IV):merit_effort_dic_cwc" =  "Intermediate class (III+IV) x Merit: Effort (WE)",
  "egpService class (I+II):merit_effort_dic_cwc" = "Service class (I+II) x Merit: Effort (WE)",
  "egpIntermediate class (III+IV):merit_talent_dic_cwc" = "Intermediate class (III+IV) x Merit: Talent (WE)",
  "egpService class (I+II):merit_talent_dic_cwc" = "Service class (I+II) x Merit: Talent (WE)",
  "egpIntermediate class (III+IV):merit_effort_dic_mean" =  "Intermediate class (III+IV) x Merit: Effort (BE)",
  "egpService class (I+II):merit_effort_dic_mean" = "Service class (I+II) x Merit: Effort (BE)",
  "egpIntermediate class (III+IV):merit_talent_dic_mean" = "Intermediate class (III+IV) x Merit: Talent (BE)",
  "egpService class (I+II):merit_talent_dic_mean" = "Service class (I+II) x Merit: Talent (BE)")

texreg::screenreg(list(m15,m16,m17,m18),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                groups = list("Social Class (Ref.= Working class (V+VI+VII))" = 4:5,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (WE)" =      10:13,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (BE)" =      14:17),
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("No",4))))



## Interactions with meritocratic variables as dichotomized variables with controls


# WE and BE Interactions merit dummy

m19 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              ideo + sex + age + 
              egp*merit_effort_dic_cwc + (1 + ola_num + merit_effort_dic_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m20 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              ideo + sex + age + 
              egp*merit_talent_dic_cwc + (1 + ola_num + merit_talent_dic_cwc| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


m21 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              ideo + sex + age + 
              egp*merit_effort_dic_mean + (1 + ola_num + merit_effort_dic_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m22 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_dic_cwc + merit_talent_dic_cwc + 
              merit_effort_dic_mean + merit_talent_dic_mean + 
              ideo + sex + age + 
              egp*merit_talent_dic_mean + (1 + ola_num + merit_talent_dic_mean| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2",  
  "egpIntermediate class (III+IV)" = "Intermediate class (III+IV)",
  "egpService class (I+II)" = "Service class (I+II)",
  merit_effort_dic_cwc = "Merit: Effort (WE)",
  merit_talent_dic_cwc = "Merit: Talent (WE)",
  merit_effort_dic_mean = "Merit: Effort (BE)",
  merit_talent_dic_mean = "Merit: Talent (BE)",
  "egpIntermediate class (III+IV):merit_effort_dic_cwc" =  "Intermediate class (III+IV) x Merit: Effort (WE)",
  "egpService class (I+II):merit_effort_dic_cwc" = "Service class (I+II) x Merit: Effort (WE)",
  "egpIntermediate class (III+IV):merit_talent_dic_cwc" = "Intermediate class (III+IV) x Merit: Talent (WE)",
  "egpService class (I+II):merit_talent_dic_cwc" = "Service class (I+II) x Merit: Talent (WE)",
  "egpIntermediate class (III+IV):merit_effort_dic_mean" =  "Intermediate class (III+IV) x Merit: Effort (BE)",
  "egpService class (I+II):merit_effort_dic_mean" = "Service class (I+II) x Merit: Effort (BE)",
  "egpIntermediate class (III+IV):merit_talent_dic_mean" = "Intermediate class (III+IV) x Merit: Talent (BE)",
  "egpService class (I+II):merit_talent_dic_mean" = "Service class (I+II) x Merit: Talent (BE)")

texreg::screenreg(list(m19,m20,m21,m22),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                groups = list("Social Class (Ref.= Working class (V+VI+VII))" = 4:5,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (WE)" =      10:13,
                              "Social Class (Ref.= Working class (V+VI+VII)) x Meritocracy (BE)" =      14:17),
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("Yes",4))))


## Interactions meritocracy x time without controls

m23 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc +
              merit_effort_mean + merit_talent_mean +
              merit_effort_cwc*ola_num + (1 + ola_num| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m24 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              merit_talent_cwc*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


m25 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              merit_effort_mean*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m26 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              merit_talent_mean*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))


ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2",  
  merit_effort_cwc = "Merit: Effort (WE)",
  merit_talent_cwc = "Merit: Talent (WE)",
  merit_effort_mean = "Merit: Effort (BE)",
  merit_talent_mean = "Merit: Talent (BE)",
  "ola_num:merit_effort_cwc" =  "Time x Merit: Effort (WE)",
  "ola_num:merit_talent_cwc" = "Time x Merit: Talent (WE)",
  "ola_num:merit_effort_mean" = "Time x Merit: Effort (BE)",
  "ola_num:merit_talent_mean" = "Time x Merit: Talent (BE)")

texreg::screenreg(list(m23,m24,m25,m26),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("No",4))))

## Interactions meritocracy x time with controls

m27 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              educ + ideo + sex + age + 
              merit_effort_cwc*ola_num + (1 + ola_num| idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m28 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              educ + ideo + sex + age + 
              merit_talent_cwc*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m29 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              educ + ideo + sex + age + 
              merit_effort_mean*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

m30 <- lmer(just_pension ~ 1 + ola_num + ola_2 + egp + 
              merit_effort_cwc + merit_talent_cwc + 
              merit_effort_mean + merit_talent_mean + 
              educ + ideo + sex + age + 
              merit_talent_mean*ola_num + (1 + ola_num | idencuesta), 
              data = df_study1,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5)))

ccoef <- list(
  "(Intercept)" = "Intercept",
  ola_num = "Wave",
  ola_2 = "Wave^2", 
  merit_effort_cwc = "Merit: Effort (WE)",
  merit_talent_cwc = "Merit: Talent (WE)",
  merit_effort_mean = "Merit: Effort (BE)",
  merit_talent_mean = "Merit: Talent (BE)",
  "ola_num:merit_effort_cwc" =  "Time x Merit: Effort (WE)",
  "ola_num:merit_talent_cwc" = "Time x Merit: Talent (WE)",
  "ola_num:merit_effort_mean" = "Time x Merit: Effort (BE)",
  "ola_num:merit_talent_mean" = "Time x Merit: Talent (BE)")

texreg::screenreg(list(m27,m28,m29,m30),
                caption.above = T,
                caption = NULL,
                stars = c(0.05, 0.01, 0.001),
                custom.coef.map = ccoef,
                digits = 3,
                custom.note = "Note: Cells contain regression coefficients with standard errors in parentheses. %stars.",
                leading.zero = T,
                use.packages = F,
                booktabs = F,
                scalebox = 0.80,
                include.loglik = FALSE,
                include.aic = FALSE,
                center = T,
                custom.gof.rows = list("Controls"=c(rep("Yes",4))))
