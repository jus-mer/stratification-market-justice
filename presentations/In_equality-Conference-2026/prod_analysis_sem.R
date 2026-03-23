# 0. Identification ---------------------------------------------------

# Title: LCA for Merit-scale code of EDUMERCO data
# Institution: JUSMER
# Responsible: Andreas Laffert

# Executive Summary: This script contains the code for run an LCA for merit scale in EDUMERCO data
# Date: March 5, 2026

# 1. Packages  -----------------------------------------------------

if (! require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               sjmisc, 
               here,
               lavaan,
               psych,
               corrplot,
               ggdist,
               sjlabelled,
               patchwork,
               RColorBrewer,
               poLCA,
               reshape,
               summarytools)

options(scipen=999)
rm(list = ls())
# 2. Data -----------------------------------------------------------------

load(url("https://github.com/jus-mer/lca-merit-pensions/raw/refs/heads/main/input/data/proc/db_proc.RData"))

glimpse(db)

## Analytical sample

db <- db %>% 
  mutate(
    income_4 = case_when(
      income %in% c(
        "Menos de $280.000 mensuales liquidos",
        "De $280.001 a $380.000 mensuales liquidos",
        "De $380.001 a $470.000 mensuales liquidos"
      ) ~ "Bajo",
      
      income %in% c(
        "De $470.001 a $610.000 mensuales liquidos",
        "De $610.001 a $730.000 mensuales liquidos"
      ) ~ "Medio-bajo",
      
      income %in% c(
        "De $730.001 a $890.000 mensuales liquidos",
        "De $890.001 a $1.100.000 mensuales liquidos"
      ) ~ "Medio-alto",
      
      income %in% c(
        "De $1.100.001 a $2.700.000 mensuales liquidos",
        "De $2.700.001 a $4.100.000 mensuales liquidos",
        "Mas de $4.100.001 mensuales liquidos"
      ) ~ "Alto",
      
      TRUE ~ NA_character_
    ),
    income_4 = factor(income_4, levels = c("Bajo", "Medio-bajo", "Medio-alto", "Alto"))
  )

db$sex <- if_else(db$sex == 1, "Male", "Female")
db$sex <- factor(db$sex, levels = c("Male", "Female"))

db <- db %>% 
  dplyr::select(-income) %>% 
  na.omit()
# 3. Analysis -------------------------------------------------------------

# 3.1 descriptive ----
vars_m <- c("perc_effort",
            "perc_talent",
            "perc_rich_parents",
            "perc_contact",
            "pref_effort",
            "pref_talent",
            "pref_rich_parents",
            "pref_contact")

t1 <- db %>% 
  dplyr::select(all_of(vars_m), ends_with("_d"))


df<-dfSummary(t1,
              plain.ascii = FALSE,
              style = "multiline",
              tmp.img.dir = "/tmp",
              graph.magnif = 0.75,
              headings = F,  # encabezado
              varnumbers = F, # num variable
              labels.col = T, # etiquetas
              na.col = T,    # missing
              graph.col = T, # plot
              valid.col = T, # n valido
              col.widths = c(20,10,10,10,10,10))

df$Variable <- NULL # delete variable column

print(df, method="render")

db <- db %>% 
  mutate(
    across(
      .cols = all_of(vars_m),
      .fns = ~as.numeric(.)
    ))

labels1 <- c("Strongly desagree" = 1, 
             "Desagree" = 2, 
             "Agree" = 3, 
             "Strongly agree" = 4)
db <- db %>% 
  mutate(
    across(
      .cols = all_of(vars_m),
      .fns = ~  sjlabelled::set_labels(., labels = labels1)
    )
  )

df <- db %>% 
  dplyr::select(all_of(vars_m)) %>% 
  drop_na()

theme_set(theme_ggdist())
colors <- RColorBrewer::brewer.pal(n = 4, name = "RdBu")


a <- df %>% 
  dplyr::select(perc_effort, 
                perc_talent,
                perc_rich_parents,
                perc_contact) %>% 
  sjPlot::plot_likert(geom.colors = colors,
                      title = c("a. Perceptions"),
                      geom.size = 0.8,
                      axis.labels = c("Effort", "Talent", "Rich parents", "Contacts"),
                      catcount = 4,
                      values  =  "sum.outside",
                      reverse.colors = F,
                      reverse.scale = T,
                      show.n = FALSE,
                      show.prc.sign = T
  ) +
  ggplot2::theme(legend.position = "none",
                 text = element_text(size = 16))

b <- df %>% 
  dplyr::select(pref_effort, 
                pref_talent,
                pref_rich_parents,
                pref_contact) %>% 
  sjPlot::plot_likert(geom.colors = colors,
                      title = c("b. Preferences"),
                      geom.size = 0.8,
                      axis.labels = c("Effort", "Talent", "Rich parents", "Contacts"),
                      catcount = 4,
                      values  =  "sum.outside",
                      reverse.colors = F,
                      reverse.scale = T,
                      show.n = FALSE,
                      show.prc.sign = T
  ) +
  ggplot2::theme(legend.position = "bottom",
                 text = element_text(size = 16))

likerplot <- a / b + plot_annotation(caption = paste0("Source: own elaboration based on Survey EDUMERCO"," (n = ",dim(df)[1],")"
))

likerplot

# 3.2 correlations ----

M <- df %>% 
  psych::polychoric()

diag(M$rho) <- NA

rownames(M$rho) <- c("A. Perception Effort",
                     "B. Perception Talent",
                     "C. Perception Rich parents",
                     "D. Perception Contacts",
                     "E. Preference Effort",
                     "F. Preference Talent",
                     "G. Preference Rich parents",
                     "H. Preference Contacts")

#set Column names of the matrix
colnames(M$rho) <-c("(A)", "(B)","(C)","(D)","(E)","(F)","(G)",
                    "(H)")

testp <- cor.mtest(M$rho, conf.level = 0.95)

#Plot the matrix using corrplot
corrplot::corrplot(M$rho,
                   method = "color",
                   addCoef.col = "black",
                   type = "upper",
                   tl.col = "black",
                   col = colorRampPalette(c("#E16462", "white", "#0D0887"))(12),
                   bg = "white",
                   na.label = "-") 

# 3.3 CFA: ----
#### CFA All countries #### 

model_base <- ('
perc_merit =~ perc_effort + perc_talent
perc_nmerit =~ perc_rich_parents + perc_contact
pref_merit =~ pref_effort + pref_talent
pref_nmerit =~ pref_rich_parents + pref_contact
')

# Estimación 
db %>% 
  dplyr::select(all_of(vars_m)) %>% 
  mardia(na.rm = TRUE, plot=TRUE)

fit_cfa <<- cfa(model = model_base, 
                data = db, 
                estimator = "WLSMV",
                ordered = T,
                std.lv = F,
                parameterization = "theta")

summary(fit_cfa, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

fitmeasures(fit_cfa, c("chisq", "pvalue", "df", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))

# 3.4 SEM ----

db_sem <- db %>% 
  dplyr::select(just_pension, all_of(vars_m), age, sex, educ, income_4, pol)

db_sem$just_pension <- as.numeric(db$just_pension)

db_sem <- db_sem %>% 
  mutate(
    across(
      .cols = all_of(vars_m),
      .fns = ~as.numeric(.)
    ))

# asegurar referencias
db_sem$income_4 <- relevel(db_sem$income_4, ref = "Bajo")
db_sem$pol <- relevel(db_sem$pol, ref = "Left")

# crear dummies
X_income <- model.matrix(~ income_4, data = db_sem)[, -1, drop = FALSE]
X_pol    <- model.matrix(~ pol, data = db_sem)[, -1, drop = FALSE]

# unir
db_sem <- cbind(db_sem, as.data.frame(X_income), as.data.frame(X_pol))

# limpiar nombres
names(db_sem) <- make.names(names(db_sem))

db_sem$sex_female <- ifelse(db_sem$sex == "Female", 1, 0)

model <- c('
  perc_merit =~ perc_effort + perc_talent
  perc_nmerit =~ perc_rich_parents + perc_contact
  pref_merit =~ pref_effort + pref_talent
  pref_nmerit =~ pref_rich_parents + pref_contact

  just_pension ~ perc_merit + perc_nmerit + pref_merit + pref_nmerit +
                 age + educ + sex_female +
                 income_4Medio.bajo + income_4Medio.alto + income_4Alto +
                 polCenter + polRight + polDoes.not.identify
')

ord_vars <- c(
  "just_pension",
  "perc_effort", "perc_talent",
  "perc_rich_parents", "perc_contact",
  "pref_effort", "pref_talent",
  "pref_rich_parents", "pref_contact"
)

fit_sem <- lavaan::sem(
  model,
  data = db_sem,
  estimator = "WLSMV",
  ordered = ord_vars
)

summary(fit_sem, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

fitmeasures(fit_sem, c("chisq", "pvalue", "df", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))


