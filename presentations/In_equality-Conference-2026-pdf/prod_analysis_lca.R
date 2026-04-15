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
  ggplot2::theme(legend.position = "none")

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
  ggplot2::theme(legend.position = "bottom")

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

# 3.3 LCA: dichotomized items ----

# --- helpers -------------------------------------------------------------

get_K <- function(fit) {
  # robust: number of classes = columns of posterior
  if (!is.null(fit$posterior)) return(as.integer(ncol(fit$posterior)))
  # fallback: rows of probs for first item
  if (!is.null(fit$probs) && length(fit$probs) > 0) return(as.integer(nrow(fit$probs[[1]])))
  stop("Cannot determine K from poLCA object (no posterior/probs).")
}

get_npar <- function(fit) {
  if (!is.null(fit$npar)) return(as.numeric(fit$npar))
  # fallback if npar missing (rare): compute from probs + class shares
  K <- get_K(fit)
  Rj <- vapply(fit$probs, function(m) if (is.matrix(m)) ncol(m) else NA_integer_, integer(1))
  if (anyNA(Rj)) stop("Cannot compute npar: missing category counts in fit$probs.")
  as.numeric(sum(K * (Rj - 1)) + (K - 1))
}

get_n_patterns_obs <- function(fit) {
  # poLCA usually stores predicted cell counts for observed patterns as a vector
  if (!is.null(fit$predcell)) return(as.integer(length(fit$predcell)))
  # fallback: sometimes observed pattern table exists
  if (!is.null(fit$observed) && is.matrix(fit$observed)) return(as.integer(nrow(fit$observed)))
  NA_integer_
}

n_patterns_from_data <- function(data_indicators) {
  # data_indicators = db_lca with only indicators (your a1..a8), already NA-omitted
  as.integer(nrow(dplyr::distinct(data_indicators)))
}


entropy_stats <- function(fit) {
  P <- fit$posterior
  P <- pmax(P, 1e-12)
  N <- nrow(P); K <- ncol(P)
  
  # relative entropy (0..1): 1 = best separation
  rel_entropy <- 1 - (-sum(P * log(P)) / (N * log(K)))
  
  class_hat <- max.col(P)
  max_p <- apply(P, 1, max)
  
  avg_diag <- vapply(1:K, function(k) {
    idx <- which(class_hat == k)
    if (length(idx) == 0) NA_real_ else mean(P[idx, k])
  }, numeric(1))
  
  tibble(
    entropy = rel_entropy,
    mean_max_p = mean(max_p),
    p10_max_p = unname(quantile(max_p, 0.10)),
    avg_diag_min = suppressWarnings(min(avg_diag, na.rm = TRUE)),
    avg_diag_mean = mean(avg_diag, na.rm = TRUE)
  )
}

min_class_share <- function(fit) {
  Ppost <- fit$posterior
  N <- nrow(Ppost); K <- ncol(Ppost)
  class_hat <- max.col(Ppost)
  
  tibble(
    min_class_prior = min(fit$P),
    min_class_map   = min(tabulate(class_hat, nbins = K)) / N
  )
}

# 3.3 LCA: dichotomized items ----

df <- db %>% 
  dplyr::select(id, ends_with("_d")) %>% 
  na.omit()

a1 <- df$perc_effort_d
a2 <- df$perc_talent_d
a3 <- df$perc_rich_parents_d
a4 <- df$perc_contact_d
a5 <- df$pref_effort_d
a6 <- df$pref_talent_d 
a7 <- df$pref_rich_parents_d
a8 <- df$pref_contact_d
id <- df$id

db_lca <- data.frame(id, a1, a2, a3, a4, a5, a6, a7, a8) %>% as_tibble()

db_lca <- db_lca %>% 
  mutate_at(.vars = 2:9, ~as.factor(.))

f <-cbind(a1, a2, a3, a4, a5, a6, a7, a8)~1

fit_lca <- function(K, data, formula, nrep = 50, maxiter = 3000, seed = 123) {
  set.seed(seed)
  poLCA(formula, data = data, nclass = K,
        nrep = nrep, maxiter = maxiter,
        verbose = FALSE, graphs = FALSE)
}

Ks <- 1:5
#fits <- lapply(Ks, fit_lca, data = db_lca, formula = f, nrep = 30)
#save(fits, file = here("output/fits_models_lca.RData"))

load(url("https://github.com/jus-mer/lca-merit-pensions/raw/refs/heads/main/output/fits_models_lca.RData"))

fit_tbl <- data.frame(
  K    = Ks,
  logLik = sapply(fits, `[[`, "llik"),
  AIC  = sapply(fits, `[[`, "aic"),
  BIC  = sapply(fits, `[[`, "bic"),
  Gsq  = sapply(fits, `[[`, "Gsq"),
  Chisq = sapply(fits, `[[`, "Chisq")
)

fit_tbl

# ---- build table (no list-cols) -----------------------------------------

N_used <- nrow(db_lca)

m_patterns <- n_patterns_from_data(db_lca)  # scalar

# Robust extractors (avoid list-columns)
get_K <- function(fit) as.integer(ncol(fit$posterior))  # classes = cols of posterior
get_ll  <- function(fit) as.numeric(fit$llik)
get_aic <- function(fit) as.numeric(fit$aic)
get_bic <- function(fit) as.numeric(fit$bic)
get_npar <- function(fit) {
  if (!is.null(fit$npar)) return(as.numeric(fit$npar))
  # fallback: compute from probs
  K <- get_K(fit)
  Rj <- vapply(fit$probs, function(m) ncol(m), integer(1))
  as.numeric(sum(K * (Rj - 1)) + (K - 1))
}

# Table of fit indices
fit_tbl <- tibble(
  K     = map_int(fits, get_K),
  N     = nrow(db_lca),
  logLik= map_dbl(fits, get_ll),
  npar  = map_dbl(fits, get_npar),
  AIC   = map_dbl(fits, get_aic),
  BIC   = map_dbl(fits, get_bic)
) %>%
  arrange(K) %>%
  mutate(
    # Differences relative to preceding (K-1) model
    dAIC = AIC - lag(AIC),
    dBIC = BIC - lag(BIC),
    
    # Likelihood-ratio test (K vs K-1)
    # Test statistic: 2*(LL_K - LL_{K-1}) ~ Chi-square with df = npar_K - npar_{K-1}
    # Note: in LCA, models are not strictly nested in the regular sense;
    # this LRT is widely used as a heuristic, but interpret with caution.
    LRT = 2 * (logLik - lag(logLik)),
    df_LRT = npar - lag(npar),
    p_LRT = ifelse(!is.na(LRT) & df_LRT > 0, pchisq(LRT, df = df_LRT, lower.tail = FALSE), NA_real_)
  )


fit_tbl_report <- fit_tbl %>%
  transmute(
    K, N,
    logLik = round(logLik, 1),
    npar = npar,
    AIC = round(AIC, 1),
    dAIC = round(dAIC, 1),
    BIC = round(BIC, 1),
    dBIC = round(dBIC, 1),
    LRT = round(LRT, 1),
    df_LRT = df_LRT,
    p_LRT = signif(p_LRT, 3)
  )

fit_tbl_report %>% 
  kableExtra::kable(format = "html",
                    align = "c",
                    booktabs = T,
                    escape = F,
                    caption = NULL) %>%
  kableExtra::kable_styling(full_width = T,
                            latex_options = "hold_position",
                            bootstrap_options=c("striped", "bordered", "condensed"),
                            font_size = 23) %>%
  kableExtra::column_spec(c(1,8), width = "3.5cm") %>% 
  kableExtra::column_spec(2:7, width = "4cm") %>% 
  kableExtra::column_spec(4, width = "5cm")

# Visual: AIC/BIC by K (lower is better)
fit_long <- fit_tbl %>%
  dplyr::select(K, AIC, BIC) %>%
  pivot_longer(-K, names_to = "criterion", values_to = "value")

ggplot(fit_long, aes(x = K, y = value, group = criterion)) +
  geom_line() +
  geom_point() +
  facet_wrap(~criterion, scales = "free_y") +
  labs(x = "Number of classes (K)", y = "Value") 

# Visual: ΔAIC and ΔBIC (negative = improvement over K-1)
delta_long <- fit_tbl %>%
  dplyr::select(K, dAIC, dBIC) %>%
  pivot_longer(-K, names_to = "delta", values_to = "value")

ggplot(delta_long, aes(x = K, y = value, group = delta)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  facet_wrap(~delta, scales = "free_y") +
  labs(x = "K", y = "Change vs K-1 (negative = better)") 


# 4. Four latent class model ----------------------------------------------

fit4 <- fits[[4]]

# 2) Extract conditional probabilities per item and class
# fit4$probs is a list: each element is K x Rj matrix.
# For dichotomous items: K x 2 (categories = colnames or "1","2")
probs_long <- bind_rows(lapply(names(fit4$probs), function(item) {
  mat <- fit4$probs[[item]]
  df <- as.data.frame(mat)
  df$class <- 1:nrow(mat)
  df$item <- item
  df
})) %>%
  pivot_longer(cols = -c(class, item), names_to = "category", values_to = "prob") %>%
  mutate(
    class = as.integer(class),
    prob = as.numeric(prob)
  )

# If levels are c("0","1") or c("1","2") you'll want the "higher"/"yes" one.
# Here I will assume endorsement is the LAST level (most common after dichotomization).
endorsement_level <- tail(levels(db_lca$a1), 1)

# 4) Build a clean table: P(endorsement | class) for each item
p_endorse <- probs_long %>%
  filter(category == "Pr(2)") %>%
  dplyr::select(item, class, p = prob) %>%
  mutate(item = factor(item, levels = names(fit4$probs)))

# 5) Plot: profile lines (one line per class)
p_endorse <- p_endorse %>% 
  mutate(
    variable = case_when(item == "a1" ~ "perc_effort",
                         item == "a2" ~ "perc_talent",
                         item == "a3" ~ "perc_rich_parents",
                         item == "a4" ~ "perc_contacts",
                         item == "a5" ~ "pref_effort",
                         item == "a6" ~ "pref_talent",
                         item == "a7" ~ "pref_rich_parents",
                         item == "a8" ~ "pref_contacts"
    ),
    variable = factor(variable,
                      levels = c("perc_effort",
                                 "perc_talent",
                                 "perc_rich_parents",
                                 "perc_contacts",
                                 "pref_effort",
                                 "pref_talent",
                                 "pref_rich_parents",
                                 "pref_contacts")))

ggplot(p_endorse, aes(x = variable, y = p, group = factor(class), colour = factor(class))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Item", y = paste0("P(Y = ", "High", " | class)"), group = "Class", colour = "Class") +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

# 6) Plot: heatmap (good for papers)
ggplot(p_endorse, aes(x = variable, y = factor(class), fill = p)) +
  geom_tile() +
  coord_cartesian(xlim = NULL) +
  labs(x = "Item", y = "Class", fill = paste0("P(Y=", "High", ")")) +
  theme(legend.position = "right",
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

# 7) Report table (wide): items x classes with probabilities
class_prev <- tibble(
  class = 1:length(fit4$P),
  pi = as.numeric(fit4$P)
) %>%
  mutate(pi_pct = 100 * pi)

report_probs <- p_endorse %>%
  mutate(class = paste0("Class_", class)) %>%
  pivot_wider(names_from = class, values_from = p) %>%
  mutate(across(starts_with("Class_"), ~ round(.x, 3)))

prev_lab <- class_prev %>%
  transmute(class = paste0("Class_", class),
            lab = paste0(class, " (", round(pi_pct, 1), "%)"))

report_probs2 <- report_probs
for (i in seq_len(nrow(prev_lab))) {
  old <- prev_lab$class[i]
  new <- prev_lab$lab[i]
  names(report_probs2)[names(report_probs2) == old] <- new
}

report_probs2[,-c(1)] %>% 
  kableExtra::kable(format = "html",
                    align = "c",
                    booktabs = T,
                    escape = F,
                    caption = NULL) %>%
  kableExtra::kable_styling(full_width = T,
                            latex_options = "hold_position",
                            bootstrap_options=c("striped", "bordered", "condensed"),
                            font_size = 23)


# 5. Regression analysis --------------------------------------------------

fit4 <- fits[[4]]

post <- as.data.frame(fit4$posterior)
names(post) <- paste0("p_class", 1:4)

class_df <- na.omit(db_lca) %>%
  dplyr::select(id) %>%
  bind_cols(post) %>%
  mutate(
    class_hat = max.col(dplyr::select(., starts_with("p_class"))),
    class_hat = factor(class_hat)
  )

db_full <- db %>%
  left_join(class_df, by = "id")

db_full$just_pension <- as.numeric(db_full$just_pension)
db_full$just_pension <- if_else(db_full$just_pension >= 3, 1, 0)

P <- fits[[4]]$posterior
max_p <- apply(P, 1, max)
#summary(max_p)
#quantile(max_p, c(.1, .25, .5, .75, .9))

# MAP (ref=2)
db_reg1 <- db_full %>% filter(!is.na(class_hat), !is.na(just_pension), !is.na(age), !is.na(sex), !is.na(educ), !is.na(income), !is.na(pol))

db_reg1 <- db_reg1 %>% 
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

db_reg1$class_hat <- relevel(db_reg1$class_hat, ref = "4")

m1 <- lm(just_pension ~ class_hat, data = db_reg1)
m2 <- lm(just_pension ~ class_hat + age, data = db_reg1)
m3 <- lm(just_pension ~ class_hat + age + sex, data = db_reg1)
m4 <- lm(just_pension ~ class_hat + age + sex + educ, data = db_reg1)
m5 <- lm(just_pension ~ class_hat + age + sex + educ + income_4, data = db_reg1)
m6 <- lm(just_pension ~ class_hat + age + sex + educ + income_4 + pol, data = db_reg1)
m7 <- lm(just_pension ~ class_hat*income_4 + age + sex + educ + pol, data = db_reg1)

# Posterior (ref implícita = clase 2)
#db_reg2 <- db_full %>% filter(!is.na(p_class2), #!is.na(just_pension))
#m_post <- lm(just_pension ~ p_class1 + p_class3 + p_class4, #data = db_reg2)
library(texreg)
screenreg(list(m1, m2, m2, m3, m4, m5,m6))
screenreg(list(m7))

coef_df <-
  broom::tidy(m1, conf.int = TRUE) %>%
  mutate(model = "MAP (dummies; ref=Class 4)") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    contrast = case_when(
      term == "class_hat1" ~ "Class 1 vs 4",
      term == "class_hat2" ~ "Class 2 vs 4",
      term == "class_hat3" ~ "Class 3 vs 4",
      TRUE ~ term
    )) %>%
  mutate(
    contrast = factor(contrast, levels = c("Class 1 vs 4", "Class 2 vs 4", "Class 3 vs 4"))
  )

ggplot(coef_df, aes(x = estimate, y = contrast)) +
  geom_vline(xintercept = 0, linewidth = 0.7, color = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), fatten = 2, 
                  size = 1, color = "#ca1137") + 
  labs(x = "Coefficient (OLS) with 95% CI", y = NULL) +
  theme(text = element_text(size = 14))

coef_df <-
  broom::tidy(m6, conf.int = TRUE) %>%
  mutate(model = "MAP (dummies; ref=Class 4)") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    contrast = case_when(
      term == "class_hat1" ~ "Class 1 vs 4",
      term == "class_hat2" ~ "Class 2 vs 4",
      term == "class_hat3" ~ "Class 3 vs 4",
      
      term == "sex2" ~ "Female (Ref.= Male)",
      term == "educ" ~ "Education (in years)",
      term == "age" ~ "Age",
      term == "income_4Medio-bajo" ~ "Income Middle-Low",
      term == "income_4Medio-alto" ~ "Income Middle-High",
      term == "income_4Alto" ~ "Income High",   
      term == "polCenter"  ~"Center",         
      term == "polRight"  ~"Right",
      term == "polDoes not identify" ~ "Does not identify",
      TRUE ~ term
    )) 

ggplot(coef_df, aes(x = estimate, y = contrast)) +
  geom_vline(xintercept = 0, linewidth = 0.7, color = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), fatten = 2, 
                  size = 1, color = "#ca1137") + 
  labs(x = "Coefficient (OLS) with 95% CI", y = NULL) +
  theme(text = element_text(size = 14))

coef_df <-
  broom::tidy(m6, conf.int = TRUE) %>%
  slice_tail(n=9) %>%
  mutate(
    contrast = case_when(
      term == "class_hat1" ~ "Class 1 vs 4",
      term == "class_hat2" ~ "Class 2 vs 4",
      term == "class_hat3" ~ "Class 3 vs 4",
      term == "sex2" ~ "Female (Ref.= Male)",
      term == "educ" ~ "Education (in years)",
      term == "age" ~ "Age",
      term == "income_4Medio-bajo" ~ "Income Middle-Low",
      term == "income_4Medio-alto" ~ "Income Middle-High",
      term == "income_4Alto" ~ "Income High",           
      TRUE ~ term
    )) 

ggplot(coef_df, aes(x = estimate, y = contrast)) +
  geom_vline(xintercept = 0, linewidth = 0.7, color = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), fatten = 2, 
                  size = 1, color = "#ca1137") + 
  labs(x = "Coefficient (OLS) with 95% CI", y = NULL) +
  theme(text = element_text(size = 14))
