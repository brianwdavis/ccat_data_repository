sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices datasets  utils    
# [6] methods   base     
# 
# other attached packages:
# ## omitted, see `renv.lock`
#
# loaded via a namespace (and not attached):
# ## omitted, see `renv.lock`


# Optionally: load the appropriate package versions from `renv.lock`
# renv::restore()

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(lubridate)
library(forcats)
library(stringr)
library(purrr)

library(lme4)
library(lmerTest)


## helper functions ----
log_mean_se <- function(x, mult = 1) 
{
  x <- stats::na.omit(log(x))
  x <- x[is.finite(x)]
  se <- mult * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  data.frame(y = exp(mean), ymin = exp(mean - se), ymax = exp(mean + se))
}

LSD_extract_lmer <- function(mod) {
  DFerror <- nobs(mod) - length(fixef(mod)) - 1
  MSerror <- mean((resid(mod))^2, na.rm = T)
  
  alpha <- 0.05
  Tprob <- qt(1 - alpha/2, DFerror)
  
  nr <- augment(mod) %>% 
    group_by(.fixed) %>% 
    tally() %>% 
    pull(n)
  mean_recip_nr <- mean(1/nr)
  
  LSD <- Tprob * sqrt(2 * MSerror*mean_recip_nr)
  return(LSD)
}

pval_to_star <- function(x) {
  cut(
    x, 
    breaks = c(-1, 0.001, 0.01, 0.05, 1),
    labels = c("***", "**", "*", "")
  ) %>% 
    as.character()
}

## custom theme ----
ccat_theme <- 
  theme(
    panel.background = element_rect(colour = "black"),
    strip.background = element_rect(colour = NA, fill = NA),
    legend.background = element_blank(),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    legend.direction = "horizontal",
    legend.key.width = unit(0.5, "lines"),
    strip.text = element_text(size = rel(1.05), face = "bold"),
    axis.text.x = element_text(size = rel(1.1), angle = 90, vjust = 0.5),
    panel.spacing = unit(0.25, "lines")
  )

# Read in data ----
# _ biomass ----
ccat_biomass <- readr::read_csv("biomass.csv")

# _ time to flowering ----
ccat_gddto25 <- readr::read_csv("time_to_flowering.csv")

# _ N and BNF ----
ccat_n <- readr::read_csv("n_content.csv")

# _ survival/emergence ----
ccat_survival <- readr::read_csv("survival_emergence.csv")


## biomass model ----
ccat_biomass_subset <-
  ccat_biomass %>% 
  filter(bm_per_plant < 50) %>%             # [*]
  group_by(Year, `PI`) %>% 
  mutate(n = length(bm_per_plant)) %>% 
  ungroup() %>% 
  filter(n > 1) %>% 
  mutate(yrrep = paste(Year, Rep), Year = as.factor(Year))

fit_ccat_bm <- lmer(
  log(bm_per_plant) ~ PI * Year + 0 + (1|yrrep), 
  data = ccat_biomass_subset
  )

# [*] Remove one outlier >4x any other observation.


## biomass plot ----
augment(fit_ccat_bm) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, exp(.fixed + .resid))) +
  stat_summary(
    fun.data = log_mean_se,
    aes(y = exp(.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = log_mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  labs(
    x = NULL,
    y = expression("Biomass at harvest, g dry wt plant" ^ -1),
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("biomass.png", width = 6, height = 4.5, dpi = 300)

## biomass plot alternate ----
augment(fit_ccat_bm) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, exp(.fixed + .resid))) +
  stat_summary(
    fun.data = log_mean_se,
    aes(y = exp(.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = log_mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 31,
    ymin = exp(.8 - LSD_extract_lmer(fit_ccat_bm) / 2),
    ymax = exp(.8 + LSD_extract_lmer(fit_ccat_bm) / 2)
  ) +
  annotate(
    geom = "text",
    x = 31,
    y = exp(.8 + 1.25 * LSD_extract_lmer(fit_ccat_bm) / 2),
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  expand_limits(x = c(1, 30)) +
  labs(
    x = NULL,
    y = expression("Biomass at harvest, g dry wt plant" ^ -1),
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(breaks = outer(0.25 * 1:4, 10 ^ (-1:1))) +
  coord_trans(y = "log10", limy = c(0.2, 17)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("biomass logscaled.png", width = 6, height = 4.5, dpi = 300)


### GDD model ----
ccat_gddto25_subset <- 
  ccat_gddto25 %>% 
  group_by(Year, PI) %>% 
  mutate(n = length(PI)) %>% 
  ungroup() %>% 
  filter(n > 1) %>% 
  mutate(yrrep = paste(Year, Rep), Year = as.factor(Year))

fit_ccat_gdd <- 
  lmer(
    gddto25 ~ PI * Year + 0 + (1|yrrep),
    data = ccat_gddto25_subset
    )


### GDD plot ----
augment(fit_ccat_gdd) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, (.fixed + .resid))) +
  stat_summary(
    fun.data = mean_se,
    aes(y = (.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = mean_se, size = 0.75,
    aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 30,
    ymin = 1700 - LSD_extract_lmer(fit_ccat_gdd) / 2,
    ymax = 1700 + LSD_extract_lmer(fit_ccat_gdd) / 2
  ) +
  annotate(
    geom = "text",
    x = 30,
    y = 1700 + 1.25 * LSD_extract_lmer(fit_ccat_gdd) / 2,
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  labs(
    x = NULL,
    y = expression("GDD"[0 * degree * C] * " to reach 25% maturity"),
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(limits = c(1400, NA)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("gdd.png", width = 6, height = 4.5, dpi = 300)




# %N model ----
fit_ccat_N <- lmer(
  log(pctN) ~ PI*Year+0 + (1|yrrep), 
  data = ccat_n %>% mutate(
    yrrep = paste(Year, Rep), 
    Year = as.factor(Year)
    )
  )

# %N plot ----
augment(fit_ccat_N) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, exp(.fixed + .resid))) +
  stat_summary(
    fun.data = log_mean_se,
    aes(y = exp(.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = log_mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  labs(
    x = NULL,
    y = "Plant N content, %",
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(limits = c(1.25, 4), breaks = seq(1, 4, by = .5)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("pctN.png", width = 6, height = 4.5, dpi = 300)

# %N plot alternate ----
augment(fit_ccat_N) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, exp(.fixed + .resid))) +
  stat_summary(
    fun.data = log_mean_se,
    aes(y = exp(.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = log_mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 30,
    ymin = exp(.8 - LSD_extract_lmer(fit_ccat_N) / 2),
    ymax = exp(.8 + LSD_extract_lmer(fit_ccat_N) / 2)
  ) +
  annotate(
    geom = "text",
    x = 30,
    y = exp(.8 + 1.25 * LSD_extract_lmer(fit_ccat_N) / 2),
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  expand_limits(x = c(1, 30)) +
  labs(
    x = NULL,
    y = "Plant N content, %",
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(breaks = seq(1, 4, by = .5)) +
  coord_trans(y = "log10",
              limy = scales::expand_range(c(1.5, 4.25), .05)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("pctN logscaled.png", width = 6, height = 4.5, dpi = 300)


# from BNF model ----
fit_ccat_bnf <- lmer(
  pctfromBNF ~ PI*Year + 0 + (1|yrrep), 
  data = ccat_n %>% mutate(
    yrrep = paste(Year, Rep),
    Year = as.factor(Year)
    )
)

# from BNF plot ----
augment(fit_ccat_bnf) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  ggplot(aes(PI, (.fixed + .resid))) +
  stat_summary(fun.data = mean_se,
               aes(y = (.fitted + .resid),
                   shape = factor(Year)),
               color = "grey50") +
  stat_summary(
    fun.data = mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 30,
    ymin = 65 - LSD_extract_lmer(fit_ccat_bnf) / 2,
    ymax = 65 + LSD_extract_lmer(fit_ccat_bnf) / 2
  ) +
  annotate(
    geom = "text",
    x = 30,
    y = 65 + 1.25 * LSD_extract_lmer(fit_ccat_bnf) / 2,
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  expand_limits(x = c(1, 31)) +
  labs(
    x = NULL,
    y = "Proportion of plant N content as BNF, %",
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  coord_cartesian(ylim = c(0, 100)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme 
ggsave("pctNfromBNF.png", width = 6, height = 4.5, dpi = 300)



# Survival model ----

fit_ccat_survival <-
  lmer(
    survival ~ PI * Year+0 + (1|yrrep), 
    data = ccat_survival %>% mutate(
      yrrep = paste(Year, Rep),
      Year = as.factor(Year)
      )
    )

# Survival plot ----
augment(fit_ccat_survival) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  group_by(PI) %>%
  ggplot(aes(PI, (.fixed + .resid))) +
  stat_summary(
    fun.data = mean_se,
    aes(y = (.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 31,
    ymin = 0.5 - LSD_extract_lmer(fit_ccat_survival) / 2,
    ymax = 0.5 + LSD_extract_lmer(fit_ccat_survival) / 2
  ) +
  annotate(
    geom = "text",
    x = 31,
    y = 0.5 + 1.25 * LSD_extract_lmer(fit_ccat_survival) / 2,
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  labs(
    x = NULL,
    y = "Winter survival, % of fall-emerged plants",
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expand_scale(mult = 0.025),
    breaks = 0:4 / 4
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  expand_limits(x = c(1, 31)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme
ggsave("survival.png", width = 6, height = 4.5, dpi = 300)


# Emergence model ----
fit_ccat_emergence <-
  lmer(
    emergence ~ PI*Year+0 + (1|yrrep), 
    data = ccat_survival %>% mutate(
      yrrep = paste(Year, Rep),
      Year = as.factor(Year))
    )

# Emergence plot ----
augment(fit_ccat_emergence) %>%
  group_by(PI) %>%
  mutate(grand = mean(.fixed)) %>%
  ungroup() %>%
  mutate(rnk = dense_rank(-grand)) %>%
  arrange(rnk) %>%
  mutate(PI = fct_inorder(PI)) %>%
  group_by(PI) %>%
  ggplot(aes(PI, (.fixed + .resid))) +
  stat_summary(
    fun.data = mean_se,
    aes(y = (.fitted + .resid),
        shape = factor(Year)),
    color = "grey50"
  ) +
  stat_summary(
    fun.data = mean_se, size = 0.75, aes(shape = "Grand mean")
  ) +
  annotate(
    geom = "errorbar",
    x = 30,
    ymin = 0.5 - LSD_extract_lmer(fit_ccat_emergence) / 2,
    ymax = 0.5 + LSD_extract_lmer(fit_ccat_emergence) / 2
  ) +
  annotate(
    geom = "text",
    x = 30,
    y = 0.5 + 1.25 * LSD_extract_lmer(fit_ccat_emergence) / 2,
    label = "bold(LSD[0.05])",
    hjust = 0,
    parse = T,
    angle = 90
  ) +
  labs(
    x = NULL,
    y = "Fall emergence, % of seeds planted",
    shape = "Harvest Year"
  ) +
  scale_shape_manual(values = c(25, 22, 24, 21)) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expand_scale(mult = 0.025)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  theme_classic() +
  ccat_theme
ggsave("emergence.png", width = 6, height = 4.5, dpi = 300)


# Tables of estimates ----

ccat_summarizer <- function(model, name, backtrans = identity) {
  
  LCL = paste0(name, "_LCL")
  UCL = paste0(name, "_UCL")
  
  augment(model) %>% 
    group_by(PI) %>% 
    summarise(
      !!name := mean(.fixed),
      !!LCL  := !!sym(name) - 1.96*sd(.resid)/sqrt(length(.resid)),
      !!UCL  := !!sym(name) + 1.96*sd(.resid)/sqrt(length(.resid))
    ) %>% 
    mutate_if(is.numeric, backtrans) %>%
    mutate_if(is.numeric, signif, 3)
}

ccat_summary_labeller <- function(dat) {
  dat[[2]] <- glue::glue("{dat[[2]]} ({dat[[3]]}, {dat[[4]]})")
  dat[,1:2]
}

ccat_table <- 
  tribble(
    ~mods,             ~nms,                  ~btrans,
    fit_ccat_emergence,"emergence_pct",        identity,
    fit_ccat_survival, "survival_pct",        identity,
    fit_ccat_bm,       "bio_g_plant",         exp, 
    fit_ccat_gdd,      "gdd_to_25pct_flower", identity, 
    fit_ccat_N,        "pctN",                exp, 
    fit_ccat_bnf,      "pctBNF",              identity
  ) %>% 
  mutate(d = pmap(list(mods, nms, btrans), ccat_summarizer)) %>% 
  pull(d) %>% 
  map(ccat_summary_labeller) %>% 
  reduce(full_join) 
  
write.csv(ccat_table, "CCAT estimates.csv", row.names = F)

ccat_lsd <- 
tibble(
  measurement = c("Emergence", "Survival", "Biomass", "GDD to 25%", "%N", "%BNF"),
  type = c("LSD", "LSD", "LSR", "LSD", "LSR", "LSD"),
  statistic = c(
    LSD_extract_lmer(fit_ccat_emergence),
    LSD_extract_lmer(fit_ccat_survival),
    exp(LSD_extract_lmer(fit_ccat_bm)),
    LSD_extract_lmer(fit_ccat_gdd),
    exp(LSD_extract_lmer(fit_ccat_N)),
    LSD_extract_lmer(fit_ccat_bnf)
  ),
  GrandMean = c(
    augment(fit_ccat_emergence) %>% pull(.fixed) %>% mean(),
    augment(fit_ccat_survival) %>% pull(.fixed) %>% mean(),
    augment(fit_ccat_bm) %>% pull(.fixed) %>% mean() %>% exp(),
    augment(fit_ccat_gdd) %>% pull(.fixed) %>% mean(),
    augment(fit_ccat_N) %>% pull(.fixed) %>% mean() %>% exp(),
    augment(fit_ccat_bnf) %>% pull(.fixed) %>% mean()
  )
) %>% 
  mutate_if(is.numeric, signif, 4)

write.csv(ccat_lsd, "CCAT LSDs.csv", row.names = F)


# Model summaries ----
fit_ccat_bm_null <- 
  lmer(
    log(bm_per_plant) ~ 1 + (1|yrrep), 
    data = ccat_biomass_subset
    )

fit_ccat_bm_null_unrandom <- 
  lm(log(bm_per_plant) ~ 1, data = ccat_biomass_subset)


fit_ccat_gdd_null <- 
  lmer(gddto25 ~ 1 + (1|yrrep), data = ccat_gddto25_subset)

fit_ccat_gdd_null_unrandom <- 
  lm(gddto25 ~ 1, data = ccat_gddto25_subset)


fit_ccat_N_null <- lmer(
  log(pctN) ~ 1 + (1|yrrep), 
  data = ccat_n %>% mutate(yrrep = paste(Year, Rep))
  )

fit_ccat_N_null_unrandom <- lm(
  log(pctN) ~ 1, 
  data = ccat_n %>% mutate(yrrep = paste(Year, Rep))
)



fit_ccat_bnf_null <- lmer(
  (pctfromBNF) ~ 1 + (1|yrrep), 
  data = ccat_n %>%   mutate(yrrep = paste(Year, Rep))
  )

fit_ccat_bnf_null_unrandom <- lm(
  (pctfromBNF) ~ 1, 
  data = ccat_n %>%   mutate(yrrep = paste(Year, Rep))
)


fit_ccat_survival_null <-
  lmer(survival ~ 1 + (1|yrrep), 
      data = ccat_survival %>% mutate(yrrep = paste(Year, Rep))
)

fit_ccat_survival_null_unrandom <-
  lm(survival ~ 1, 
      data = ccat_survival %>% mutate(yrrep = paste(Year, Rep))
  )


fit_ccat_emergence_null <-
  lmer(emergence ~ 1 + (1|yrrep), 
      data = ccat_survival %>% mutate(yrrep = paste(Year, Rep))
)

fit_ccat_emergence_null_unrandom <-
  lm(emergence ~ 1, 
      data = ccat_survival %>% mutate(yrrep = paste(Year, Rep))
  )






anova_with_var_decomp <- function(model, nullmodel, unrandom, name) {
  r_m <- residuals(model)
  r_n <- residuals(nullmodel)
  r_u <- residuals(unrandom)
  
  pct_fixed = 1 - (var(r_m) / var(r_n))
  pct_random = 1 - (var(r_n) / var(r_u))
  
  rmse <- sqrt(mean(r_m^2, na.rm = T))
  mu <- mean(fitted(model) + r_m, na.rm = T)
  
  anova(model) %>% 
    as_tibble(rownames = "term") %>% 
    mutate(
      statistic = sprintf("%1.2f", `F value`),
      statistic = replace(statistic, `F value` < 0.01, "<0.01"),
      statistic = paste0(statistic, pval_to_star(`Pr(>F)`))
    ) %>% 
    select(term, statistic) %>% 
    bind_rows(
      tibble(
        term = c(
          "variance_expl_fixed", "variance_expl_random",
          "rmse", "cv"
          ),
        statistic = sprintf(
          "%1.3f", 
          c(pct_fixed, pct_random, rmse, rmse/mu)
          )
      )
    ) %>% 
    magrittr::set_colnames(c("term", name))
}


random_chsqs <- function(nullmodel, unmixedmodel) {
  anova(nullmodel, unmixedmodel) %>% 
    as_tibble() %>% 
    select(ChiSq = Chisq, p = `Pr(>Chisq)`) %>% 
    slice(2)
}



model_summaries <- 
  tibble(
    name = c(
      "emergence_pct",
      "survival_pct",       
      "bio_g_plant",        
      "gdd_to_25pct_flower",
      "pctN",               
      "pctBNF"
    ),
    model = list(
      fit_ccat_emergence,
      fit_ccat_survival,
      fit_ccat_bm,      
      fit_ccat_gdd,     
      fit_ccat_N,       
      fit_ccat_bnf     
    ),
    null = list(
      fit_ccat_emergence_null,
      fit_ccat_survival_null,
      fit_ccat_bm_null,      
      fit_ccat_gdd_null,     
      fit_ccat_N_null,       
      fit_ccat_bnf_null  
    ),
    unrandom = list(
      fit_ccat_emergence_null_unrandom,
      fit_ccat_survival_null_unrandom,
      fit_ccat_bm_null_unrandom,      
      fit_ccat_gdd_null_unrandom,     
      fit_ccat_N_null_unrandom,       
      fit_ccat_bnf_null_unrandom 
    )
  ) %>% 
  mutate(
    d = pmap(
      list(model, null, unrandom, name), 
      anova_with_var_decomp
    )
  ) %>% 
  pull(d) %>% 
  bind_cols() %>% 
  select(-matches("term.$"))

write.csv(model_summaries, "CCAT model summaries.csv", row.names = F)

repeffect_chisqs <- list(
  emergence = list(fit_ccat_emergence_null, 
                   fit_ccat_emergence_null_unrandom),
  survival = list(fit_ccat_survival_null,
                  fit_ccat_survival_null_unrandom),
  biomass = list(fit_ccat_bm_null, 
                 fit_ccat_bm_null_unrandom),
  gdd = list(fit_ccat_gdd_null, 
             fit_ccat_gdd_null_unrandom),
  N = list(fit_ccat_N_null,
           fit_ccat_N_null_unrandom),
  bnf = list(fit_ccat_bnf_null,
             fit_ccat_bnf_null_unrandom)
) %>% 
  map_dfr(~random_chsqs(.[[1]], .[[2]]), .id = "metric")

write.csv(repeffect_chisqs, "CCAT model chisq tests for random effects.csv", row.names = F)
