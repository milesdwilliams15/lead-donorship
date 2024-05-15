##############################################################
# Replication code for:                                      #
# "Elusive Collaboration? The Determinants of Lead Donorship #
#  in International Development"                             #
##############################################################


# SUMMARY -----------------------------------------------------------------

# This .R script file contains the code necessary to replicate
# all figures, tables, and analyses in the manuscript. Some of
# the analysis is computationally intensive, so it may take a
# while to run. 
#
# If you notice any errors, contact me at: 
# mdwilliams@denison.edu
#
# R version used:
# R version 4.2.1 (2022-06-23 ucrt) -- "Funny-Looking Kid"
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
#
# Clean environment prior to running:

rm(list = ls())


# LOAD PACKAGES -----------------------------------------------------------

# First, you need to check that the necessary packages are 
# available in your local installation of R.
#
# Packages and versions used on my machine:
#
#         Package    Version Built
# 1      coolorrr 0.0.0.9000 4.2.1
# 2      estimatr      1.0.0 4.2.1
# 3  geomtextpath      0.1.1 4.2.1
# 4        GGally      2.1.2 4.2.3
# 5    kableExtra      1.3.4 4.2.1
# 6          metR     0.14.0 4.2.3
# 7          mgcv      1.9-0 4.2.3
# 8  RColorBrewer      1.1-3 4.2.0
# 9        texreg     1.38.6 4.2.1
# 10    tidyverse      2.0.0 4.2.3
# 
# All packages except {coolorrr} are available on the CRAN and
# can be installed using the standard library("packageName")
# convention. To install {coolorrr} you should use:
#
# install.packages("devtools")
# devtools::install_github("milesdwilliams15/coolorrr")
#
# As of this writing, the code will work with all new 
# installations of the above packages. This may change in the
# future as new versions of the packages are deployed.
#
# If you are missing any packages, the following code will
# catch this, offer an error message, and report out the 
# packages you are missing.
#

c(
  "coolorrr",
  "estimatr",
  "geomtextpath",
  "GGally",
  "kableExtra",
  "metR",
  "mgcv",
  "RColorBrewer",
  "texreg",
  "tidyverse",
  "furrr",
  "ranger",
  "neuralnet"
) -> needed_packages

installed.packages()[, 1] -> your_packages

all(
  needed_packages %in%
    your_packages
) -> all_here

if(!all_here) {
  stop(
    paste0(
      "You are missing the following packages:\n\n",
      paste0(
        needed_packages[
          !(needed_packages %in% your_packages)
        ],
        collapse = ", "
      ),
      "\nPlease install them first, and then\n",
      "rerun the script."
    )
  )
}

# If all the packages are available then they can now be 
# loaded into the environment.

needed_packages |>
  lapply(\(x) do.call("library", list(x))) |>
  invisible()

# Other settings
set_theme()
set_palette()


# LOAD DATA ---------------------------------------------------------------

dt <- read_rds(
  here::here("final_data.rds")
)


# CREATE COMPOSITE MEASURES -----------------------------------------------

# Make the helper function to compute the SSC measure of donor-recipient
# ties and recipient need.

# A simple function to convert values to z-scores
stand <- function(x) (x - mean(x)) / sd(x)

# The objective function
hat_Z <- function(X, w) {
  # standardize values
  X <- apply(X, 2, function(x) stand(x))
  
  # the linear combination
  Z <- X %*% w
  
  # the the squared covariances
  covs.sqrd <- apply(X, 2, function(x) cov(x, Z)^2)
  
  # return the negative of the sum of the squared covariances
  opt <- -sum(covs.sqrd)
  return(opt)
}

# Function that executes numerical optimizer for
# the SSC objective.
find_Z <- function(X) {
  out <- optim(
    fn = hat_Z,
    par = rep(0, len = ncol(X)),
    X = X
  )
  Z <- X %*% out$par
  Z <- stand(Z)
  return(Z)
}

# Create the SSC derived measures of ties and need and add to the
# analysis dataset.
dt |>
  mutate(
    ties2 = find_Z(
      cbind(log(dist), ihs_trade, colony, ally)
    ),
    need2 = find_Z(
      cbind(log(income), log(pop), ihs_disaster, civilwar, fh_total)
    )
  ) -> dt

# Add PCA and MFA derived versions of these measures.
prcomp(
  with(dt, cbind(log(dist), ihs_trade, colony, ally)),
  scale = T
)$x[, 1] -> pca_ties
prcomp(
  with(dt, cbind(log(income), log(pop), ihs_disaster, fh_total, civilwar)),
  scale = T
)$x[, 1] -> pca_need
factanal(
  with(dt, cbind(log(dist), ihs_trade, colony, ally)),
  factors = 1,
  scores = "regression"
)$scores[, 1] -> mfa_ties
factanal(
  with(dt, cbind(log(income), log(pop), ihs_disaster, fh_total, civilwar)),
  factors = 1,
  scores = "regression"
)$scores[, 1] -> mfa_need

dt |>
  mutate(
    pca_ties = pca_ties, 
    pca_need = -pca_need, # PCA measure reverses direction 
    mfa_ties = mfa_ties, 
    mfa_need = mfa_need  
  ) -> dt

# The dataset should now be populated with the relevant composite measures.


# FIGURE 1 ----------------------------------------------------------------

expand_grid(
  p1 = seq(0, 1, len = 10),
  p2 = seq(0, 1, len = 10)
) |>
  mutate(
    y = (1 - p1) * p2
  ) -> plot_data
Derivate(
  y ~ p1 + p2, 
  data = plot_data
) -> dy
plot_data |>
  mutate(
    dy.p1 = dy$y.dp1,
    dy.p2 = dy$y.dp2
  ) |>
  ggplot() +
  aes(x = p1, y = p2, z = y) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white"
  ) +
  geom_textabline(
    slope = -1, intercept = 1,
    label = "Increasing Likelihood of Lead Donorship",
    color = "white",
    text_only = T,
    fontface = "bold",
    vjust = 1
  ) +
  ggpal("diverging", "fill", 0.5) +
  scale_x_continuous(
    breaks = 0:1,
    labels = c("High Opportunity\nLow Need",
               "Low Opportunity\nHigh Need")
  ) +
  scale_y_continuous(
    breaks = 0:1,
    labels = c("Even Returns",
               "Concentrated Returns")
  ) +
  theme(
    axis.text.x = element_text(
      hjust = 0:1
    ),
    axis.text.y = element_text(
      angle = 90, hjust = 0:1
    ),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  )


# FIGURE 2 ----------------------------------------------------------------

expand_grid(
  p1 = seq(0, 1, len = 10),
  p2 = seq(0, 1, len = 10)
) |>
  mutate(
    y = p1 * p2
  ) -> plot_data
Derivate(
  y ~ p1 + p2, 
  data = plot_data
) -> dy
plot_data |>
  mutate(
    dy.p1 = dy$y.dp1,
    dy.p2 = dy$y.dp2
  ) |>
  ggplot() +
  aes(x = p1, y = p2, z = y) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white"
  ) +
  geom_textabline(
    slope = 1, intercept = 0,
    label = "Increasing Likelihood of Lead Donorship",
    color = "white",
    text_only = T,
    fontface = "bold",
    vjust = 1
  ) +
  ggpal("diverging", "fill", 0.5) +
  scale_x_continuous(
    breaks = 0:1,
    labels = c("High Opportunity\nLow Need",
               "Low Opportunity\nHigh Need")
  ) +
  scale_y_continuous(
    breaks = 0:1,
    labels = c("Even Returns",
               "Concentrated Returns")
  ) +
  theme(
    axis.text.x = element_text(
      hjust = 0:1
    ),
    axis.text.y = element_text(
      angle = 90, hjust = 0:1
    ),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  )


# MEASURE LEAD DONORSHIP --------------------------------------------------

# A rolling average function
roll_ave <- function(x) {
  x <- xn <- c(NA, NA, x)
  for(i in 3:length(x)) {
    x[i] <- mean(xn[(i-2):i], na.rm=T)
  }
  x <- x[3:length(x)]
  return(x)
}

# A modified function to help measure lead donorship
# per the longitudinal criteria in Steinwand 2015
stein_ave <- function(x) {
  x <- xn <- c(NA, NA, NA, NA, x, NA, NA, NA, NA)
  for(i in 5:(length(x)-4)) {
    x[i] <- (mean(xn[(i - 4):(i + 4)], na.rm=T)>=(5/9))+0
    x2 <- 1-any(cumsum(diff(cumsum(1-xn[(i - 4):(i + 4)])))>=2)
    x[i] <- (x[i]*x2)+0
  }
  x <- x[5:(length(x)-4)]
  x
}

# Add measure of lead donorship to the data
#
# Criteria:
# 1.  The Herfindahl Index (HI) for a given recipient in a given year is
#     greater than the median HI for the data sample;
# 2.  The difference between the total amount of aid given by the top
#     donor relative to the next largest donor to a recipient in a year is
#     greater than the sample median of this difference;
# 3.  The share of aid to a recipient in a given year from the top donor
#     is greater than the sample median of aid shares.
# 4.  For a given year, I calculate the proportion of times a donor meets
#     the criteria specified in points 1-3 in that year and the two
#     previous.
dt |>
  group_by(recipient, year) |>
  mutate(
    aid_share = aid / sum(aid),
    HI = sum(aid_share^2),
    top_donor = (aid==max(aid)),
    top_share = max(aid_share),
    aid_1st = max(aid),
    aid_2nd = sort(aid, decreasing = T)[2],
    aid_diff = aid_1st - aid_2nd
  ) |>
  ungroup() |>
  
  ## within year lead donorship
  mutate(
    lead_donor_year = ifelse(
      (HI > median(HI)) &
        (aid_diff > median(aid_diff)) &
        (aid_share > median(aid_share)) &
        top_donor,
      1, 0
    )
  ) |>
  # make sliding scale of lead donorship
  arrange(year) |>
  group_by(donor, recipient) |>
  mutate(
    # measure proposed in paper with more flexible
    # longitudinal criteria
    lead_donor = roll_ave(lead_donor_year) %>%
      ifelse(. == 0.5, 1/3, .), ## fix issue where new lead donors get 
                                ## score of 1/2 in second year
    # version with Steinwand's (2015) longitudinal
    # criteria
    lead_donor_stein = stein_ave(lead_donor_year)
  ) |>
  ungroup() -> dt


# FIGURE 3 ----------------------------------------------------------------

dt |>
  mutate(
    lead_donor_stein = ifelse(
      year %in% 2000:2010,
      replace_na(lead_donor_stein, 0),
      NA_integer_
    )
  ) -> dt
dt |>
  group_by(recipient, year) |>
  summarize(
    lead_donor = max(lead_donor),
    HI = mean(HI),
    lead_donor_stein = max(lead_donor_stein)
  ) |>
  group_by(year) |>
  summarize(
    lead_donor = mean(lead_donor),
    HI = mean(HI),
    lead_donor_stein = mean(lead_donor_stein)
  ) -> year_dt
year_dt |>
  rename(
    "Lead Donorship" = lead_donor,
    "Lead Donorship (Steinwand)" = lead_donor_stein
  ) |>
  ungroup() |>
  pivot_longer(
    cols = c(HI, `Lead Donorship`, `Lead Donorship (Steinwand)`)
  ) %>%
  ggplot() +
  aes(
    x = year,
    y = value,
    linetype = name
  ) +
  geom_line(size = 0.75) +
  scale_x_continuous(
    n.breaks = 10
  ) +
  labs(
    x = "Year",
    y = "Yearly Average per Recipient",
    linetype = NULL
  ) 


# FIGURE 4 ----------------------------------------------------------------

dt |>
  select(ties2, need2) %>%
  rename(
    Ties = ties2,
    Need = need2
  ) |>
  pivot_longer(
    everything()
  ) |>
  ggplot() +
  aes(x = value) +
  geom_histogram(
    color = "black",
    fill = "gray"
  ) +
  facet_wrap(~ name, scales = "free") +
  labs(
    x = NULL,
    y = "Count"
  ) +
  scale_y_continuous(
    labels = scales::comma
  )


# FIGURE 5 ----------------------------------------------------------------

ties_dt <- dt %>%
  select(ties2, ihs_trade, dist, ally, colony) |>
  mutate(dist = log(dist))
names(ties_dt) <- c(
  "Ties", "Trade", "Distance",
  "Allies", "Colony"
)
ggcorr(
  ties_dt,
  label = T
) +
  ggpal("diverging", "fill") +
  theme(legend.position = "none")


# FIGURE 6 ----------------------------------------------------------------

need_dt <- dt |>
  select(need2, income, pop, ihs_disaster, fh_total, civilwar) %>%
  mutate(income = log(income),
         pop = log(pop))
names(need_dt) <- c(
  "Need", "Income", "Population", 
  "Disaster", "Democracy", "Conflict"
)
ggcorr(
  need_dt,
  label = T
) +
  ggpal("diverging", "fill") +
  theme(legend.position = "none")


# TABLE 1 W/ REGRESSION ANALYSIS ------------------------------------------

# The following five regression models appear in Table 1 in the 
# manuscript.

ols_fit1 <- lm_robust(
  ihs_aid ~ 
    ## dyadic factors
    ihs_trade + ihs_dist + colony + ally +
    ## recipient factors
    log(income) + log(pop) + ihs_disaster + civilwar + fh_total,
  data = dt,
  fixed_effects = ~ year + donor,
  se_type = "stata",
  clusters = paste0(recipient, donor)
)
ols_fit2 <- lm_robust(
  ihs_aid ~ 
    ## dyadic factors
    ties2 +
    ## recipient factors
    log(income) + log(pop) + ihs_disaster + civilwar + fh_total,
  data = dt,
  fixed_effects = ~ year + donor,
  se_type = "stata",
  clusters = paste0(recipient, donor)
)
ols_fit3 <- lm_robust(
  ihs_aid ~ 
    ## dyadic factors
    ihs_trade + ihs_dist + colony + ally +
    ## recipient factors
    need2,
  data = dt,
  fixed_effects = ~ year + donor,
  se_type = "stata",
  clusters = paste0(recipient, donor)
)
ols_fit4 <- lm_robust(
  ihs_aid ~ 
    ## dyadic factors
    ties2 +
    ## recipient factors
    need2,
  data = dt,
  fixed_effects = ~ year + donor,
  se_type = "stata",
  clusters = paste0(recipient, donor)
)
ols_fit5 <- lm_robust(
  ihs_aid ~ 
    ## dyadic factors
    ties2 +
    ## recipient factors
    need2 +
    ## interaction
    ties2:need2,
  data = dt,
  fixed_effects = ~ year + donor,
  se_type = "stata",
  clusters = paste0(recipient, donor)
)

# The regression table output created via texreg(). Provides
# regression table in Latex output.
texreg(
  list(ols_fit1, ols_fit2, ols_fit3, ols_fit4, ols_fit5),
  include.ci = F,
  caption = "OLS Estimates for Dyadic Aid Commitments",
  caption.above = T,
  custom.header = list("Bilateral ODA Commitments" = 1:5),
  custom.coef.names = c(
    "Trade (asinh)",
    "Distance (asinh)",
    "Colony",
    "Alliance",
    "Income (log)",
    "Population (log)",
    "Disaster (asinh)",
    "Civil War",
    "Democracy",
    "Ties",
    "Need",
    "Ties $\\times$ Need"
  ),
  custom.gof.rows = list(
    "Year FE" = rep("Yes", len = 5),
    "Donor FE" = rep("Yes", len = 5),
    "N" = rep(scales::comma(nrow(dt)), len = 5)
  ),
  include.nobs = F,
  include.nclusts = F,
)


# TABLE 2 W/ GAM ANALYSIS -------------------------------------------------

## Aggregate the data for the recipient-year analysis
dt |>
  group_by(
    year, recipient
  ) |>
  mutate(
    ## The SSC version of the measure
    ties_share = (ties2 - min(ties2)) / 
      sum(ties2 - min(ties2)),
    ## The PCA version of the measure
    pca_share = (pca_ties - min(pca_ties)) / 
      sum(pca_ties - min(pca_ties)),
    ## The MFA version of the measure
    mfa_share = (mfa_ties - min(mfa_ties)) / 
      sum(mfa_ties - min(mfa_ties))
  ) |>
  reframe(
    lead_donor = max(lead_donor, na.rm=T),
    lead_donor_stein = max(lead_donor_stein, na.rm=T),
    aid_hhi = log(sum(aid_share^2, na.rm=T)),
    ties_ssp = log(sum(ties_share^2, na.rm=T)),
    pca_ssp = log(sum(pca_share^2, na.rm = T)),
    mfa_ssp = log(sum(mfa_share^2, na.rm = T)),
    need = mean(need2, na.rm=T),
    pca_need = mean(pca_need, na.rm = T),
    mfa_need = mean(mfa_need, na.rm = T)
  ) |>
  mutate(
    num_dyad = as.numeric(as.factor(recipient))
  )-> agg_dt

## Estimate GAMs
gam_fit1 <- gam(
  lead_donor ~ s(ties_ssp) + s(need) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)
gam_fit2 <- gam(
  lead_donor ~ s(ties_ssp, need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)

## Make Table 2
texreg(
  list(Separate = gam_fit1, Joint = gam_fit2),
  custom.coef.map = list(
    "EDF: s(ties_ssp)" = "Ties Concentration",
    "EDF: s(need)" = "Recipient Need",
    "EDF: s(ties_ssp,need)" = "Joint Ties and Need"
  ),
  custom.header = list("Lead Donorship" = 1:2),
  caption = "Logistic GAM Estimates",
  caption.above = T,
  include.aic = F,
  include.bic = F,
  include.log = F,
  include.nsmooth = F,
  include.nobs = F,
  include.nclusts = F,
  digits = 3,
  custom.gof.rows = list(
    "Year FE" = rep("Yes", len = 2),
    "Recipient RE" = rep("Yes", len = 2),
    "N" = c(scales::comma(nrow(agg_dt)), 
            scales::comma(nrow(agg_dt))),
    "Recipients" = c(scales::comma(length(unique(agg_dt$recipient))),
                     scales::comma(length(unique(agg_dt$recipient))))
  )
)


# FIGURE 7 ----------------------------------------------------------------

## Produce a heat map of the model predictions

ties_range <- agg_dt$ties_ssp |> range()
need_range <- agg_dt$need |> range()
pred_data <- expand_grid(
  ties_ssp = seq(ties_range[1], ties_range[2], len = 20),
  need = seq(need_range[1], need_range[2], len = 20),
  num_dyad = 1,
  year = 2001
)
pred_data %>%
  mutate(
    pred = predict(gam_fit2, newdata = ., type = "response")
  ) -> pred_data
dy <- Derivate(pred ~ need + ties_ssp, data = pred_data)
pred_data |>
  mutate(
    dy.p1 = dy$pred.dneed,
    dy.p2 = dy$pred.dties_ssp
  ) -> pred_data
ggplot(pred_data) +
  aes(
    x = need,
    y = ties_ssp,
    z = pred
  ) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white",
    size = 0.1,
    skip = 1
  ) +
  ggpal(
    type = "diverging",
    aes = "fill",
    breaks = c(0.25, 0.5, 0.75),
    labels = scales::percent,
    midpoint = 0.5
  ) +
  labs(
    x = expression(
      "Recipient Need (" * nu[r][t] * ")"
    ),
    y = expression(
      "Concentration of Ties (" * tau[r][t] * ")" 
    ),
    fill = "Likelihood of\nLead Donorship",
    caption = expression(""%->%" Direction of Greater Lead Donorship")
  ) +
  theme(
    legend.position = "top"
  )


# TABLE 3 W/ GAM ANALYSIS -------------------------------------------------

## As a robustness check, fit models using measure of lead donorship
## based on Steinwand's more restrictive longitudinal criteria.
gam_fit3 <- gam(
  lead_donor_stein ~ s(ties_ssp) + s(need) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt |> 
    filter(between(lead_donor_stein, 0, 1)),
  method = "REML",
  family = quasibinomial
)
gam_fit4 <- gam(
  lead_donor_stein ~ s(ties_ssp, need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt |> filter(between(lead_donor_stein, 0, 1)),
  method = "REML",
  family = quasibinomial
)

## Make Table 4
nagg_dt <- agg_dt |>
  filter(between(lead_donor_stein, 0, 1))
texreg(
  list(Separate = gam_fit3, Joint = gam_fit4),
  custom.coef.map = list(
    "EDF: s(ties_ssp)" = "Ties Concentration",
    "EDF: s(need)" = "Recipient Need",
    "EDF: s(ties_ssp,need)" = "Joint Ties and Need"
  ),
  custom.header = list("Lead Donorship (Steinwand)" = 1:2),
  caption = "Logistic GAM Estimates",
  caption.above = T,
  include.aic = F,
  include.bic = F,
  include.log = F,
  include.nsmooth = F,
  include.nobs = F,
  include.nclusts = F,
  digits = 3,
  custom.gof.rows = list(
    "Year FE" = rep("Yes", len = 2),
    "Recipient RE" = rep("Yes", len = 2),
    "N" = c(scales::comma(nrow(nagg_dt)), 
            scales::comma(nrow(nagg_dt))),
    "Recipients" = c(scales::comma(length(unique(nagg_dt$recipient))),
                     scales::comma(length(unique(nagg_dt$recipient))))
  )
)


# FIGURE 8 ----------------------------------------------------------------

## Visualize model predictions with more restrictive measure

ties_range <- (agg_dt |> filter(
  between(lead_donor_stein, 0, 1)
))$ties_ssp |> range()
need_range <- (agg_dt |> filter(
  between(lead_donor_stein, 0, 1)
))$need |> range()
pred_data <- expand_grid(
  ties_ssp = seq(ties_range[1], ties_range[2], len = 20),
  need = seq(need_range[1], need_range[2], len = 20),
  num_dyad = 1,
  year = 2001
)
pred_data %>%
  mutate(
    pred = predict(gam_fit4, newdata = ., type = "response")
  ) -> pred_data
dy <- Derivate(pred ~ need + ties_ssp, data = pred_data)
pred_data %>%
  mutate(
    dy.p1 = dy$pred.dneed,
    dy.p2 = dy$pred.dties_ssp
  ) -> pred_data
ggplot(pred_data) +
  aes(
    x = need,
    y = ties_ssp,
    z = pred
  ) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white",
    size = 0.1,
    skip = 1
  ) +
  ggpal(
    type = "diverging",
    aes = "fill",
    breaks = c(0.25, 0.5, 0.75),
    labels = scales::percent,
    midpoint = 0.5
  ) +
  labs(
    x = expression(
      "Recipient Need (" * nu[r][t] * ")"
    ),
    y = expression(
      "Concentration of Ties (" * tau[r][t] * ")" 
    ),
    fill = "Likelihood of\nLead Donorship",
    caption = expression(""%->%" Direction of Greater Lead Donorship")
  ) +
  theme(
    legend.position = "top"
  )



# FIGURE 9 W/ ANALYSIS ----------------------------------------------------

## Replicate the analysis with PCA based measures
pca_fit1 <- gam(
  lead_donor ~ s(pca_ssp) + s(pca_need) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)
pca_fit2 <- gam(
  lead_donor ~ s(pca_ssp, pca_need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)

ties_range <- agg_dt$pca_ssp |> range()
need_range <- agg_dt$pca_need |> range()
pred_data <- expand_grid(
  pca_ssp = seq(ties_range[1], ties_range[2], len = 20),
  pca_need = seq(need_range[1], need_range[2], len = 20),
  num_dyad = 1,
  year = 2001
)
pred_data %>%
  mutate(
    pred = predict(pca_fit2, newdata = ., type = "response")
  ) -> pred_data
dy <- metR::Derivate(pred ~ pca_need + pca_ssp, data = pred_data)
pred_data %>%
  mutate(
    dy.p1 = dy$pred.dpca_need,
    dy.p2 = dy$pred.dpca_ssp
  ) -> pred_data
ggplot(pred_data) +
  aes(
    x = pca_need,
    y = pca_ssp,
    z = pred
  ) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white",
    size = 0.1,
    skip = 1
  ) +
  ggpal(
    type = "diverging",
    aes = "fill",
    breaks = c(0.25, 0.5, 0.75),
    labels = scales::percent,
    midpoint = 0.5
  ) +
  labs(
    x = expression(
      "Recipient Need (" * nu[r][t] * ")"
    ),
    y = expression(
      "Concentration of Ties (" * tau[r][t] * ")" 
    ),
    fill = "Likelihood of\nLead Donorship",
    caption = expression(""%->%" Direction of Greater Lead Donorship")
  ) +
  theme(
    legend.position = "top"
  ) 

# FIGURE 10 W/ ANALYSIS ---------------------------------------------------

# Replicate the analysis with MFA based measures
mfa_fit1 <- gam(
  lead_donor ~ s(mfa_ssp) + s(mfa_need) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)
mfa_fit2 <- gam(
  lead_donor ~ s(mfa_ssp, mfa_need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)

ties_range <- agg_dt$mfa_ssp |> range()
need_range <- agg_dt$mfa_need |> range()
pred_data <- expand_grid(
  mfa_ssp = seq(ties_range[1], ties_range[2], len = 20),
  mfa_need = seq(need_range[1], need_range[2], len = 20),
  num_dyad = 1,
  year = 2001
)
pred_data %>%
  mutate(
    pred = predict(mfa_fit2, newdata = ., type = "response")
  ) -> pred_data
dy <- Derivate(pred ~ mfa_need + mfa_ssp, data = pred_data)
pred_data %>%
  mutate(
    dy.p1 = dy$pred.dmfa_need,
    dy.p2 = dy$pred.dmfa_ssp
  ) -> pred_data
ggplot(pred_data) +
  aes(
    x = mfa_need,
    y = mfa_ssp,
    z = pred
  ) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white",
    size = 0.1,
    skip = 1
  ) +
  ggpal(
    type = "diverging",
    aes = "fill",
    breaks = c(0.25, 0.5, 0.75),
    labels = scales::percent,
    midpoint = 0.5
  ) +
  labs(
    x = expression(
      "Recipient Need (" * nu[r][t] * ")"
    ),
    y = expression(
      "Concentration of Ties (" * tau[r][t] * ")" 
    ),
    fill = "Likelihood of\nLead Donorship",
    caption = expression(""%->%" Direction of Greater Lead Donorship")
  ) +
  theme(
    legend.position = "top"
  )


# FIGURE 11 W/ ANALYSIS ---------------------------------------------------

# Fit with lagged version of ties
agg_dt |>
  group_by(recipient) |>
  mutate(
    lag_ties_ssp = lag(ties_ssp, order_by = year, n = 2L)
  ) |>
  ungroup() -> agg_dt
gam(
  lead_donor ~ s(lag_ties_ssp, need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
) -> gam_fit5

ties_range <- agg_dt$lag_ties_ssp |> range(na.rm = T)
need_range <- agg_dt$need |> range()
pred_data <- expand_grid(
  lag_ties_ssp = seq(ties_range[1], ties_range[2], len = 20),
  need = seq(need_range[1], need_range[2], len = 20),
  num_dyad = 1,
  year = 2001
)
pred_data %>%
  mutate(
    pred = predict(gam_fit5, newdata = ., type = "response")
  ) -> pred_data
dy <- Derivate(pred ~ need + lag_ties_ssp, data = pred_data)
pred_data %>%
  mutate(
    dy.p1 = dy$pred.dneed,
    dy.p2 = dy$pred.dlag_ties_ssp
  ) -> pred_data
ggplot(pred_data) +
  aes(
    x = need,
    y = lag_ties_ssp,
    z = pred
  ) +
  geom_contour_fill() +
  geom_arrow(
    aes(dx = dy.p1, dy = dy.p2),
    color = "white",
    size = 0.1,
    skip = 1
  ) +
  ggpal(
    type = "diverging",
    aes = "fill",
    breaks = c(0.25, 0.5, 0.75),
    labels = scales::percent,
    midpoint = 0.5
  ) +
  labs(
    x = expression(
      "Recipient Need (" * nu[r][t] * ")"
    ),
    y = expression(
      "Concentration of Ties (" * tau[r][t] * ")" 
    ),
    fill = "Likelihood of\nLead Donorship",
    caption = expression(""%->%" Direction of Greater Lead Donorship")
  ) +
  theme(
    legend.position = "top"
  )


# FIGURE 12 ---------------------------------------------------------------

agg_dt |>
  mutate(
    "Ties Concentration" = cut(
      ties_ssp,
      breaks = quantile(ties_ssp),
      labels = paste0(
        "Ties: ", c("1st", "2nd", "3rd", "4th")
      ),
      include.lowest = T,
      ordered_result = T
    ) |>
      factor(levels = paste0(
        "Ties: ", c("4th", "3rd", "2nd", "1st")
      )),
    "Need" = cut(
      need,
      breaks = quantile(need),
      labels = paste0(
        "Need: ", c("1st", "2nd", "3rd", "4th")
      ),
      include.lowest = T,
      ordered_result = T
    )
  ) |>
  group_by(
    `Ties Concentration`,
    `Need`
  ) |>
  count(lead_donor) |>
  ggplot() +
  aes(
    x = lead_donor,
    y = n
  ) +
  geom_col(
    fill = "gray",
    color = "black"
  ) +
  geom_point() +
  geom_line() +
  geom_text(
    aes(label = n),
    vjust = -0.3,
    fontface = "bold"
  ) +
  labs(
    x = "Strength of Lead Donorship",
    y = "Count"
  ) +
  scale_x_continuous(
    n.breaks = 4,
    labels = c("0", "1/3", "2/3", "1")
  ) +
  ylim(c(0, 165)) +
  facet_grid(
    `Ties Concentration` ~
      `Need`
  ) 

# FIGURE A1 (online appendix a) -------------------------------------------

dt |>
  mutate(
    donor_code = countrycode::countrycode(
      donor, "country.name", "iso3c"
    )
  ) |>
  filter(
    recipient %in% c("Mexico", "Mozambique")
  ) |>
  group_by(recipient, year) |>
  transmute(
    donor_code = donor_code,
    aid_share = aid / sum(aid),
    aid_1st = aid_share == sort(aid_share, decreasing = T)[1],
    aid_2nd = aid_share == sort(aid_share, decreasing = T)[2]
  ) |>
  filter(
    aid_1st | aid_2nd
  ) |>
  ungroup() |>
  group_split(aid_2nd) -> plt_dt

ggplot(plt_dt[[1]]) +
  aes(year, aid_share, label = donor_code) +
  geom_col(
    color = 'black',
    fill = 'firebrick'
  ) +
  geom_text(vjust = .5, hjust = 1, angle = 90,
            color = "white", fontface = "bold") +
  geom_text(
    aes(label = round(aid_share*100)),
    vjust = -.5
  ) +
  geom_col(
    data = plt_dt[[2]],
    aes(year, aid_share),
    fill = "lightgrey",
    alpha = 0.4,
    width = 0.75
  ) +
  facet_wrap(~ recipient) + 
  scale_x_continuous(
    breaks = 1996:2014
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Aid Share"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 


# FIGURE A2 (online appendix a) -------------------------------------------

dt97 <- filter(dt, year > 1997)
dt97 |>
  mutate(
    donor_code = countrycode::countrycode(
      donor, "country.name", "iso3c"
    )
  ) |>
  filter(
    recipient %in% c("Mexico", "Mozambique")
  ) |>
  group_by(recipient, year) |>
  transmute(
    donor_code = donor_code,
    lead_donor = lead_donor,
    aid_share = aid / sum(aid),
    aid_1st = aid_share == sort(aid_share, decreasing = T)[1]
  ) |>
  filter(
    aid_1st
  ) |>
  ungroup() |>
  mutate(
    lead_donor = case_when(
      lead_donor == 0 ~ "0",
      lead_donor == 1/3 ~ "1/3",
      lead_donor == 2/3 ~ "2/3",
      TRUE ~ "1"
    ) |>
      factor(
        levels = c("0", "1/3", "2/3", "1")
      )
  ) -> plt_dt
ggplot(plt_dt) +
  aes(year, aid_share, label = donor_code,
      fill = lead_donor) +
  geom_col() +
  geom_text(vjust = .5, hjust = 1, angle = 90,
            color = "white", fontface = "bold") +
  geom_text(
    aes(label = round(aid_share*100)),
    vjust = -.5
  ) +
  facet_wrap(~ recipient) + 
  scale_x_continuous(
    breaks = 1998:2014
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_fill_brewer(
    palette = "Reds"
  ) +
  labs(
    x = NULL,
    y = "Aid Share",
    fill = "Lead Donorship"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = c(0.8, 0.75),
    legend.background = element_rect(color = "black")
  ) 


# TABLE B1 W/ ANALYSIS (online appendix b) --------------------------------

# Code to test the predictive performance of SSC, PCA, and MFA measures


full_data <- dt %>%
  select(
    donor, recipient, year, ihs_aid, need2, ties2, pca_need, pca_ties,
    mfa_need, mfa_ties
  ) %>%
  mutate(
    dyad = paste0(donor, "-", recipient)
  ) %>%
  group_by(donor, year) %>%
  mutate(
    across(ihs_aid:mfa_ties, ~ .x - mean(.x))
  ) %>%
  ungroup %>%
  mutate(
    year = as.factor(year)
  )

# determine which dyads will be in the training vs test sets
dyads <- unique(full_data$dyad)
dyad_num <- 1:length(dyads)
set.seed(111)
dyad_train <- sample(dyad_num, round(0.7 * length(dyad_num)))

# split the data
train <- full_data %>% 
  filter(
    dyad %in% dyads[dyad_train]
  )
test  <- full_data %>%
  filter(
    dyad %in% dyads[-dyad_train]
  )

## TRAIN THE MODELS AND VALIDATE WITH OOS PREDICTIONS
forms <- list(
  ihs_aid ~ need2 + ties2,
  ihs_aid ~ pca_need + pca_ties,
  ihs_aid ~ mfa_need + mfa_ties
)

## fit, predict, and compare using test data
# linear
lm_fits <- map(
  forms,
  ~ lm(., data = train)
)

# generalized additive model
gm_fits <- map(
  list(
    ihs_aid ~ s(need2) + s(ties2),
    ihs_aid ~ s(pca_need) + s(pca_ties),
    ihs_aid ~ s(mfa_need) + s(mfa_ties)
  ),
  ~ gam(., data = train)
)

# random forests
rf_fits <- map(
  forms,
  ~ ranger(., data = train)
)

# neural nets
nn_fits <- map(
  forms,
  ~ neuralnet(., data = train, 
              hidden = 2, 
              linear.output = T, 
              threshold = 2)
)

# predictions
lm_preds <- map(
  lm_fits,
  ~ predict(., newdata = test)
)
gm_preds <- map(
  gm_fits,
  ~ predict(., newdata = test)
)
rf_preds <- map(
  rf_fits,
  ~ predict(., data = test)$predictions
)
nn_preds <- map(
  nn_fits,
  ~ predict(., test)[, 1]
)

# percent variance explained
mse <- function(x, y) mean((x - y)^2)
lm_gof <- map_dfc(
  lm_preds,
  ~ mse(., test$ihs_aid)
)
gm_gof <- map_dfc(
  gm_preds,
  ~ mse(., test$ihs_aid)
)
rf_gof <- map_dfc(
  rf_preds,
  ~ mse(., test$ihs_aid)
)
nn_gof <- map_dfc(
  nn_preds,
  ~ mse(., test$ihs_aid)
)
diff_test <- function(
    preds
) {
  x <- (preds[[1]] - test$ihs_aid)^2
  y <- (preds[[2]] - test$ihs_aid)^2
  z <- (preds[[3]] - test$ihs_aid)^2
  t.test(
    y - x
  ) -> y_test
  t.test(
    z - x
  ) -> z_test
  return(
    data.frame(y_stat = y_test$statistic,
               z_stat = z_test$statistic,
               y_diff = y_test$p.value,
               z_diff = z_test$p.value)
  )
}
lm_diff <- diff_test(lm_preds)
gm_diff <- diff_test(gm_preds)
rf_diff <- diff_test(rf_preds)
nn_diff <- diff_test(nn_preds)
col_names <- c("SSC", "PCA", "MFA")
out <- bind_rows(
  lm_gof,
  gm_gof,
  rf_gof,
  nn_gof
)
names(out) <- col_names
#names(out_se) <- col_names
out <- mutate(out, Model = c("Linear", "GAM", "Random Forests", "Neural Net")) %>%
  select(Model, everything())
pvals <- bind_rows(
  lm_diff,
  gm_diff,
  rf_diff,
  nn_diff
)
out$PCA <- paste0(
  round(out$PCA, 3), 
  " (", round(pvals$y_stat, 2), ")",
  gtools::stars.pval(pvals$y_diff)
)
out$MFA <- paste0(
  round(out$MFA, 3),
  " (", round(pvals$z_stat, 2), ")",
  gtools::stars.pval(pvals$z_diff)
)

## predict and compare using training data (within sample predictions)
# predictions
lm_preds <- map(
  lm_fits,
  ~ predict(., newdata = train)
)
gm_preds <- map(
  gm_fits,
  ~ predict(., newdata = train)
)
rf_preds <- map(
  rf_fits,
  ~ predict(., data = train)$predictions
)
nn_preds <- map(
  nn_fits,
  ~ predict(., train)[, 1]
)

# percent variance explained
mse <- function(x, y) mean((x - y)^2)
se_mse <- function(x, y) sd((x - y)^2) / sqrt(length(x))
lm_gof <- map_dfc(
  lm_preds,
  ~ mse(., train$ihs_aid)
)
gm_gof <- map_dfc(
  gm_preds,
  ~ mse(., train$ihs_aid)
)
rf_gof <- map_dfc(
  rf_preds,
  ~ mse(., train$ihs_aid)
)
nn_gof <- map_dfc(
  nn_preds,
  ~ mse(., train$ihs_aid)
)
diff_test <- function(
    preds
) {
  x <- (preds[[1]] - train$ihs_aid)^2
  y <- (preds[[2]] - train$ihs_aid)^2
  z <- (preds[[3]] - train$ihs_aid)^2
  t.test(
    y - x
  ) -> y_test
  t.test(
    z - x
  ) -> z_test
  return(
    data.frame(y_stat = y_test$statistic,
               z_stat = z_test$statistic,
               y_diff = y_test$p.value,
               z_diff = z_test$p.value)
  )
}
lm_diff <- diff_test(lm_preds)
gm_diff <- diff_test(gm_preds)
rf_diff <- diff_test(rf_preds)
nn_diff <- diff_test(nn_preds)
col_names <- c("SSC", "PCA", "MFA")
out_train <- bind_rows(
  lm_gof,
  gm_gof,
  rf_gof,
  nn_gof
)
names(out_train) <- col_names
out_train <- mutate(out_train, Model = c("Linear", "GAM", "Random Forests", "Neural Net")) %>%
  select(Model, everything())
pvals <- bind_rows(
  lm_diff,
  gm_diff,
  rf_diff,
  nn_diff
)
out_train$PCA <- paste0(
  round(out_train$PCA, 3),
  " (", round(pvals$y_stat, 2), ")",
  gtools::stars.pval(pvals$y_diff)
)
out_train$MFA <- paste0(
  round(out_train$MFA, 3),
  " (", round(pvals$z_stat, 2), ")",
  gtools::stars.pval(pvals$z_diff)
)

## show the results
new_out <- left_join(out, out_train, by = 'Model')
colnames(new_out) <- c(
  'Model',
  'SSC',
  'PCA',
  'MFA',
  'SSC',
  'PCA',
  'MFA'
)

# Table B1
kable(
  "latex",
  new_out,
  booktabs = T,
  digits = 3,
  caption = "Prediction MSE",
  align = c('l', rep('c', len = 6))
) %>%
  add_footnote(
    label = c("* p < 0.05; ** p < 0.01; *** p < 0.001.",
              "Sig. tests are from pair-wise t-tests (t-stat in parentheses) of the squared error in",
              "predictions using SSC derived measures relative to one of the alternatives. Entries not", 
              "in parentheses are the prediction mean squared error."),
    notation = "none"
  ) %>%
  add_header_above(
    c(" ", "Test Predictions" = 3, "Train Predictions" = 3)
  ) %>%
  column_spec(
    2:7,
    width = ".7in"
  )



# TABLE C1 (online appendix c) --------------------------------------------

mf <- cbind(
  model.frame(ols_fit1),
  model.frame(ols_fit5)[, -1]
)
mf |>
  pivot_longer(
    cols = everything()
  ) |>
  filter(
    !is.infinite(value)
  ) |>
  group_by(name) |>
  summarize(
    Mean = mean(value, na.rm = T),
    Median = median(value, na.rm = T),
    SD = sd(value, na.rm = T),
    Min. = min(value, na.rm = T),
    Max. = max(value, na.rm = T),
    N = sum(!is.infinite(value)) |>
      scales::comma()
  ) |>
  mutate(
    name = c(
      "Alliance",
      "Civil War", 
      "Colony",
      "Democracy",
      "ODA (asinh)",
      "Distaster (asinh)",
      "Distance (asinh)",
      "Trade (asinh)",
      "Income (log)",
      "Population (log)",
      "Recipient Need",
      "Donor-Recipient Ties"
    )
  ) |>
  rename(
    " " = name
  ) |>
  kable(
    "latex",
    caption = "Summary Statistics for Dyadic Aid Analysis",
    digits = 3,
    booktabs = T,
    linesep = ""
  ) 


# TABLE C2 (online appendix c) --------------------------------------------

agg_dt |>
  select(
    lead_donor,
    lead_donor_stein,
    ties_ssp,
    need
  ) |>
  pivot_longer(
    cols = everything()
  ) |>
  filter(
    !is.infinite(value)
  ) |>
  group_by(name) |>
  summarize(
    Mean = mean(value, na.rm = T),
    Median = median(value, na.rm = T),
    SD = sd(value, na.rm = T),
    Min. = min(value, na.rm = T),
    Max. = max(value, na.rm = T),
    N = sum(!is.infinite(value)) |>
      scales::comma()
  ) |>
  mutate(
    name = c("Lead Donorship",
             "Lead Donorship (Steinwand 2015 version)",
             "Recipient Need",
             "Donor-Recipient Ties Concentration")
  ) |>
  rename(
    " " = name
  ) |>
  kable(
    "latex",
    caption = "Summary Statistics for Lead Donorship Analysis",
    digits = 3,
    booktabs = T
  )
