
# -------------------------------------------------------------------------
# Robustness checks using PCA and MFA constructed measures
# -------------------------------------------------------------------------

## PACKAGES ----

library(tidyverse)
library(furrr)
library(ranger)
library(neuralnet)
library(kableExtra)
library(coolorrr)
set_theme()
set_palette()
plan(multisession, workers = 3)

## DATA ----
dt <- read_csv(
  here::here("01_data", "final_data.csv")
)

## CONSTRUCT THE MEASURES ----

# make SSC ties and need measures

stand <- function(x) (x - mean(x)) / sd(x)
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

# implement the objective (returns ssc)
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

dt <- dt |>
  mutate(
    ties2 = find_Z(
      cbind(log(dist), ihs_trade, colony, ally)
    ),
    need2 = find_Z(
      cbind(log(income), log(pop), ihs_disaster, civilwar, fh_total)
    )
  )

# use PCA and MFA
pca_ties <- prcomp(
  with(dt, cbind(log(dist), ihs_trade, colony, ally)),
  scale = T
)$x[, 1]
pca_need <- prcomp(
  with(dt, cbind(log(income), log(pop), ihs_disaster, fh_total, civilwar)),
  scale = T
)$x[, 1]
fca_ties <- factanal(
  with(dt, cbind(log(dist), ihs_trade, colony, ally)),
  factors = 1,
  scores = "regression"
)$scores[, 1]
fca_need <- factanal(
  with(dt, cbind(log(income), log(pop), ihs_disaster, fh_total, civilwar)),
  factors = 1,
  scores = "regression"
)$scores[, 1]

dt <- dt |>
  mutate(
    pca_ties = pca_ties,
    pca_need = -pca_need,
    mfa_ties = fca_ties,
    mfa_need = fca_need
  )

## SPLIT THE DATA INTO TRAINING AND TEST SETS ----
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
new_out
save(
  new_out,
  file = here::here("03_report", "pred_check.R")
)


## REPLICATE ANALYSIS WITH ALTERNATIVE MEASURES ----

roll_ave <- function(x) {
  x <- xn <- c(NA, NA, x)
  for(i in 3:length(x)) {
    x[i] <- mean(xn[(i-2):i], na.rm=T)
  }
  x <- x[3:length(x)]
  return(x)
}
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

# define lead donorship bilaterally ----
dt %>%
  # mutate(
  #   ties = - ties,
  #   need = - need
  # ) %>%
  group_by(recipient, year) %>%
  mutate(
    aid_share = aid / sum(aid),
    HI = sum(aid_share^2),
    top_donor = (aid==max(aid)),
    top_share = max(aid_share),
    aid_1st = max(aid),
    aid_2nd = sort(aid, decreasing = T)[2],
    aid_diff = aid_1st - aid_2nd,
    ## aid_share = aid / sum(aid),
    ## aid_share_gt = aid_share > median(aid_share),
    ## top_donor = (aid==max(aid))+0,
    ## fst_aid = max(aid),
    ## snd_aid = sort(aid, decreasing = T)[2],
    ties_share = (ties - min(ties)) / sum(ties - min(ties))
  ) %>%
  ungroup %>%
  
  ## Steinwand's measure
  mutate(
    lead_donor_year = ifelse(
      (HI > median(HI)) &
        (aid_diff > median(aid_diff)) &
        (aid_share > median(aid_share)) &
        top_donor,
      1, 0
    )
  ) %>%
  arrange(year) %>%
  group_by(donor, recipient) %>%
  mutate(
    lead_donor = roll_ave(lead_donor_year),
    lead_donor_stein = stein_ave(lead_donor_year)
  ) %>%
  ungroup -> dt

# show distributions ----

vals_to_keep <- dt %>%
  select(ties2, need2, pca_ties, pca_need, mfa_ties, mfa_need) 
colnames(vals_to_keep) <- c(
  "1. Ties (SSC)", "2. Need (SSC)",
  "3. Ties (PCA)", "4. Need (PCA)",
  "5. Ties (MFA)", "6. Need (MFA)"
)
vals_to_keep |>
  pivot_longer(
    everything()
  ) %>%
  ggplot() +
  aes(x = value) +
  geom_histogram(
    color = "black",
    fill = "gray"
  ) +
  facet_wrap(~ name, scales = "free",
             ncol = 2) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  scale_y_continuous(
    labels = scales::comma
  )

# show correlations ----
ties_dt <- dt %>%
  select(ties2, pca_ties, mfa_ties, ihs_trade, dist, ally, colony) |>
  mutate(dist = log(dist))
names(ties_dt) <- c(
  "SSC", "PCA", "MFA", 
  "Trade", "Distance",
  "Allies", "Colony"
)
ggcorr(
  ties_dt,
  label = T
) +
  ggpal("diverging", "fill") +
  theme(legend.position = "none")

need_dt <- dt %>%
  select(need2, pca_need, mfa_need, income, pop, ihs_disaster, fh_total, civilwar) %>%
  mutate(income = log(income),
         pop = log(pop))
names(need_dt) <- c(
  "SSC", "PCA", "MFA", "Income", "Population", 
  "Disaster", "Democracy", "Conflict"
)
ggcorr(
  need_dt,
  label = T
) +
  ggpal("diverging", "fill") +
  theme(legend.position = "none")

# Aggregate the data ----
hhi <- function(x) {
  tot <- sum(x, na.rm=T)
  pct <- x / tot
  ssp <- sum(pct^2, na.rm=T)
  ssp # return
}
dt %>%
  group_by(
    year, recipient
  ) %>%
  mutate(
    ties_share = (ties2 - min(ties2)) / sum(ties2 - min(ties2)),
    pca_share = (pca_ties - min(pca_ties)) / 
      sum(pca_ties - min(pca_ties)),
    mfa_share = (mfa_ties - min(mfa_ties)) / 
      sum(mfa_ties - min(mfa_ties))
  ) %>%
  reframe(
    lead_donor = max(lead_donor, na.rm=T),
    lead_donor_stein = max(lead_donor_stein, na.rm=T),
    aid_hhi = log(sum(aid_share^2, na.rm=T)),
    total_aid = sum(sinh(ihs_aid)),
    ties_ssp = log(sum(ties_share^2, na.rm=T)),
    need = mean(need2, na.rm=T),
    pca_ssp = log(sum(pca_share^2, na.rm=T)),
    pca_need = mean(pca_need, na.rm = T),
    mfa_ssp = log(sum(mfa_share^2, na.rm = T)),
    mfa_need = mean(mfa_need, na.rm = T)
  ) %>%
  mutate(
    num_dyad = as.numeric(as.factor(recipient))
  )-> agg_dt

## Estimate GAM
library(mgcv)

## with PCA based measures
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
ggsave(
  "pca_results.png",
  height = 4.5,
  width = 4
)

## with MFA based measures
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
ggsave(
  "mfa_results.png",
  height = 4.5,
  width = 4
)

## using fragmentation as the outcome
frag_fit1 <- gam(
  aid_hhi ~ s(ties_ssp) + s(need) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML",
  family = quasibinomial
)
frag_fit2 <- gam(
  aid_hhi ~ s(ties_ssp, need, k = 15) + 
    as.factor(year) + 
    s(num_dyad, bs = "re"),
  data = agg_dt,
  method = "REML"
)

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
    pred = predict(frag_fit2, newdata = .)
  ) -> pred_data
dy <- metR::Derivate(pred ~ need + ties_ssp, data = pred_data)
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
    type = "sequential",
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
  ) +
  geom_point(
    data = agg_dt |> mutate(pred = 0)
  )
ggsave(
  "frag_results.png",
  height = 4.5,
  width = 4
)

agg_dt |>
  group_by(
    need > median(need),
    ties_ssp > median(ties_ssp)
  ) |>
  summarize(
    x = median(aid_hhi)
  )
