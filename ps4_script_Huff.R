library(tidyverse)
library(fpp3)
library(BVAR)
library(future)

# ---- Settings ----
quick <- FALSE            # TRUE = fast test run (don't report these numbers)
plan(multisession)        # parallel fitting for fable; harmless if it doesn't help

# =========================================================
# Part 2: ARIMA vs. ETS
# =========================================================
# 2.1a
liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv", na = ".",
  col_types = cols(observation_date = col_date(format = "%Y-%m-%d"),
                   MRTSSM4453USN = col_double())
)
liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month <= yearmonth("2019 Dec"))
has_gaps(liquor)

# 2.1b
p1 <- liquor |> autoplot(sales) +
  labs(title = "U.S. liquor store sales", x = "Month", y = "Millions of dollars")
p2 <- liquor |> gg_season(log(sales)) +
  labs(title = "Seasonal plot: log liquor sales", x = "Month", y = "log(sales)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
p1; p2
ggsave("results/autoplot_sales.png", p1, width = 8, height = 4.5)
ggsave("results/season_plot.png",    p2, width = 8, height = 5)

# 2.2a
liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)),
        ets   = ETS(log(sales)))
liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

# 2.2b  (dof = 1 + 2 + 2 + 1 = 6 for ARIMA(1,0,2)(2,1,1)[12])
liquor_fit |> select(arima) |> gg_tsresiduals()
liquor_fit |> select(arima) |> augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

# 2.3a  (slow: 10+ minutes; saved after the first run)
liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = if (quick) 12 else 1) |>
  filter(.id < max(.id))
liquor_cv |> as_tibble() |> distinct(.id) |> nrow()   # 216 origins when quick = FALSE

cv_file <- if (quick) "results/liquor_cv_fit_quick.rds" else "results/liquor_cv_fit.rds"
if (file.exists(cv_file)) {
  liquor_cv_fit <- readRDS(cv_file)
} else {
  liquor_cv_fit <- liquor_cv |>
    model(snaive = SNAIVE(sales),
          ets    = ETS(log(sales)),
          arima  = ARIMA(log(sales)))
  saveRDS(liquor_cv_fit, cv_file)
}

# 2.3b
liquor_cv_fc <- liquor_cv_fit |> forecast(h = 1)
range(liquor_cv_fc$Month)
liquor_acc <- liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))
liquor_acc
saveRDS(liquor_acc, "results/liquor_acc.rds")

# 2.4
liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales), by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

dm_arima_ets    <- forecast::dm.test(liquor_errs$arima, liquor_errs$ets,    h = 1, power = 2)
dm_arima_snaive <- forecast::dm.test(liquor_errs$arima, liquor_errs$snaive, h = 1, power = 2)
dm_arima_ets
dm_arima_snaive
saveRDS(list(arima_ets = dm_arima_ets, arima_snaive = dm_arima_snaive),
        "results/dm_liquor.rds")

# =========================================================
# Part 3: Forecast combination
# =========================================================
liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)
combo_acc <- liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))
combo_acc
saveRDS(combo_acc, "results/combo_acc.rds")

cor(liquor_errs$arima, liquor_errs$ets)

# =========================================================
# Part 4: VAR vs. BVAR
# =========================================================
# 4.1a-b
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2)) |>
  mutate(across(-sasdate, as.numeric))     # makes sure series are numeric

growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),
            cons   = growth(PCECC96),
            inv    = growth(GPDIC1),
            emp    = growth(PAYEMS),
            infl   = growth(GDPCTPI),
            unrate = difference(UNRATE),
            ffr    = difference(FEDFUNDS),
            gs10   = difference(GS10)) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()
nrow(y)      # 240
head(y)

# 4.1c
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]
length(origins)     # 120 forecasts per variable (1990 Q1 to 2019 Q4)

# 4.2: OLS VARs
var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)
  M <- ncol(z)
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))
rmse_by_var(fc_var1)
rmse_by_var(fc_var4)

# 4.3a-b: BVAR (about 1-2 minutes per lambda)
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))
    fc <- rbind(fc, x %*% coef(fit_T))
  }
  fc
}

if (file.exists("results/fc_bvar_03.rds")) {
  fc_bvar_03 <- readRDS("results/fc_bvar_03.rds")
} else {
  set.seed(422)
  fc_bvar_03 <- bvar_cv(0.3)
  saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")
}

if (file.exists("results/fc_bvar_005.rds")) {
  fc_bvar_005 <- readRDS("results/fc_bvar_005.rds")
} else {
  set.seed(422)
  fc_bvar_005 <- bvar_cv(0.05)
  saveRDS(fc_bvar_005, "results/fc_bvar_005.rds")
}

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3`  = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)
saveRDS(var_results, "results/var_results.rds")
var_results |> knitr::kable(digits = 3)

# 4.3c: DM tests, GDP (column 1) and consumption (column 2)
e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]
dm_gdp <- forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

e_var4_c <- actual[, "cons"] - fc_var4[, 2]
e_bvar_c <- actual[, "cons"] - fc_bvar_03[, 2]
dm_cons <- forecast::dm.test(e_var4_c, e_bvar_c, h = 1, power = 2)

dm_gdp; dm_cons
saveRDS(list(gdp = dm_gdp, cons = dm_cons), "results/dm_var.rds")

# =========================================================
# Part 5: Factor model with FRED-MD
# =========================================================
# 5.1a
md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,
         "2" = d(x),
         "3" = d(d(x)),
         "4" = log(x),
         "5" = 100 * d(log(x)),
         "6" = 100 * d(d(log(x))),
         "7" = d(x / dplyr::lag(x) - 1))
}

# 5.1b
fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
fm <- fm[, colSums(is.na(fm)) == 0]
dim(fm)     # 718 x 122

# 5.1c
tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")]

# 5.2a
Z <- fm |> select(-date) |> as.matrix()
pc <- prcomp(scale(Z))
summary(pc)$importance[, 1:5]

# 5.2b
F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))
tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

# 5.2c
sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]

# 5.3  (several minutes: PCA is re-run at every origin)
ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))

e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}
ip_rmse <- c(RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)))
ip_rmse     # expect 0.617 and 0.589
saveRDS(ip_rmse, "results/ip_rmse.rds")

# 5.4  (also slow)
direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)
  g  <- lv$g
  
  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                                       varestimator = "bartlett")$p.value)
}

cpi_race <- direct_race("CPIAUCSL")
cpi_race
saveRDS(cpi_race, "results/cpi_race.rds")
# percent RMSE reduction: 100 * (cpi_race$RMSE_AR - cpi_race$RMSE_FA) / cpi_race$RMSE_AR

plan(sequential)