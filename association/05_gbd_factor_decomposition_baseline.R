library(tidyverse)
library(arrow)
library(data.table)

rm(list = ls())
gc()

find_project_root <- function() {
  path <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, ".here"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root (missing .here file)")
    }
    path <- parent
  }
}

DATA_ROOT <- Sys.getenv("AU_FIRE_ROOT", unset = file.path(find_project_root(), "data"))
setwd(file.path(DATA_ROOT, "burden"))

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

INPUT_DIR <-
  "gbd_3factor_decomposition_total_uncertainty_adjacent"

INPUT_FILE <- file.path(
  INPUT_DIR,
  "gbd_factor_inputs_by_sim_model_period.parquet"
)

OUT_DIR <-
  "gbd_3factor_decomposition_total_uncertainty_baseline"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

BASELINE_PERIOD <- "2000-2020"

FUTURE_PERIODS <- c(
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100"
)

LEVEL_ORDER <- c(
  "Urban",
  "P_Indigenous_q5",
  "P_Bachelor_q5",
  "Irsad_q5",
  "state",
  "Australia"
)

# ==============================================================================
# 1. READ FACTOR INPUT
# ==============================================================================

factor_dt <- as.data.table(
  read_parquet(
    INPUT_FILE
  )
)

cat(
  "Rows in factor input:",
  nrow(factor_dt),
  "\n"
)

cat(
  "Models:",
  uniqueN(
    factor_dt$mod
  ),
  "\n"
)

cat(
  "Simulation IDs:",
  uniqueN(
    factor_dt$sim
  ),
  "\n"
)

# ==============================================================================
# 2. GBD / SHAPLEY DECOMPOSITION FUNCTION
# ==============================================================================

decompose_pair <- function(
    dat,
    period0,
    period1
) {

  x0 <- dat[
    period == period0,
    .(
      sim,
      scenario,
      mod,
      group_id,
      level,
      group,

      P0 = pop,
      M0 = mortality_rate,
      A0 = af,
      B0 = burden
    )
  ]

  x1 <- dat[
    period == period1,
    .(
      sim,
      scenario,
      mod,
      group_id,

      P1 = pop,
      M1 = mortality_rate,
      A1 = af,
      B1 = burden
    )
  ]

  x <- merge(
    x0,
    x1,
    by = c(
      "sim",
      "scenario",
      "mod",
      "group_id"
    )
  )

  x[
    ,
    comparison :=
      paste0(
        period0,
        " -> ",
        period1
      )
  ]

  # Population contribution
  x[
    ,
    contribution_population :=
      (P1 - P0) *
      (
        M0 * A0 / 3 +
        (M1 * A0 + M0 * A1) / 6 +
        M1 * A1 / 3
      )
  ]

  # Mortality-rate contribution
  x[
    ,
    contribution_mortality :=
      (M1 - M0) *
      (
        P0 * A0 / 3 +
        (P1 * A0 + P0 * A1) / 6 +
        P1 * A1 / 3
      )
  ]

  # Pollution / AF contribution
  x[
    ,
    contribution_pollution :=
      (A1 - A0) *
      (
        P0 * M0 / 3 +
        (P1 * M0 + P0 * M1) / 6 +
        P1 * M1 / 3
      )
  ]

  x[
    ,
    total_change :=
      B1 - B0
  ]

  x[
    ,
    decomposition_sum :=
      contribution_population +
      contribution_pollution +
      contribution_mortality
  ]

  x[
    ,
    decomposition_error :=
      total_change -
      decomposition_sum
  ]

  cat(
    period0,
    " -> ",
    period1,
    " | max error = ",
    max(
      abs(
        x$decomposition_error
      ),
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )

  out <- melt(
    x,
    id.vars = c(
      "sim",
      "scenario",
      "mod",
      "group_id",
      "level",
      "group",
      "comparison",
      "B0",
      "B1",
      "total_change"
    ),
    measure.vars = c(
      "contribution_population",
      "contribution_pollution",
      "contribution_mortality"
    ),
    variable.name =
      "component",
    value.name =
      "contribution"
  )

  out[
    ,
    component := fcase(
      component ==
        "contribution_population",
      "Population",

      component ==
        "contribution_pollution",
      "Pollution (AF)",

      component ==
        "contribution_mortality",
      "Mortality rate"
    )
  ]

  out[]
}

# ==============================================================================
# 3. BASELINE -> EACH FUTURE PERIOD
# ==============================================================================

decomp_list <- vector(
  "list",
  length(
    FUTURE_PERIODS
  )
)

for (
  j in seq_along(
    FUTURE_PERIODS
  )
) {

  p1 <- FUTURE_PERIODS[j]

  cat(
    "\nDecomposing ",
    BASELINE_PERIOD,
    " -> ",
    p1,
    "\n",
    sep = ""
  )

  decomp_list[[j]] <-
    decompose_pair(
      factor_dt,
      BASELINE_PERIOD,
      p1
    )
}

decomp_dt <- rbindlist(
  decomp_list
)

rm(
  decomp_list,
  factor_dt
)

gc()

# ==============================================================================
# 4. CONTRIBUTION % FOR EVERY GCM × SIMULATION DRAW
# ==============================================================================

decomp_dt[
  ,
  contribution_pct :=
    fifelse(
      abs(total_change) > 1e-12,
      100 *
        contribution /
        total_change,
      NA_real_
    )
]

write_parquet(
  as.data.frame(decomp_dt),
  file.path(
    OUT_DIR,
    "GBD_baseline_all_model_sim_draws.parquet"
  )
)

# ==============================================================================
# 5. POINT ESTIMATE
#
# est -> average across GCMs
# ==============================================================================

point <- decomp_dt[
  sim == "est",
  .(
    estimate =
      mean(
        contribution,
        na.rm = TRUE
      ),

    baseline_burden_est =
      mean(
        B0,
        na.rm = TRUE
      ),

    target_burden_est =
      mean(
        B1,
        na.rm = TRUE
      ),

    total_change_est =
      mean(
        total_change,
        na.rm = TRUE
      ),

    n_models_est =
      uniqueN(mod)
  ),
  by = .(
    scenario,
    level,
    group,
    comparison,
    component
  )
]

point[
  ,
  estimate_pct :=
    fifelse(
      abs(total_change_est) > 1e-12,
      100 *
        estimate /
        total_change_est,
      NA_real_
    )
]

# ==============================================================================
# 6. TOTAL UNCERTAINTY
#
# No model average.
# No simulation average.
# Direct quantiles over all model × simulation draws.
# ==============================================================================

ci <- decomp_dt[
  sim != "est",
  .(
    low =
      quantile(
        contribution,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    high =
      quantile(
        contribution,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    pct_low =
      quantile(
        contribution_pct,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    pct_high =
      quantile(
        contribution_pct,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    n_draw = .N,
    n_sim = uniqueN(sim),
    n_models = uniqueN(mod)
  ),
  by = .(
    scenario,
    level,
    group,
    comparison,
    component
  )
]

# ==============================================================================
# 7. TOTAL BURDEN CHANGE CI
# ==============================================================================

change_draws <- unique(
  decomp_dt[
    ,
    .(
      sim,
      scenario,
      mod,
      level,
      group,
      comparison,
      B0,
      B1,
      total_change
    )
  ]
)

change_ci <- change_draws[
  sim != "est",
  .(
    baseline_burden_low =
      quantile(
        B0,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    baseline_burden_high =
      quantile(
        B0,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    target_burden_low =
      quantile(
        B1,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    target_burden_high =
      quantile(
        B1,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    total_change_low =
      quantile(
        total_change,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    total_change_high =
      quantile(
        total_change,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    n_total_draw = .N
  ),
  by = .(
    scenario,
    level,
    group,
    comparison
  )
]

# ==============================================================================
# 8. MERGE
# ==============================================================================

final_long <- merge(
  point,
  ci,
  by = c(
    "scenario",
    "level",
    "group",
    "comparison",
    "component"
  ),
  all.x = TRUE
)

final_long <- merge(
  final_long,
  change_ci,
  by = c(
    "scenario",
    "level",
    "group",
    "comparison"
  ),
  all.x = TRUE
)

# ==============================================================================
# 9. ORDER
# ==============================================================================

final_long[
  ,
  level :=
    factor(
      level,
      levels = LEVEL_ORDER
    )
]

final_long[
  ,
  component :=
    factor(
      component,
      levels = c(
        "Population",
        "Pollution (AF)",
        "Mortality rate"
      )
    )
]

setorder(
  final_long,
  scenario,
  level,
  group,
  comparison,
  component
)

# ==============================================================================
# 10. VERIFY EXACT POINT DECOMPOSITION
# ==============================================================================

check <- final_long[
  ,
  .(
    contribution_sum =
      sum(
        estimate
      ),

    total_change =
      unique(
        total_change_est
      )
  ),
  by = .(
    scenario,
    level,
    group,
    comparison
  )
]

check[
  ,
  error :=
    total_change -
    contribution_sum
]

cat(
  "\nMaximum point-estimate decomposition error:",
  max(
    abs(
      check$error
    ),
    na.rm = TRUE
  ),
  "\n"
)

# ==============================================================================
# 11. CI STRINGS
# ==============================================================================

final_long[
  ,
  estimate_ci :=
    sprintf(
      "%.2f (%.2f, %.2f)",
      estimate,
      low,
      high
    )
]

final_long[
  ,
  pct_ci :=
    sprintf(
      "%.1f%% (%.1f%%, %.1f%%)",
      estimate_pct,
      pct_low,
      pct_high
    )
]

# ==============================================================================
# 12. SAVE LONG
# ==============================================================================

write_csv(
  as_tibble(final_long),
  file.path(
    OUT_DIR,
    "GBD_3factor_baseline_total_uncertainty_long.csv"
  )
)

# ==============================================================================
# 13. WIDE TABLE
# ==============================================================================

final_wide <- as_tibble(
  final_long
) %>%
  mutate(
    component_short =
      case_when(
        component ==
          "Population" ~
          "population",

        component ==
          "Pollution (AF)" ~
          "pollution_af",

        component ==
          "Mortality rate" ~
          "mortality_rate"
      )
  ) %>%
  select(
    scenario,
    level,
    group,
    comparison,

    baseline_burden_est,
    baseline_burden_low,
    baseline_burden_high,

    target_burden_est,
    target_burden_low,
    target_burden_high,

    total_change_est,
    total_change_low,
    total_change_high,

    component_short,

    estimate,
    low,
    high,

    estimate_pct,
    pct_low,
    pct_high,

    estimate_ci,
    pct_ci,

    n_draw,
    n_sim,
    n_models
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from =
      component_short,
    values_from = c(
      estimate,
      low,
      high,
      estimate_pct,
      pct_low,
      pct_high,
      estimate_ci,
      pct_ci,
      n_draw,
      n_sim,
      n_models
    ),
    names_glue =
      "{component_short}_{.value}"
  )

write_csv(
  final_wide,
  file.path(
    OUT_DIR,
    "GBD_3factor_baseline_total_uncertainty_wide.csv"
  )
)

# ==============================================================================
# 14. SAVE EACH STRATIFICATION
# ==============================================================================

for (lv in LEVEL_ORDER) {

  tmp <- final_wide %>%
    filter(
      level == lv
    )

  write_csv(
    tmp,
    file.path(
      OUT_DIR,
      paste0(
        "GBD_3factor_baseline_",
        lv,
        ".csv"
      )
    )
  )
}

# ==============================================================================
# 15. AUSTRALIA
# ==============================================================================

aus_result <- final_wide %>%
  filter(
    level == "Australia"
  )

aus_result %>%
  select(
    scenario,
    comparison,

    baseline_burden_est,
    target_burden_est,
    total_change_est,

    population_estimate_ci,
    pollution_af_estimate_ci,
    mortality_rate_estimate_ci,

    population_pct_ci,
    pollution_af_pct_ci,
    mortality_rate_pct_ci
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat("\nDone.\n")
