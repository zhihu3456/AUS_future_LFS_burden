suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(data.table)
  library(lubridate)
  library(tsModel)
  library(future)
  library(future.apply)
  library(parallelly)
})

rm(list = ls())
gc()


# ============================================================
# 1. Settings
# ============================================================

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

ROOT <- Sys.getenv("AU_FIRE_ROOT", unset = file.path(find_project_root(), "data"))

POP_DIR <- file.path(
  ROOT,
  "population"
)

HIST_EXP_DIR <- file.path(
  POP_DIR,
  "exposure_worldpop/hist_yearly"
)

SSP_EXP_DIR <- file.path(
  POP_DIR,
  "exposure_worldpop/ssp_yearly"
)

POP_CDR_FILE <- file.path(
  POP_DIR,
  "sa2_pop_cdr_cali_2000_2100.feather"
)


# ============================================================
# Raw burden output
# ============================================================

OUT_ROOT <- file.path(
  ROOT,
  "burden/sa2_burden_by_demography"
)

OUT_DIRS <- c(
  spatial_cdr =
    file.path(
      OUT_ROOT,
      "spatial_cdr"
    )
)


# ============================================================
# SA2 yearly summary output
# ============================================================

SUMMARY_ROOT <- file.path(
  ROOT,
  "burden/burden_summary_year_by_demography"
)

SUMMARY_DIRS <- c(
  spatial_cdr =
    file.path(
      SUMMARY_ROOT,
      "spatial_cdr"
    )
)


walk(
  c(
    OUT_DIRS,
    SUMMARY_DIRS
  ),
  ~ dir.create(
    .x,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


DEMOGRAPHY_MODES <- "spatial_cdr"


SCENARIOS <- c(
  "ssp126",
  "ssp245",
  "ssp370"
)


YEARS_HIST <- 2000:2020
YEARS_FUT  <- 2021:2100
YEARS_ALL  <- 2000:2100


SIM_CHUNK <- 1001L

SAVE_ALL_SIMS <- TRUE

AMR_SCALE <- 100000


MAX_WORKERS <- 25L
RESERVE_CORES <- 2L


Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)


options(
  future.globals.maxSize =
    20 * 1024^3
)


# ============================================================
# 2. Coefficients
# ============================================================

beta_tbl <- read_feather(
  file.path(
    ROOT,
    "burden/betasim_short_long.feather"
  )
) %>%
  
  as_tibble() %>%
  
  mutate(
    sim =
      as.character(
        sim
      )
  ) %>%
  
  select(
    sim,
    beta_pm25l,
    beta_o3l,
    beta_pm25s,
    beta_o3s
  ) %>%
  
  arrange(
    sim != "est"
  ) 


# ============================================================
# 3. Historical DOY mortality distribution
# ============================================================

doy <- read_feather(
  file.path(
    POP_DIR,
    "sa2_hist_MorDoy_ratio.feather"
  )
) %>%
  
  as_tibble() %>%
  
  filter(
    doy %in%
      1:365
  ) %>%
  
  group_by(
    locationID
  ) %>%
  
  mutate(
    
    death_ratio =
      death_ratio /
      sum(
        death_ratio,
        na.rm = TRUE
      )
  ) %>%
  
  ungroup()


# ============================================================
# 4. Population + CDR
# ============================================================

scenario_map <- c(
  ssp1 = "ssp126",
  ssp2 = "ssp245",
  ssp3 = "ssp370"
)


demo <- read_feather(
  POP_CDR_FILE
) %>%
  
  as_tibble() %>%
  
  mutate(
    
    scenario =
      recode(
        as.character(
          scenario
        ),
        !!!scenario_map
      ),
    
    year =
      as.integer(
        year
      )
  ) %>%
  
  filter(
    scenario %in%
      SCENARIOS,
    year %in%
      YEARS_ALL
  )


# ============================================================
# 4. Spatial CDR population / mortality
# ============================================================

pop_cdr <- demo %>%
  
  transmute(
    
    locationID,
    
    scenario,
    
    year,
    
    demography =
      "spatial_cdr",
    
    pop,
    
    cdr_used =
      cdr_spatial,
    
    deaths_annual =
      deaths_annual_spatial
  )


rm(demo)
gc()


# ============================================================
# 5. SA2 selection
# ============================================================

sa2 <- sort(
  unique(
    pop_cdr$locationID
  )
)


# ------------------------------------------------------------
# Change this on each computer
# ------------------------------------------------------------

use_sa2 <- sa2[501:2288]

# Examples:
#
# use_sa2 <- sa2[501:1000]
# use_sa2 <- sa2[1001:1500]
# use_sa2 <- sa2[1501:2000]
# use_sa2 <- sa2[2001:length(sa2)]


# ============================================================
# 6. Calendar
# ============================================================

calendar_365 <- tibble(
  
  date =
    seq.Date(
      as.Date(
        "2000-01-01"
      ),
      as.Date(
        "2100-12-31"
      ),
      by =
        "day"
    )
) %>%
  
  filter(
    format(
      date,
      "%m-%d"
    ) !=
      "02-29"
  ) %>%
  
  mutate(
    
    year =
      year(
        date
      )
  ) %>%
  
  group_by(
    year
  ) %>%
  
  mutate(
    
    doy =
      row_number()
  ) %>%
  
  ungroup()


# ============================================================
# 7. Helpers
# ============================================================

mean_or_na <- function(x) {
  
  if (
    all(
      is.na(
        x
      )
    )
  ) {
    
    NA_real_
    
  } else {
    
    mean(
      x,
      na.rm = TRUE
    )
  }
}


# ============================================================
# 7.1 Daily deaths
# ============================================================

build_daily_deaths <- function(
    tmploc,
    pop_cdr_loc,
    doy_loc
) {
  
  out <- tidyr::crossing(
    
    scenario =
      SCENARIOS,
    
    calendar_365
    
  ) %>%
    
    mutate(
      locationID =
        tmploc
    ) %>%
    
    left_join(
      
      pop_cdr_loc,
      
      by =
        c(
          "locationID",
          "scenario",
          "year"
        )
    ) %>%
    
    left_join(
      
      doy_loc %>%
        
        select(
          locationID,
          doy,
          death_ratio
        ),
      
      by =
        c(
          "locationID",
          "doy"
        )
    ) %>%
    
    group_by(
      scenario,
      year
    ) %>%
    
    mutate(
      
      death_ratio =
        death_ratio /
        sum(
          death_ratio,
          na.rm = TRUE
        ),
      
      deaths_daily =
        deaths_annual *
        death_ratio
    ) %>%
    
    ungroup() %>%
    
    arrange(
      scenario,
      date
    ) %>%
    
    group_by(
      scenario
    ) %>%
    
    group_modify(
      ~ {
        
        .x$deathexp <- rowMeans(
          
          as.matrix(
            
            tsModel::Lag(
              .x$deaths_daily,
              -seq(
                0,
                2
              )
            )
          )
        )
        
        .x
      }
    ) %>%
    
    ungroup()
  
  
  out
}


# ============================================================
# 7.2 Historical exposure
# ============================================================

load_hist_exposure <- function(
    tmploc,
    hist_ds
) {
  
  hist_ds %>%
    
    filter(
      type ==
        "sa2",
      locationID ==
        tmploc
    ) %>%
    
    select(
      
      locationID,
      
      date,
      
      firepm25 =
        firepm25_pop,
      
      fireo3 =
        fireo3_pop
    ) %>%
    
    collect() %>%
    
    as_tibble() %>%
    
    mutate(
      
      date =
        as.Date(
          date
        ),
      
      year =
        year(
          date
        ),
      
      firepm25 =
        pmax(
          firepm25,
          0
        ),
      
      fireo3 =
        pmax(
          fireo3,
          0
        )
    ) %>%
    
    filter(
      
      year %in%
        YEARS_HIST,
      
      format(
        date,
        "%m-%d"
      ) !=
        "02-29"
    ) %>%
    
    arrange(
      date
    )
}


# ============================================================
# 7.3 Future exposure
# ============================================================

load_ssp_exposure <- function(
    tmploc,
    ssp_ds
) {
  
  ssp_ds %>%
    
    filter(
      
      type ==
        "sa2",
      
      locationID ==
        tmploc,
      
      year >=
        2021,
      
      ssp %in%
        SCENARIOS
    ) %>%
    
    select(
      
      locationID,
      
      date,
      
      year,
      
      scenario =
        ssp,
      
      mod,
      
      firepm25 =
        firepm25_pop,
      
      fireo3 =
        fireo3_pop
    ) %>%
    
    collect() %>%
    
    as_tibble() %>%
    
    mutate(
      
      date =
        as.Date(
          date
        ),
      
      year =
        as.integer(
          year
        ),
      
      scenario =
        as.character(
          scenario
        ),
      
      mod =
        as.character(
          mod
        ),
      
      firepm25 =
        pmax(
          firepm25,
          0
        ),
      
      fireo3 =
        pmax(
          fireo3,
          0
        )
    ) %>%
    
    filter(
      
      year %in%
        YEARS_FUT,
      
      format(
        date,
        "%m-%d"
      ) !=
        "02-29"
    ) %>%
    
    arrange(
      scenario,
      mod,
      date
    )
}


# ============================================================
# 7.4 Long-term burden
# ============================================================

calc_long_burden <- function(
    annual_df,
    beta_tbl,
    chunk_size =
      SIM_CHUNK
) {
  
  dt <- as.data.table(
    annual_df
  )
  
  
  n <-
    nrow(
      dt
    )
  
  
  chunks <- split(
    
    seq_len(
      nrow(
        beta_tbl
      )
    ),
    
    ceiling(
      seq_len(
        nrow(
          beta_tbl
        )
      ) /
        chunk_size
    )
  )
  
  
  rbindlist(
    
    lapply(
      
      chunks,
      
      function(idx) {
        
        k <-
          length(
            idx
          )
        
        
        af_pm25 <- 1 - exp(
          
          -outer(
            dt$firepm25,
            beta_tbl$beta_pm25l[
              idx
            ],
            "*"
          )
        )
        
        
        af_o3 <- 1 - exp(
          
          -outer(
            dt$fireo3,
            beta_tbl$beta_o3l[
              idx
            ],
            "*"
          )
        )
        
        
        an_pm25 <-
          af_pm25 *
          dt$deaths_annual
        
        
        an_o3 <-
          af_o3 *
          dt$deaths_annual
        
        
        ans <- dt[
          
          rep(
            seq_len(
              n
            ),
            times =
              k
          ),
          
          .(
            locationID,
            scenario,
            mod,
            year,
            pop,
            deaths_annual
          )
        ]
        
        
        ans[
          ,
          sim :=
            rep(
              beta_tbl$sim[
                idx
              ],
              each =
                n
            )
        ]
        
        
        ans[
          ,
          an_firepm25l :=
            as.vector(
              an_pm25
            )
        ]
        
        
        ans[
          ,
          an_fireo3l :=
            as.vector(
              an_o3
            )
        ]
        
        
        ans
      }
    ),
    
    use.names =
      TRUE
  )
}


# ============================================================
# 7.5 Short-term burden
# ============================================================

calc_short_burden <- function(
    daily_df,
    beta_tbl,
    chunk_size =
      SIM_CHUNK
) {
  
  dt <- as.data.table(
    daily_df
  )
  
  
  setorder(
    dt,
    date
  )
  
  
  yrs <- sort(
    unique(
      dt$year
    )
  )
  
  
  year_group <- factor(
    dt$year,
    levels =
      yrs
  )
  
  
  chunks <- split(
    
    seq_len(
      nrow(
        beta_tbl
      )
    ),
    
    ceiling(
      seq_len(
        nrow(
          beta_tbl
        )
      ) /
        chunk_size
    )
  )
  
  
  rbindlist(
    
    lapply(
      
      chunks,
      
      function(idx) {
        
        k <-
          length(
            idx
          )
        
        
        an_pm25 <- (
          
          1 -
            exp(
              -outer(
                dt$firepm25,
                beta_tbl$beta_pm25s[
                  idx
                ],
                "*"
              )
            )
          
        ) *
          dt$deathexp
        
        
        an_o3 <- (
          
          1 -
            exp(
              -outer(
                dt$fireo3,
                beta_tbl$beta_o3s[
                  idx
                ],
                "*"
              )
            )
          
        ) *
          dt$deathexp
        
        
        an_pm25_year <- rowsum(
          
          an_pm25,
          
          group =
            year_group,
          
          reorder =
            FALSE,
          
          na.rm =
            TRUE
        )
        
        
        an_o3_year <- rowsum(
          
          an_o3,
          
          group =
            year_group,
          
          reorder =
            FALSE,
          
          na.rm =
            TRUE
        )
        
        
        data.table(
          
          sim =
            rep(
              beta_tbl$sim[
                idx
              ],
              each =
                length(
                  yrs
                )
            ),
          
          year =
            rep(
              yrs,
              times =
                k
            ),
          
          an_firepm25s =
            as.vector(
              an_pm25_year
            ),
          
          an_fireo3s =
            as.vector(
              an_o3_year
            )
        )
      }
    ),
    
    use.names =
      TRUE
  )
}


# ============================================================
# 8. YEARLY SUMMARY HELPERS
# ============================================================

AN_COLS <- c(
  "an_pm25l",
  "an_o3l",
  "an_pm25s",
  "an_o3s",
  "an_total"
)


# ============================================================
# 8.1 Add AF and AMR
# ============================================================

add_derived_metrics <- function(x) {
  
  x <- as.data.table(
    x
  )
  
  
  for (
    suf in
    c(
      "est",
      "lci",
      "uci"
    )
  ) {
    
    for (
      p in
      c(
        "pm25l",
        "o3l",
        "pm25s",
        "o3s",
        "total"
      )
    ) {
      
      an_col <-
        paste0(
          "an_",
          p,
          "_",
          suf
        )
      
      
      af_col <-
        paste0(
          "af_",
          p,
          "_",
          suf
        )
      
      
      amr_col <-
        paste0(
          "amr_",
          p,
          "_",
          suf
        )
      
      
      x[
        ,
        (af_col) :=
          fifelse(
            deaths_annual >
              0,
            get(
              an_col
            ) /
              deaths_annual,
            NA_real_
          )
      ]
      
      
      x[
        ,
        (amr_col) :=
          fifelse(
            pop >
              0,
            get(
              an_col
            ) /
              pop *
              AMR_SCALE,
            NA_real_
          )
      ]
    }
  }
  
  
  x
}


# ============================================================
# 8.2 SA2 yearly model-specific + ensemble summary
# ============================================================

summarise_sa2_year <- function(
    raw_burden,
    demography_name
) {
  
  dt <- as.data.table(
    copy(
      raw_burden
    )
  )
  
  
  # ----------------------------------------------------------
  # Rename burden columns
  # ----------------------------------------------------------
  
  setnames(
    dt,
    
    old =
      c(
        "an_firepm25l",
        "an_fireo3l",
        "an_firepm25s",
        "an_fireo3s"
      ),
    
    new =
      c(
        "an_pm25l",
        "an_o3l",
        "an_pm25s",
        "an_o3s"
      )
  )
  
  
  # ----------------------------------------------------------
  # Total burden at simulation level
  # ----------------------------------------------------------
  
  dt[
    ,
    an_total :=
      an_pm25l +
      an_o3l +
      an_pm25s +
      an_o3s
  ]
  
  
  model_by <- c(
    "locationID",
    "scenario",
    "year",
    "mod"
  )
  
  
  ensemble_by <- c(
    "locationID",
    "scenario",
    "year"
  )
  
  
  # ==========================================================
  # A. MODEL-SPECIFIC point estimate
  # ==========================================================
  
  model_est <- dt[
    
    sim ==
      "est",
    
    c(
      
      list(
        pop =
          pop[
            1L
          ],
        
        deaths_annual =
          deaths_annual[
            1L
          ]
      ),
      
      lapply(
        .SD,
        function(z) {
          z[
            1L
          ]
        }
      )
    ),
    
    by =
      model_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    model_est,
    AN_COLS,
    paste0(
      AN_COLS,
      "_est"
    )
  )
  
  
  # ==========================================================
  # B. MODEL-SPECIFIC lower CI
  # ==========================================================
  
  model_lci <- dt[
    
    sim !=
      "est",
    
    lapply(
      .SD,
      quantile,
      probs =
        0.025,
      na.rm =
        TRUE,
      names =
        FALSE
    ),
    
    by =
      model_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    model_lci,
    AN_COLS,
    paste0(
      AN_COLS,
      "_lci"
    )
  )
  
  
  # ==========================================================
  # C. MODEL-SPECIFIC upper CI
  # ==========================================================
  
  model_uci <- dt[
    
    sim !=
      "est",
    
    lapply(
      .SD,
      quantile,
      probs =
        0.975,
      na.rm =
        TRUE,
      names =
        FALSE
    ),
    
    by =
      model_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    model_uci,
    AN_COLS,
    paste0(
      AN_COLS,
      "_uci"
    )
  )
  
  
  model_res <- Reduce(
    
    function(x, y) {
      
      merge(
        x,
        y,
        by =
          model_by,
        all =
          TRUE,
        sort =
          FALSE
      )
    },
    
    list(
      model_est,
      model_lci,
      model_uci
    )
  )
  
  
  model_res[
    ,
    estimate_scope :=
      "model_specific"
  ]
  
  
  # ==========================================================
  # D. ENSEMBLE point estimate
  #
  # Mean of model-specific EST values
  # ==========================================================
  
  ensemble_est <- dt[
    
    sim ==
      "est",
    
    c(
      
      list(
        pop =
          pop[
            1L
          ],
        
        deaths_annual =
          deaths_annual[
            1L
          ]
      ),
      
      lapply(
        .SD,
        mean,
        na.rm =
          TRUE
      )
    ),
    
    by =
      ensemble_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    ensemble_est,
    AN_COLS,
    paste0(
      AN_COLS,
      "_est"
    )
  )
  
  
  # ==========================================================
  # E. ENSEMBLE lower CI
  #
  # Pool ALL model × simulation non-est values
  # ==========================================================
  
  ensemble_lci <- dt[
    
    sim !=
      "est",
    
    lapply(
      .SD,
      quantile,
      probs =
        0.025,
      na.rm =
        TRUE,
      names =
        FALSE
    ),
    
    by =
      ensemble_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    ensemble_lci,
    AN_COLS,
    paste0(
      AN_COLS,
      "_lci"
    )
  )
  
  
  # ==========================================================
  # F. ENSEMBLE upper CI
  # ==========================================================
  
  ensemble_uci <- dt[
    
    sim !=
      "est",
    
    lapply(
      .SD,
      quantile,
      probs =
        0.975,
      na.rm =
        TRUE,
      names =
        FALSE
    ),
    
    by =
      ensemble_by,
    
    .SDcols =
      AN_COLS
  ]
  
  
  setnames(
    ensemble_uci,
    AN_COLS,
    paste0(
      AN_COLS,
      "_uci"
    )
  )
  
  
  ensemble_res <- Reduce(
    
    function(x, y) {
      
      merge(
        x,
        y,
        by =
          ensemble_by,
        all =
          TRUE,
        sort =
          FALSE
      )
    },
    
    list(
      ensemble_est,
      ensemble_lci,
      ensemble_uci
    )
  )
  
  
  ensemble_res[
    ,
    `:=`(
      mod =
        "ensemble",
      
      estimate_scope =
        "ensemble"
    )
  ]
  
  
  # ==========================================================
  # G. Combine
  # ==========================================================
  
  res <- rbindlist(
    
    list(
      model_res,
      ensemble_res
    ),
    
    use.names =
      TRUE,
    
    fill =
      TRUE
  )
  
  
  # ==========================================================
  # H. AF + AMR
  # ==========================================================
  
  res <- add_derived_metrics(
    res
  )
  
  
  res[
    ,
    demography :=
      demography_name
  ]
  
  
  setcolorder(
    
    res,
    
    c(
      
      "demography",
      
      "estimate_scope",
      
      "locationID",
      
      "scenario",
      
      "mod",
      
      "year",
      
      "pop",
      
      "deaths_annual",
      
      setdiff(
        names(
          res
        ),
        c(
          "demography",
          "estimate_scope",
          "locationID",
          "scenario",
          "mod",
          "year",
          "pop",
          "deaths_annual"
        )
      )
    )
  )
  
  
  setorder(
    res,
    estimate_scope,
    scenario,
    mod,
    year
  )
  
  
  res
}


# ============================================================
# 9. Atomic writer
# ============================================================

write_atomic <- function(
    x,
    outfile
) {
  
  tmp <- paste0(
    outfile,
    ".partial"
  )
  
  
  if (
    file.exists(
      tmp
    )
  ) {
    
    unlink(
      tmp
    )
  }
  
  
  arrow::write_parquet(
    
    x,
    
    sink =
      tmp,
    
    compression =
      "zstd"
  )
  
  
  if (
    file.exists(
      outfile
    )
  ) {
    
    unlink(
      outfile
    )
  }
  
  
  file.rename(
    tmp,
    outfile
  )
}


# ============================================================
# 10. Process one SA2
# ============================================================

process_one_location <- function(
    tmploc
) {
  
  message(
    "\nSTART | ",
    tmploc
  )
  
  
  # ==========================================================
  # Exposure
  # ==========================================================
  
  hist_ds <- arrow::open_dataset(
    HIST_EXP_DIR
  )
  
  
  ssp_ds <- arrow::open_dataset(
    SSP_EXP_DIR
  )
  
  
  hist_exp <- load_hist_exposure(
    tmploc,
    hist_ds
  )
  
  
  ssp_exp <- load_ssp_exposure(
    tmploc,
    ssp_ds
  )
  
  
  rm(
    hist_ds,
    ssp_ds
  )
  
  
  models <- sort(
    unique(
      ssp_exp$mod
    )
  )
  
  
  # ==========================================================
  # Annual exposure
  # ==========================================================
  
  hist_annual <- hist_exp %>%
    
    group_by(
      year
    ) %>%
    
    summarise(
      
      firepm25 =
        mean_or_na(
          firepm25
        ),
      
      fireo3 =
        mean_or_na(
          fireo3
        ),
      
      .groups =
        "drop"
    )
  
  
  hist_annual <- crossing(
    
    hist_annual,
    
    scenario =
      SCENARIOS,
    
    mod =
      models
  )
  
  
  fut_annual <- ssp_exp %>%
    
    group_by(
      scenario,
      mod,
      year
    ) %>%
    
    summarise(
      
      firepm25 =
        mean_or_na(
          firepm25
        ),
      
      fireo3 =
        mean_or_na(
          fireo3
        ),
      
      .groups =
        "drop"
    )
  
  
  annual_exp <- bind_rows(
    hist_annual,
    fut_annual
  ) %>%
    
    mutate(
      locationID =
        tmploc
    )
  
  
  combo <- ssp_exp %>%
    
    distinct(
      scenario,
      mod
    ) %>%
    
    arrange(
      scenario,
      mod
    )
  
  
  doy_loc <- doy %>%
    
    filter(
      locationID ==
        tmploc
    )
  
  
  # ==========================================================
  # Run the spatial_cdr demographic assumption
  # ==========================================================
  
  for (
    dm in
    DEMOGRAPHY_MODES
  ) {
    
    # --------------------------------------------------------
    # Raw burden output
    # --------------------------------------------------------
    
    outfile <- file.path(
      
      OUT_DIRS[[dm]],
      
      paste0(
        tmploc,
        ".parquet"
      )
    )
    
    
    # --------------------------------------------------------
    # Yearly summary output
    # --------------------------------------------------------
    
    summaryfile <- file.path(
      
      SUMMARY_DIRS[[dm]],
      
      paste0(
        tmploc,
        ".parquet"
      )
    )
    
    
    # ========================================================
    # If both already exist -> skip
    # ========================================================
    
    if (
      file.exists(
        outfile
      ) &&
      file.exists(
        summaryfile
      )
    ) {
      
      message(
        "SKIP | ",
        tmploc,
        " | ",
        dm
      )
      
      next
    }
    
    
    # ========================================================
    # If raw burden exists but summary does not:
    #
    # directly create summary from existing raw file
    # ========================================================
    
    if (
      file.exists(
        outfile
      ) &&
      !file.exists(
        summaryfile
      ) &&
      SAVE_ALL_SIMS
    ) {
      
      message(
        "SUMMARY ONLY | ",
        tmploc,
        " | ",
        dm
      )
      
      
      raw_existing <- arrow::read_parquet(
        outfile
      )
      
      
      summary_res <- summarise_sa2_year(
        
        raw_existing,
        
        dm
      )
      
      
      write_atomic(
        summary_res,
        summaryfile
      )
      
      
      rm(
        raw_existing,
        summary_res
      )
      
      gc()
      
      next
    }
    
    
    message(
      "RUN | ",
      tmploc,
      " | ",
      dm
    )
    
    
    # ========================================================
    # Demographic data
    # ========================================================
    
    pop_cdr_loc <- pop_cdr %>%
      
      filter(
        locationID ==
          tmploc,
        demography ==
          dm
      ) %>%
      
      select(
        locationID,
        scenario,
        year,
        pop,
        deaths_annual
      )
    
    
    # ========================================================
    # Daily deaths
    # ========================================================
    
    death_daily <- build_daily_deaths(
      
      tmploc,
      
      pop_cdr_loc,
      
      doy_loc
    )
    
    
    # ========================================================
    # LONG-TERM
    # ========================================================
    
    annual_df <- annual_exp %>%
      
      left_join(
        
        pop_cdr_loc,
        
        by =
          c(
            "locationID",
            "scenario",
            "year"
          )
      )
    
    
    long_res <- calc_long_burden(
      
      annual_df,
      
      beta_tbl
    )
    
    
    rm(
      annual_df
    )
    
    gc()
    
    
    # ========================================================
    # SHORT-TERM historical
    # ========================================================
    
    short_hist <- rbindlist(
      
      lapply(
        
        SCENARIOS,
        
        function(sc) {
          
          tmp <- hist_exp %>%
            
            select(
              date,
              year,
              firepm25,
              fireo3
            ) %>%
            
            left_join(
              
              death_daily %>%
                
                filter(
                  scenario ==
                    sc,
                  year %in%
                    YEARS_HIST
                ) %>%
                
                select(
                  date,
                  year,
                  deathexp
                ),
              
              by =
                c(
                  "date",
                  "year"
                )
            )
          
          
          calc_short_burden(
            
            tmp,
            
            beta_tbl
          ) %>%
            
            as_tibble() %>%
            
            mutate(
              
              locationID =
                tmploc,
              
              scenario =
                sc
            ) %>%
            
            crossing(
              mod =
                models
            ) %>%
            
            as.data.table()
        }
      ),
      
      use.names =
        TRUE
    )
    
    
    # ========================================================
    # SHORT-TERM future
    # ========================================================
    
    short_future <- rbindlist(
      
      lapply(
        
        seq_len(
          nrow(
            combo
          )
        ),
        
        function(ii) {
          
          sc <-
            combo$scenario[
              ii
            ]
          
          
          mm <-
            combo$mod[
              ii
            ]
          
          
          tmp <- ssp_exp %>%
            
            filter(
              scenario ==
                sc,
              mod ==
                mm
            ) %>%
            
            select(
              date,
              year,
              firepm25,
              fireo3
            ) %>%
            
            left_join(
              
              death_daily %>%
                
                filter(
                  scenario ==
                    sc,
                  year %in%
                    YEARS_FUT
                ) %>%
                
                select(
                  date,
                  year,
                  deathexp
                ),
              
              by =
                c(
                  "date",
                  "year"
                )
            )
          
          
          ans <- calc_short_burden(
            
            tmp,
            
            beta_tbl
          )
          
          
          ans[
            ,
            `:=`(
              locationID =
                tmploc,
              
              scenario =
                sc,
              
              mod =
                mm
            )
          ]
          
          
          ans
        }
      ),
      
      use.names =
        TRUE
    )
    
    
    short_res <- rbindlist(
      
      list(
        short_hist,
        short_future
      ),
      
      use.names =
        TRUE
    )
    
    
    # ========================================================
    # Combine long + short
    # ========================================================
    
    final_burden <- merge(
      
      long_res,
      
      short_res,
      
      by =
        c(
          "sim",
          "locationID",
          "scenario",
          "mod",
          "year"
        ),
      
      all =
        TRUE
    )
    
    
    setcolorder(
      
      final_burden,
      
      c(
        "sim",
        "locationID",
        "scenario",
        "mod",
        "year",
        "pop",
        "deaths_annual",
        "an_firepm25l",
        "an_fireo3l",
        "an_firepm25s",
        "an_fireo3s"
      )
    )
    
    
    # ========================================================
    # ★ Create SA2 YEARLY SUMMARY immediately
    #
    # This must happen BEFORE simulations are optionally removed.
    # ========================================================
    
    summary_res <- summarise_sa2_year(
      
      final_burden,
      
      dm
    )
    
    
    # ========================================================
    # Optional compression of raw burden
    # ========================================================
    
    if (
      !SAVE_ALL_SIMS
    ) {
      
      est <- final_burden[
        sim ==
          "est"
      ]
      
      
      ci <- final_burden[
        
        sim !=
          "est",
        
        .(
          
          an_firepm25l_low =
            quantile(
              an_firepm25l,
              0.025,
              na.rm = TRUE
            ),
          
          an_firepm25l_high =
            quantile(
              an_firepm25l,
              0.975,
              na.rm = TRUE
            ),
          
          an_fireo3l_low =
            quantile(
              an_fireo3l,
              0.025,
              na.rm = TRUE
            ),
          
          an_fireo3l_high =
            quantile(
              an_fireo3l,
              0.975,
              na.rm = TRUE
            ),
          
          an_firepm25s_low =
            quantile(
              an_firepm25s,
              0.025,
              na.rm = TRUE
            ),
          
          an_firepm25s_high =
            quantile(
              an_firepm25s,
              0.975,
              na.rm = TRUE
            ),
          
          an_fireo3s_low =
            quantile(
              an_fireo3s,
              0.025,
              na.rm = TRUE
            ),
          
          an_fireo3s_high =
            quantile(
              an_fireo3s,
              0.975,
              na.rm = TRUE
            )
        ),
        
        by =
          .(
            locationID,
            scenario,
            mod,
            year
          )
      ]
      
      
      final_burden <- merge(
        
        est,
        
        ci,
        
        by =
          c(
            "locationID",
            "scenario",
            "mod",
            "year"
          )
      )
    }
    
    
    # ========================================================
    # Save RAW burden
    # ========================================================
    
    write_atomic(
      
      final_burden,
      
      outfile
    )
    
    
    # ========================================================
    # Save YEARLY summary
    # ========================================================
    
    write_atomic(
      
      summary_res,
      
      summaryfile
    )
    
    
    message(
      "DONE | ",
      tmploc,
      " | ",
      dm
    )
    
    
    rm(
      pop_cdr_loc,
      death_daily,
      long_res,
      short_hist,
      short_future,
      short_res,
      final_burden,
      summary_res
    )
    
    
    gc()
  }
  
  
  rm(
    hist_exp,
    ssp_exp,
    hist_annual,
    fut_annual,
    annual_exp,
    combo,
    doy_loc
  )
  
  
  gc()
  
  
  tmploc
}


# ============================================================
# 11. Parallel run
# ============================================================

available_cpu <- as.integer(
  parallelly::availableCores()
)


n_workers <- min(
  
  MAX_WORKERS,
  
  max(
    1L,
    available_cpu -
      RESERVE_CORES
  ),
  
  length(
    use_sa2
  )
)


message(
  "Workers: ",
  n_workers
)


future::plan(
  
  future::multisession,
  
  workers =
    n_workers
)



status <- future.apply::future_lapply(
  
  use_sa2,
  
  function(tmploc) {
    
    tryCatch(
      
      {
        
        process_one_location(
          tmploc
        )
        
      },
      
      error =
        function(e) {
          
          message(
            "ERROR | ",
            tmploc,
            " | ",
            conditionMessage(
              e
            )
          )
          
          
          NA_character_
        }
    )
  },
  
  future.seed =
    FALSE,
  
  future.packages =
    c(
      "tidyverse",
      "arrow",
      "data.table",
      "lubridate",
      "tsModel"
    )
)


future::plan(
  future::sequential
)


gc()