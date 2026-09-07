suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(data.table)
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


# ------------------------------------------------------------
# Raw SA2 burden directories
# ------------------------------------------------------------

BURDEN_ROOT <- file.path(
  ROOT,
  "burden/sa2_burden_by_demography"
)


DEMOGRAPHY_MODES <- "spatial_cdr"


BURDEN_DIRS <- c(
  spatial_cdr =
    file.path(
      BURDEN_ROOT,
      "spatial_cdr"
    )
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

SUMMARY_ROOT <- file.path(
  ROOT,
  "burden/burden_summary_year_grouped_by_demography"
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
    SUMMARY_ROOT,
    SUMMARY_DIRS
  ),
  ~ dir.create(
    .x,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


# ------------------------------------------------------------
# Other settings
# ------------------------------------------------------------

YEAR_MIN <- 2000L
YEAR_MAX <- 2100L

AMR_SCALE <- 100000

GROUP_BATCH_SIZE <- 10L


# ============================================================
# 2. SA2 -> State lookup
# ============================================================

sa2group <- read_rds(
  file.path(
    ROOT,
    "burden/sa2_used_statistics.rds"
  )
) %>%
  
  select(
    locationID,
    state
  ) %>%
  
  distinct(
    locationID,
    .keep_all = TRUE
  ) %>%
  
  mutate(
    state =
      as.character(
        state
      )
  )


# ============================================================
# 3. Column settings
# ============================================================

raw_columns <- c(
  "sim",
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


AN_COLS <- c(
  "an_pm25l",
  "an_o3l",
  "an_pm25s",
  "an_o3s",
  "an_total"
)


# ============================================================
# 4. Standardise burden names
# ============================================================

standardise_raw_burden <- function(x) {
  
  x <- as.data.table(
    x
  )
  
  
  setnames(
    
    x,
    
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
      ),
    
    skip_absent =
      TRUE
  )
  
  
  x[
    ,
    an_total :=
      an_pm25l +
      an_o3l +
      an_pm25s +
      an_o3s
  ]
  
  
  x
}


# ============================================================
# 5. AF + AMR
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
      
      an_col <- paste0(
        "an_",
        p,
        "_",
        suf
      )
      
      
      af_col <- paste0(
        "af_",
        p,
        "_",
        suf
      )
      
      
      amr_col <- paste0(
        "amr_",
        p,
        "_",
        suf
      )
      
      
      x[
        ,
        (af_col) :=
          fifelse(
            deaths_annual > 0,
            get(an_col) /
              deaths_annual,
            NA_real_
          )
      ]
      
      
      x[
        ,
        (amr_col) :=
          fifelse(
            pop > 0,
            get(an_col) /
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
# 6. Model-specific + ensemble summary
#
# Input must already be:
#
# sim × scenario × mod × year
#
# for ONE state or Australia.
# ============================================================

summarise_sim_table <- function(
    dt,
    level,
    group_value,
    demography
) {
  
  dt <- as.data.table(
    dt
  )
  
  
  # ==========================================================
  # Model-specific EST
  # ==========================================================
  
  model_est <- dt[
    
    sim ==
      "est",
    
    c(
      
      list(
        pop =
          pop[1L],
        
        deaths_annual =
          deaths_annual[1L]
      ),
      
      lapply(
        .SD,
        function(z) {
          z[1L]
        }
      )
    ),
    
    by =
      .(
        scenario,
        mod,
        year
      ),
    
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
  # Model-specific LCI
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
      .(
        scenario,
        mod,
        year
      ),
    
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
  # Model-specific UCI
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
      .(
        scenario,
        mod,
        year
      ),
    
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
          c(
            "scenario",
            "mod",
            "year"
          ),
        
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
  # Ensemble EST
  #
  # Mean of model-specific est values
  # ==========================================================
  
  ensemble_est <- dt[
    
    sim ==
      "est",
    
    c(
      
      list(
        pop =
          pop[1L],
        
        deaths_annual =
          deaths_annual[1L]
      ),
      
      lapply(
        .SD,
        mean,
        na.rm =
          TRUE
      )
    ),
    
    by =
      .(
        scenario,
        year
      ),
    
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
  # Ensemble LCI
  #
  # Pool all GCM × simulations
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
      .(
        scenario,
        year
      ),
    
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
  # Ensemble UCI
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
      .(
        scenario,
        year
      ),
    
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
          c(
            "scenario",
            "year"
          ),
        
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
  # Combine
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
  # AF / AMR
  # ==========================================================
  
  res <- add_derived_metrics(
    res
  )
  
  
  # ==========================================================
  # Add group information
  # ==========================================================
  
  res[
    ,
    `:=`(
      
      demography =
        demography,
      
      level =
        level,
      
      group_value =
        group_value
    )
  ]
  
  
  setcolorder(
    
    res,
    
    c(
      "demography",
      "level",
      "group_value",
      "estimate_scope",
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
          "level",
          "group_value",
          "estimate_scope",
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
# 7. Read + aggregate one batch
#
# SA2 files are summed at:
#
# sim × scenario × mod × year
# ============================================================

aggregate_file_batch <- function(
    batch_files
) {
  
  x <- rbindlist(
    
    lapply(
      
      batch_files,
      
      function(f) {
        
        arrow::read_parquet(
          
          f,
          
          col_select =
            raw_columns
          
        ) %>%
          
          as.data.table()
      }
    ),
    
    use.names =
      TRUE,
    
    fill =
      FALSE
  )
  
  
  x <- x[
    year >=
      YEAR_MIN &
      year <=
      YEAR_MAX
  ]
  
  
  x <- standardise_raw_burden(
    x
  )
  
  
  out <- x[
    ,
    .(
      
      pop =
        sum(
          pop,
          na.rm = TRUE
        ),
      
      deaths_annual =
        sum(
          deaths_annual,
          na.rm = TRUE
        ),
      
      an_pm25l =
        sum(
          an_pm25l,
          na.rm = TRUE
        ),
      
      an_o3l =
        sum(
          an_o3l,
          na.rm = TRUE
        ),
      
      an_pm25s =
        sum(
          an_pm25s,
          na.rm = TRUE
        ),
      
      an_o3s =
        sum(
          an_o3s,
          na.rm = TRUE
        )
      
    ),
    
    by =
      .(
        sim,
        scenario,
        mod,
        year
      )
  ]
  
  
  rm(
    x
  )
  
  gc()
  
  
  out
}


# ============================================================
# 8. Aggregate arbitrary set of SA2s
# ============================================================

aggregate_locations <- function(
    location_ids,
    file_lookup
) {
  
  files <- unname(
    file_lookup[
      location_ids
    ]
  )
  
  
  files <- files[
    !is.na(
      files
    ) &
      file.exists(
        files
      )
  ]
  
  
  batches <- split(
    
    files,
    
    ceiling(
      seq_along(
        files
      ) /
        GROUP_BATCH_SIZE
    )
  )
  
  
  accumulator <- NULL
  
  
  for (
    b in
    seq_along(
      batches
    )
  ) {
    
    tmp <- aggregate_file_batch(
      batches[[b]]
    )
    
    
    if (
      is.null(
        accumulator
      )
    ) {
      
      accumulator <-
        tmp
      
    } else {
      
      accumulator <- rbindlist(
        
        list(
          accumulator,
          tmp
        ),
        
        use.names =
          TRUE
        
      )[
        ,
        .(
          
          pop =
            sum(
              pop,
              na.rm = TRUE
            ),
          
          deaths_annual =
            sum(
              deaths_annual,
              na.rm = TRUE
            ),
          
          an_pm25l =
            sum(
              an_pm25l,
              na.rm = TRUE
            ),
          
          an_o3l =
            sum(
              an_o3l,
              na.rm = TRUE
            ),
          
          an_pm25s =
            sum(
              an_pm25s,
              na.rm = TRUE
            ),
          
          an_o3s =
            sum(
              an_o3s,
              na.rm = TRUE
            )
          
        ),
        
        by =
          .(
            sim,
            scenario,
            mod,
            year
          )
      ]
    }
    
    
    rm(
      tmp
    )
    
    gc()
  }
  
  
  accumulator[
    ,
    an_total :=
      an_pm25l +
      an_o3l +
      an_pm25s +
      an_o3s
  ]
  
  
  accumulator
}


# ============================================================
# 9. Run one demographic assumption
# ============================================================

run_one_demography <- function(
    dm
) {
  
  message(
    "\n============================================================"
  )
  
  message(
    "DEMOGRAPHY: ",
    dm
  )
  
  message(
    "============================================================"
  )
  
  
  burden_dir <-
    BURDEN_DIRS[[dm]]
  
  
  summary_dir <-
    SUMMARY_DIRS[[dm]]
  
  
  # ----------------------------------------------------------
  # Input SA2 burden files
  # ----------------------------------------------------------
  
  files <- list.files(
    
    burden_dir,
    
    pattern =
      "^AUS_.*\\.parquet$",
    
    full.names =
      TRUE
  )
  
  
  ids <- tools::file_path_sans_ext(
    basename(
      files
    )
  )
  
  
  file_lookup <- setNames(
    files,
    ids
  )
  
  
  # ==========================================================
  # 9.1 STATE summaries
  # ==========================================================
  
  states <- sort(
    unique(
      sa2group$state
    )
  )
  
  
  state_list <- lapply(
    
    states,
    
    function(st) {
      
      message(
        "State | ",
        st
      )
      
      
      location_ids <- sa2group %>%
        
        filter(
          state ==
            st
        ) %>%
        
        pull(
          locationID
        )
      
      
      agg <- aggregate_locations(
        
        location_ids =
          location_ids,
        
        file_lookup =
          file_lookup
      )
      
      
      res <- summarise_sim_table(
        
        dt =
          agg,
        
        level =
          "state",
        
        group_value =
          st,
        
        demography =
          dm
      )
      
      
      rm(
        agg
      )
      
      gc()
      
      
      res
    }
  )
  
  
  state_res <- rbindlist(
    
    state_list,
    
    use.names =
      TRUE,
    
    fill =
      TRUE
  )
  
  
  arrow::write_parquet(
    
    state_res,
    
    file.path(
      summary_dir,
      "burden_summary_year_state.parquet"
    ),
    
    compression =
      "zstd"
  )
  
  
  rm(
    state_list
  )
  
  gc()
  
  
  # ==========================================================
  # 9.2 AUSTRALIA summary
  # ==========================================================
  
  message(
    "Australia"
  )
  
  
  australia_agg <- aggregate_locations(
    
    location_ids =
      sa2group$locationID,
    
    file_lookup =
      file_lookup
  )
  
  
  australia_res <- summarise_sim_table(
    
    dt =
      australia_agg,
    
    level =
      "Australia",
    
    group_value =
      "Australia",
    
    demography =
      dm
  )
  
  
  arrow::write_parquet(
    
    australia_res,
    
    file.path(
      summary_dir,
      "burden_summary_year_Australia.parquet"
    ),
    
    compression =
      "zstd"
  )
  
  
  rm(
    australia_agg
  )
  
  gc()
  
  
  # ==========================================================
  # 9.3 State + Australia combined
  # ==========================================================
  
  combined <- rbindlist(
    
    list(
      state_res,
      australia_res
    ),
    
    use.names =
      TRUE,
    
    fill =
      TRUE
  )
  
  
  arrow::write_parquet(
    
    combined,
    
    file.path(
      summary_dir,
      "burden_summary_year_state_Australia.parquet"
    ),
    
    compression =
      "zstd"
  )
  
  
  combined
}


# ============================================================
# 10. Run all three demographic assumptions
# ============================================================

all_res <- lapply(
  
  DEMOGRAPHY_MODES,
  
  run_one_demography
)


names(
  all_res
) <-
  DEMOGRAPHY_MODES


# ============================================================
# 11. Combine ALL demographic assumptions
# ============================================================

all_demography_res <- rbindlist(
  
  all_res,
  
  use.names =
    TRUE,
  
  fill =
    TRUE
)


arrow::write_parquet(
  
  all_demography_res,
  
  file.path(
    SUMMARY_ROOT,
    "burden_summary_year_state_Australia_ALL_demography.parquet"
  ),
  
  compression =
    "zstd"
)


gc()