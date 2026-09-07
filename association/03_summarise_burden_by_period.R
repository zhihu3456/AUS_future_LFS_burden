# ============================================================
# Summarise SA2 burden outputs by multi-year period
# NO DUCKDB VERSION
#
# Packages:
#   tidyverse
#   arrow
#   future
#   future.apply
#   parallelly
#
# Workflow:
#
# 1. For each SA2:
#      annual data
#        -> sim x model x scenario x period
#        -> save intermediate period-level simulation file
#        -> calculate SA2 model-specific / ensemble summary
#
# 2. For grouped/state/Australia:
#      read all intermediate files using Arrow Dataset
#        -> spatial aggregation within each
#           sim x model x scenario x period x group
#        -> collect aggregated data
#        -> calculate AF / AMR
#        -> calculate model-specific / ensemble CI
#
# IMPORTANT:
# CI is calculated AFTER spatial aggregation.
# ============================================================


library(tidyverse)
library(arrow)
library(future)
library(future.apply)
library(parallelly)

rm(list = ls())
gc()


# ============================================================
# 1. Paths and settings
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

DATA_ROOT <- Sys.getenv("AU_FIRE_ROOT", unset = file.path(find_project_root(), "data"))
setwd(DATA_ROOT)

burden_dir <- "burden/sa2_burden_by_demography/spatial_cdr"

summary_dir <- paste0(
  "burden/burden_summary_period_by_demography/",
  "spatial_cdr"
)

sa2_summary_dir <- file.path(
  summary_dir,
  "sa2"
)

# Intermediate:
# one period-level sim file per SA2
period_raw_dir <- file.path(
  summary_dir,
  "period_raw"
)


# Parallel workers when processing individual SA2 files
MAX_WORKERS <- 10L
RESERVE_CORES <- 2L


AMR_SCALE <- 100000

OVERWRITE <- FALSE

REQUIRE_ALL_SA2 <- TRUE


dir.create(
  summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  sa2_summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  period_raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Period definitions
# ============================================================

period_lookup <- tibble(
  period = c(
    "2000-2020",
    "2021-2040",
    "2041-2060",
    "2061-2080",
    "2081-2100"
  ),
  period_start = c(
    2000L,
    2021L,
    2041L,
    2061L,
    2081L
  ),
  period_end = c(
    2020L,
    2040L,
    2060L,
    2080L,
    2100L
  ),
  expected_n_years = c(
    21L,
    20L,
    20L,
    20L,
    20L
  )
)


get_period <- function(year) {
  
  case_when(
    year >= 2000 & year <= 2020 ~ "2000-2020",
    year >= 2021 & year <= 2040 ~ "2021-2040",
    year >= 2041 & year <= 2060 ~ "2041-2060",
    year >= 2061 & year <= 2080 ~ "2061-2080",
    year >= 2081 & year <= 2100 ~ "2081-2100",
    TRUE ~ NA_character_
  )
}


# ============================================================
# 3. Input files
# ============================================================

files <- list.files(
  burden_dir,
  pattern = "^AUS_.*\\.parquet$",
  full.names = TRUE
)


if (length(files) == 0) {
  
  stop(
    "No SA2 burden Parquet files found."
  )
}


message(
  "SA2 burden files found: ",
  length(files)
)


# ============================================================
# 4. SA2 grouping information
# ============================================================

sa2group <- read_rds(
  "burden/sa2_used_statistics.rds"
) %>%
  dplyr::select(
    locationID,
    Urban,
    P_Indigenous_q5,
    P_Bachelor_q5,
    Irsad_q5,
    state
  ) %>%
  distinct(
    locationID,
    .keep_all = TRUE
  ) %>%
  mutate(
    across(
      c(
        Urban,
        P_Indigenous_q5,
        P_Bachelor_q5,
        Irsad_q5,
        state
      ),
      as.character
    )
  )


if (anyDuplicated(sa2group$locationID)) {
  
  stop(
    "sa2group contains duplicated locationID values."
  )
}


missing_group <- sa2group %>%
  summarise(
    across(
      c(
        Urban,
        P_Indigenous_q5,
        P_Bachelor_q5,
        Irsad_q5,
        state
      ),
      ~ sum(is.na(.x))
    )
  )


print(missing_group)


# ============================================================
# 5. Check burden files
# ============================================================

file_location_ids <- tools::file_path_sans_ext(
  basename(files)
)


missing_burden_files <- setdiff(
  sa2group$locationID,
  file_location_ids
)


extra_burden_files <- setdiff(
  file_location_ids,
  sa2group$locationID
)


message(
  "SA2 in grouping table: ",
  nrow(sa2group)
)

message(
  "SA2 burden files: ",
  length(file_location_ids)
)

message(
  "Missing burden SA2 files: ",
  length(missing_burden_files)
)

message(
  "Extra burden SA2 files: ",
  length(extra_burden_files)
)


if (
  REQUIRE_ALL_SA2 &&
  length(missing_burden_files) > 0
) {
  
  stop(
    paste0(
      "Grouped/state/Australia summaries require all SA2 files. ",
      "Missing ",
      length(missing_burden_files),
      " SA2 burden files."
    )
  )
}


# Only burden files represented in grouping table

files <- files[
  file_location_ids %in%
    sa2group$locationID
]


# ============================================================
# 6. Metric definitions
# ============================================================

an_names <- c(
  "an_pm25l",
  "an_o3l",
  "an_pm25s",
  "an_o3s",
  "an_total"
)


an_annual_mean_names <- c(
  "an_pm25l_annual_mean",
  "an_o3l_annual_mean",
  "an_pm25s_annual_mean",
  "an_o3s_annual_mean",
  "an_total_annual_mean"
)


af_names <- c(
  "af_pm25l",
  "af_o3l",
  "af_pm25s",
  "af_o3s",
  "af_total"
)


amr_names <- c(
  "amr_pm25l",
  "amr_o3l",
  "amr_pm25s",
  "amr_o3s",
  "amr_total"
)


metric_names <- c(
  an_names,
  an_annual_mean_names,
  af_names,
  amr_names
)


# ============================================================
# 7. Calculate metrics AFTER period/spatial aggregation
# ============================================================

calculate_period_metrics <- function(df) {
  
  df %>%
    mutate(
      
      # --------------------------------------------------------
      # Total AN
      # --------------------------------------------------------
      
      an_total =
        an_pm25l +
        an_o3l +
        an_pm25s +
        an_o3s,
      
      
      # --------------------------------------------------------
      # Mean annual AN
      # --------------------------------------------------------
      
      an_pm25l_annual_mean =
        if_else(
          n_years > 0,
          an_pm25l / n_years,
          NA_real_
        ),
      
      an_o3l_annual_mean =
        if_else(
          n_years > 0,
          an_o3l / n_years,
          NA_real_
        ),
      
      an_pm25s_annual_mean =
        if_else(
          n_years > 0,
          an_pm25s / n_years,
          NA_real_
        ),
      
      an_o3s_annual_mean =
        if_else(
          n_years > 0,
          an_o3s / n_years,
          NA_real_
        ),
      
      an_total_annual_mean =
        if_else(
          n_years > 0,
          an_total / n_years,
          NA_real_
        ),
      
      
      # --------------------------------------------------------
      # AF
      # --------------------------------------------------------
      
      af_pm25l =
        if_else(
          deaths_period > 0,
          an_pm25l / deaths_period,
          NA_real_
        ),
      
      af_o3l =
        if_else(
          deaths_period > 0,
          an_o3l / deaths_period,
          NA_real_
        ),
      
      af_pm25s =
        if_else(
          deaths_period > 0,
          an_pm25s / deaths_period,
          NA_real_
        ),
      
      af_o3s =
        if_else(
          deaths_period > 0,
          an_o3s / deaths_period,
          NA_real_
        ),
      
      af_total =
        if_else(
          deaths_period > 0,
          an_total / deaths_period,
          NA_real_
        ),
      
      
      # --------------------------------------------------------
      # AMR
      # --------------------------------------------------------
      
      amr_pm25l =
        if_else(
          pop_person_years > 0,
          an_pm25l /
            pop_person_years *
            AMR_SCALE,
          NA_real_
        ),
      
      amr_o3l =
        if_else(
          pop_person_years > 0,
          an_o3l /
            pop_person_years *
            AMR_SCALE,
          NA_real_
        ),
      
      amr_pm25s =
        if_else(
          pop_person_years > 0,
          an_pm25s /
            pop_person_years *
            AMR_SCALE,
          NA_real_
        ),
      
      amr_o3s =
        if_else(
          pop_person_years > 0,
          an_o3s /
            pop_person_years *
            AMR_SCALE,
          NA_real_
        ),
      
      amr_total =
        if_else(
          pop_person_years > 0,
          an_total /
            pop_person_years *
            AMR_SCALE,
          NA_real_
        )
    )
}


# ============================================================
# 8. Safe quantile
# ============================================================

safe_quantile <- function(x, prob) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      NA_real_
    )
  }
  
  as.numeric(
    quantile(
      x,
      probs = prob,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
}


# ============================================================
# 9. Model-specific summary
# ============================================================

make_model_specific_summary <- function(
    df,
    group_vars
) {
  
  df %>%
    group_by(
      across(
        all_of(group_vars)
      )
    ) %>%
    summarise(
      
      period_start =
        max(
          period_start,
          na.rm = TRUE
        ),
      
      period_end =
        max(
          period_end,
          na.rm = TRUE
        ),
      
      n_years =
        max(
          n_years,
          na.rm = TRUE
        ),
      
      pop_person_years =
        max(
          pop_person_years,
          na.rm = TRUE
        ),
      
      pop_mean_annual =
        max(
          pop_mean_annual,
          na.rm = TRUE
        ),
      
      deaths_period =
        max(
          deaths_period,
          na.rm = TRUE
        ),
      
      deaths_mean_annual =
        max(
          deaths_mean_annual,
          na.rm = TRUE
        ),
      
      
      across(
        all_of(metric_names),
        
        list(
          
          est = ~ {
            
            x <- .x[
              sim == "est"
            ]
            
            if (length(x) == 0) {
              
              NA_real_
              
            } else {
              
              max(
                x,
                na.rm = TRUE
              )
            }
          },
          
          
          lci = ~ {
            
            safe_quantile(
              .x[
                sim != "est"
              ],
              0.025
            )
          },
          
          
          uci = ~ {
            
            safe_quantile(
              .x[
                sim != "est"
              ],
              0.975
            )
          }
          
        ),
        
        .names = "{.col}_{.fn}"
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      estimate_scope =
        "model_specific",
      .before = 1
    )
}


# ============================================================
# 10. Ensemble summary
#
# est:
# mean of model-specific est values
#
# CI:
# pool all non-est model x simulation values
# ============================================================

make_ensemble_summary <- function(
    df,
    group_vars
) {
  
  df %>%
    group_by(
      across(
        all_of(group_vars)
      )
    ) %>%
    summarise(
      
      period_start =
        max(
          period_start,
          na.rm = TRUE
        ),
      
      period_end =
        max(
          period_end,
          na.rm = TRUE
        ),
      
      n_years =
        max(
          n_years,
          na.rm = TRUE
        ),
      
      pop_person_years =
        max(
          pop_person_years,
          na.rm = TRUE
        ),
      
      pop_mean_annual =
        max(
          pop_mean_annual,
          na.rm = TRUE
        ),
      
      deaths_period =
        max(
          deaths_period,
          na.rm = TRUE
        ),
      
      deaths_mean_annual =
        max(
          deaths_mean_annual,
          na.rm = TRUE
        ),
      
      
      across(
        all_of(metric_names),
        
        list(
          
          est = ~ {
            
            x <- .x[
              sim == "est"
            ]
            
            if (
              length(x) == 0 ||
              all(is.na(x))
            ) {
              
              NA_real_
              
            } else {
              
              mean(
                x,
                na.rm = TRUE
              )
            }
          },
          
          
          lci = ~ {
            
            safe_quantile(
              .x[
                sim != "est"
              ],
              0.025
            )
          },
          
          
          uci = ~ {
            
            safe_quantile(
              .x[
                sim != "est"
              ],
              0.975
            )
          }
          
        ),
        
        .names = "{.col}_{.fn}"
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      estimate_scope =
        "ensemble",
      mod =
        "ensemble",
      .before = 1
    )
}


# ============================================================
# 11. Summarise one SA2
# ============================================================

summarise_one_sa2_period <- function(
    infile
) {
  
  location_id <-
    tools::file_path_sans_ext(
      basename(infile)
    )
  
  
  sa2_outfile <- file.path(
    sa2_summary_dir,
    paste0(
      location_id,
      ".parquet"
    )
  )
  
  
  raw_outfile <- file.path(
    period_raw_dir,
    paste0(
      location_id,
      ".parquet"
    )
  )
  
  
  if (
    file.exists(sa2_outfile) &&
    file.exists(raw_outfile) &&
    !OVERWRITE
  ) {
    
    return(
      tibble(
        locationID = location_id,
        success = TRUE,
        skipped = TRUE,
        outfile = sa2_outfile,
        error = NA_character_
      )
    )
  }
  
  
  tryCatch(
    
    {
      
      # ========================================================
      # Read only required columns
      # ========================================================
      
      dat <- arrow::read_parquet(
        infile,
        col_select = c(
          sim,
          locationID,
          scenario,
          mod,
          year,
          pop,
          deaths_annual,
          an_firepm25l,
          an_fireo3l,
          an_firepm25s,
          an_fireo3s
        )
      ) %>%
        as_tibble() %>%
        filter(
          year >= 2000,
          year <= 2100
        ) %>%
        mutate(
          period =
            get_period(year)
        ) %>%
        filter(
          !is.na(period)
        )
      
      
      # ========================================================
      # Aggregate annual -> period
      #
      # retain sim/model!
      # ========================================================
      
      period_dat <- dat %>%
        group_by(
          sim,
          locationID,
          scenario,
          mod,
          period
        ) %>%
        summarise(
          
          n_years =
            n_distinct(year),
          
          pop_person_years =
            sum(
              pop,
              na.rm = TRUE
            ),
          
          deaths_period =
            sum(
              deaths_annual,
              na.rm = TRUE
            ),
          
          an_pm25l =
            sum(
              an_firepm25l,
              na.rm = TRUE
            ),
          
          an_o3l =
            sum(
              an_fireo3l,
              na.rm = TRUE
            ),
          
          an_pm25s =
            sum(
              an_firepm25s,
              na.rm = TRUE
            ),
          
          an_o3s =
            sum(
              an_fireo3s,
              na.rm = TRUE
            ),
          
          .groups = "drop"
        ) %>%
        mutate(
          
          pop_mean_annual =
            pop_person_years /
            n_years,
          
          deaths_mean_annual =
            deaths_period /
            n_years
        ) %>%
        left_join(
          period_lookup %>%
            select(
              period,
              period_start,
              period_end
            ),
          by = "period"
        )
      
      
      # ========================================================
      # Add grouping information
      #
      # This is important because Arrow can later aggregate
      # directly from period_raw files.
      # ========================================================
      
      group_info <- sa2group %>%
        filter(
          locationID ==
            location_id
        )
      
      
      if (
        nrow(group_info) != 1
      ) {
        
        stop(
          paste0(
            "Grouping information not unique for ",
            location_id
          )
        )
      }
      
      
      period_dat <- period_dat %>%
        left_join(
          group_info,
          by = "locationID"
        )
      
      
      # ========================================================
      # Save intermediate sim-level period data
      # ========================================================
      
      arrow::write_parquet(
        period_dat,
        raw_outfile,
        compression = "zstd"
      )
      
      
      # ========================================================
      # Metrics for this SA2
      # ========================================================
      
      metrics <-
        calculate_period_metrics(
          period_dat
        )
      
      
      # ========================================================
      # Model-specific
      # ========================================================
      
      model_specific <-
        make_model_specific_summary(
          metrics,
          group_vars = c(
            "locationID",
            "scenario",
            "mod",
            "period"
          )
        )
      
      
      # ========================================================
      # Ensemble
      # ========================================================
      
      ensemble <-
        make_ensemble_summary(
          metrics,
          group_vars = c(
            "locationID",
            "scenario",
            "period"
          )
        )
      
      
      # ========================================================
      # Combine
      # ========================================================
      
      result <- bind_rows(
        model_specific,
        ensemble
      ) %>%
        arrange(
          locationID,
          scenario,
          period_start,
          estimate_scope,
          mod
        )
      
      
      arrow::write_parquet(
        result,
        sa2_outfile,
        compression = "zstd"
      )
      
      
      rm(
        dat,
        period_dat,
        metrics,
        model_specific,
        ensemble,
        result
      )
      
      gc()
      
      
      tibble(
        locationID = location_id,
        success = TRUE,
        skipped = FALSE,
        outfile = sa2_outfile,
        error = NA_character_
      )
      
    },
    
    
    error = function(e) {
      
      tibble(
        locationID = location_id,
        success = FALSE,
        skipped = FALSE,
        outfile = sa2_outfile,
        error = conditionMessage(e)
      )
    }
  )
}


# ============================================================
# 12. Run SA2 files in parallel
# ============================================================

available_cpu <-
  as.integer(
    parallelly::availableCores()
  )


n_workers <- min(
  MAX_WORKERS,
  max(
    1L,
    available_cpu -
      RESERVE_CORES
  ),
  length(files)
)


message(
  "Available CPUs: ",
  available_cpu
)

message(
  "SA2 workers: ",
  n_workers
)


future::plan(
  future::multisession,
  workers = n_workers
)


sa2_status <-
  future.apply::future_lapply(
    
    files,
    
    summarise_one_sa2_period,
    
    future.seed = FALSE,
    
    future.scheduling = 1
    
  ) %>%
  bind_rows()


future::plan(
  future::sequential
)


arrow::write_parquet(
  sa2_status,
  file.path(
    summary_dir,
    "sa2_period_summary_status.parquet"
  ),
  compression = "zstd"
)


print(
  sa2_status %>%
    count(
      success,
      skipped
    )
)


if (
  any(!sa2_status$success)
) {
  
  warning(
    paste0(
      "Some SA2 period summary files failed. ",
      "Check sa2_period_summary_status.parquet"
    )
  )
}


# ============================================================
# STOP if SA2 processing failed
# ============================================================

# ============================================================
# Retry failed SA2 files sequentially
# ============================================================

failed_ids <- sa2_status %>%
  filter(!success) %>%
  pull(locationID)


if (length(failed_ids) > 0) {
  
  message(
    "\nFailed SA2 files: ",
    length(failed_ids)
  )
  
  print(failed_ids)
  
  
  # ----------------------------------------------------------
  # Match failed locationID back to original input files
  # ----------------------------------------------------------
  
  file_lookup <- setNames(
    files,
    tools::file_path_sans_ext(
      basename(files)
    )
  )
  
  
  failed_files <- file_lookup[
    failed_ids
  ]
  
  
  # ----------------------------------------------------------
  # Make sure we are no longer using parallel workers
  # ----------------------------------------------------------
  
  future::plan(
    future::sequential
  )
  
  
  # ----------------------------------------------------------
  # Retry each failed SA2
  # Up to 5 attempts
  # ----------------------------------------------------------
  
  retry_results <- vector(
    "list",
    length(failed_files)
  )
  
  
  for (i in seq_along(failed_files)) {
    
    infile <- failed_files[[i]]
    
    location_id <- names(failed_files)[i]
    
    
    message(
      "\n============================================"
    )
    
    message(
      "Retrying: ",
      location_id
    )
    
    message(
      "============================================"
    )
    
    
    # --------------------------------------------------------
    # Remove potentially incomplete outputs from failed run
    # --------------------------------------------------------
    
    sa2_outfile <- file.path(
      sa2_summary_dir,
      paste0(
        location_id,
        ".parquet"
      )
    )
    
    
    raw_outfile <- file.path(
      period_raw_dir,
      paste0(
        location_id,
        ".parquet"
      )
    )
    
    
    if (file.exists(sa2_outfile)) {
      
      file.remove(
        sa2_outfile
      )
    }
    
    
    if (file.exists(raw_outfile)) {
      
      file.remove(
        raw_outfile
      )
    }
    
    
    # --------------------------------------------------------
    # Retry loop
    # --------------------------------------------------------
    
    result <- NULL
    
    
    for (attempt in 1:5) {
      
      message(
        "Attempt ",
        attempt,
        "/5"
      )
      
      
      result <- summarise_one_sa2_period(
        infile
      )
      
      
      if (isTRUE(result$success)) {
        
        message(
          "Success: ",
          location_id
        )
        
        break
      }
      
      
      message(
        "Failed: ",
        result$error
      )
      
      
      if (attempt < 5) {
        
        wait_time <- attempt * 5
        
        message(
          "Waiting ",
          wait_time,
          " seconds before retry..."
        )
        
        Sys.sleep(
          wait_time
        )
      }
    }
    
    
    retry_results[[i]] <- result
  }
  
  
  retry_status <- bind_rows(
    retry_results
  )
  
  
  print(
    retry_status
  )
  
  
  # ----------------------------------------------------------
  # Replace failed rows in original status
  # ----------------------------------------------------------
  
  sa2_status <- sa2_status %>%
    
    filter(
      !locationID %in%
        failed_ids
    ) %>%
    
    bind_rows(
      retry_status
    )
  
  
  # ----------------------------------------------------------
  # Save updated status
  # ----------------------------------------------------------
  
  arrow::write_parquet(
    sa2_status,
    file.path(
      summary_dir,
      "sa2_period_summary_status.parquet"
    ),
    compression = "zstd"
  )
}


# ============================================================
# Final check
# ============================================================

message(
  "\nFinal SA2 processing status:"
)


print(
  sa2_status %>%
    count(
      success,
      skipped
    )
)


if (
  any(!sa2_status$success)
) {
  
  failed_final <- sa2_status %>%
    filter(
      !success
    )
  
  
  print(
    failed_final,
    width = Inf
  )
  
  
  stop(
    paste0(
      sum(!sa2_status$success),
      " SA2 file(s) still failed after retry. ",
      "Grouped aggregation stopped."
    )
  )
  
} else {
  
  message(
    "\nAll SA2 files successfully processed."
  )
  
}


gc()


# ============================================================
# 13. Open intermediate files using Arrow Dataset
#
# These files are already:
#
# sim x location x scenario x model x period
#
# Therefore dramatically smaller than original annual files.
# ============================================================

period_ds <-
  arrow::open_dataset(
    period_raw_dir,
    format = "parquet"
  )


# ============================================================
# 14. Arrow spatial aggregation function
#
# Arrow does the large spatial SUM before collect().
#
# Therefore R only receives:
#
# sim x model x scenario x period x group
#
# which is much smaller.
# ============================================================

aggregate_group_arrow <- function(
    dataset,
    level_name,
    group_variable = NULL
) {
  
  message(
    "\nAggregating: ",
    level_name
  )
  
  
  if (
    identical(
      level_name,
      "Australia"
    )
  ) {
    
    tmp <- dataset %>%
      
      group_by(
        sim,
        scenario,
        mod,
        period
      ) %>%
      
      summarise(
        
        n_years =
          max(
            n_years,
            na.rm = TRUE
          ),
        
        pop_person_years =
          sum(
            pop_person_years,
            na.rm = TRUE
          ),
        
        deaths_period =
          sum(
            deaths_period,
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
        
      ) %>%
      
      collect() %>%
      
      mutate(
        level =
          "Australia",
        
        group_value =
          "Australia"
      )
    
  } else {
    
    tmp <- dataset %>%
      
      group_by(
        sim,
        scenario,
        mod,
        period,
        group_value =
          !!rlang::sym(
            group_variable
          )
      ) %>%
      
      summarise(
        
        n_years =
          max(
            n_years,
            na.rm = TRUE
          ),
        
        pop_person_years =
          sum(
            pop_person_years,
            na.rm = TRUE
          ),
        
        deaths_period =
          sum(
            deaths_period,
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
        
      ) %>%
      
      collect() %>%
      
      mutate(
        level =
          level_name
      )
  }
  
  
  # ----------------------------------------------------------
  # Add period metadata and calculate means
  # ----------------------------------------------------------
  
  tmp <- tmp %>%
    
    mutate(
      
      pop_mean_annual =
        pop_person_years /
        n_years,
      
      deaths_mean_annual =
        deaths_period /
        n_years
      
    ) %>%
    
    left_join(
      period_lookup %>%
        select(
          period,
          period_start,
          period_end
        ),
      by = "period"
    )
  
  
  # ----------------------------------------------------------
  # AF / AMR must be calculated after spatial aggregation
  # ----------------------------------------------------------
  
  tmp <-
    calculate_period_metrics(
      tmp
    )
  
  
  tmp
}


# ============================================================
# 15. Final summary function
# ============================================================

summarise_group_final <- function(
    aggregated_data
) {
  
  
  # ----------------------------------------------------------
  # model-specific
  # ----------------------------------------------------------
  
  model_specific <-
    make_model_specific_summary(
      
      aggregated_data,
      
      group_vars = c(
        "level",
        "group_value",
        "scenario",
        "mod",
        "period"
      )
      
    )
  
  
  # ----------------------------------------------------------
  # ensemble
  # ----------------------------------------------------------
  
  ensemble <-
    make_ensemble_summary(
      
      aggregated_data,
      
      group_vars = c(
        "level",
        "group_value",
        "scenario",
        "period"
      )
      
    )
  
  
  bind_rows(
    model_specific,
    ensemble
  ) %>%
    
    arrange(
      level,
      group_value,
      scenario,
      period_start,
      estimate_scope,
      mod
    )
}


# ============================================================
# 16. State
# ============================================================

state_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name = "state",
    group_variable = "state"
  )


state_summary <-
  summarise_group_final(
    state_raw
  )


arrow::write_parquet(
  state_summary,
  file.path(
    summary_dir,
    "burden_summary_period_state.parquet"
  ),
  compression = "zstd"
)


rm(
  state_raw,
  state_summary
)

gc()


# ============================================================
# 17. Urban
# ============================================================

urban_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name = "Urban",
    group_variable = "Urban"
  )


urban_summary <-
  summarise_group_final(
    urban_raw
  )


arrow::write_parquet(
  urban_summary,
  file.path(
    summary_dir,
    "burden_summary_period_Urban.parquet"
  ),
  compression = "zstd"
)


rm(
  urban_raw,
  urban_summary
)

gc()


# ============================================================
# 18. Indigenous quintile
# ============================================================

indigenous_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name =
      "P_Indigenous_q5",
    group_variable =
      "P_Indigenous_q5"
  )


indigenous_summary <-
  summarise_group_final(
    indigenous_raw
  )


arrow::write_parquet(
  indigenous_summary,
  file.path(
    summary_dir,
    "burden_summary_period_P_Indigenous_q5.parquet"
  ),
  compression = "zstd"
)


rm(
  indigenous_raw,
  indigenous_summary
)

gc()


# ============================================================
# 19. Bachelor quintile
# ============================================================

bachelor_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name =
      "P_Bachelor_q5",
    group_variable =
      "P_Bachelor_q5"
  )


bachelor_summary <-
  summarise_group_final(
    bachelor_raw
  )


arrow::write_parquet(
  bachelor_summary,
  file.path(
    summary_dir,
    "burden_summary_period_P_Bachelor_q5.parquet"
  ),
  compression = "zstd"
)


rm(
  bachelor_raw,
  bachelor_summary
)

gc()


# ============================================================
# 20. IRSAD quintile
# ============================================================

irsad_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name =
      "Irsad_q5",
    group_variable =
      "Irsad_q5"
  )


irsad_summary <-
  summarise_group_final(
    irsad_raw
  )


arrow::write_parquet(
  irsad_summary,
  file.path(
    summary_dir,
    "burden_summary_period_Irsad_q5.parquet"
  ),
  compression = "zstd"
)


rm(
  irsad_raw,
  irsad_summary
)

gc()


# ============================================================
# 21. Australia
# ============================================================

aus_raw <-
  aggregate_group_arrow(
    period_ds,
    level_name = "Australia"
  )


aus_summary <-
  summarise_group_final(
    aus_raw
  )


arrow::write_parquet(
  aus_summary,
  file.path(
    summary_dir,
    "burden_summary_period_Australia.parquet"
  ),
  compression = "zstd"
)


rm(
  aus_raw,
  aus_summary
)

gc()


# ============================================================
# 22. Combine grouped summaries
#
# These files are small, so reading them together is fine.
# ============================================================

grouped_files <- c(
  
  file.path(
    summary_dir,
    "burden_summary_period_state.parquet"
  ),
  
  file.path(
    summary_dir,
    "burden_summary_period_Urban.parquet"
  ),
  
  file.path(
    summary_dir,
    "burden_summary_period_P_Indigenous_q5.parquet"
  ),
  
  file.path(
    summary_dir,
    "burden_summary_period_P_Bachelor_q5.parquet"
  ),
  
  file.path(
    summary_dir,
    "burden_summary_period_Irsad_q5.parquet"
  ),
  
  file.path(
    summary_dir,
    "burden_summary_period_Australia.parquet"
  )
)


grouped_summary <-
  lapply(
    grouped_files,
    arrow::read_parquet
  ) %>%
  bind_rows() %>%
  as_tibble()


arrow::write_parquet(
  grouped_summary,
  file.path(
    summary_dir,
    "burden_summary_period_grouped_all.parquet"
  ),
  compression = "zstd"
)


# ============================================================
# 23. Diagnostics
# ============================================================

message(
  "\nPeriod lookup:"
)

print(
  period_lookup
)


message(
  "\nGrouped summary dimensions:"
)

print(
  grouped_summary %>%
    count(
      level,
      estimate_scope
    )
)


# ============================================================
# Model count
# ============================================================

model_check <-
  grouped_summary %>%
  
  filter(
    estimate_scope ==
      "model_specific"
  ) %>%
  
  distinct(
    level,
    group_value,
    scenario,
    period,
    mod
  ) %>%
  
  count(
    level,
    group_value,
    scenario,
    period,
    name = "n_models"
  )


message(
  "\nModel counts:"
)

print(
  model_check %>%
    count(
      n_models
    )
)


# ============================================================
# Period-year count
# ============================================================

year_check <-
  grouped_summary %>%
  
  distinct(
    level,
    group_value,
    scenario,
    estimate_scope,
    mod,
    period,
    n_years
  ) %>%
  
  left_join(
    period_lookup,
    by = "period"
  ) %>%
  
  mutate(
    year_count_ok =
      n_years ==
      expected_n_years
  )


message(
  "\nPeriod year-count check:"
)

print(
  year_check %>%
    count(
      period,
      n_years,
      expected_n_years,
      year_count_ok
    )
)


if (
  any(
    !year_check$year_count_ok,
    na.rm = TRUE
  )
) {
  
  warning(
    paste0(
      "Some period summaries do not contain ",
      "the expected number of years. ",
      "Inspect year_check."
    )
  )
}


# ============================================================
# 24. Example
# ============================================================

# Australia ensemble
#
# aus_period <- grouped_summary %>%
#   filter(
#     level == "Australia",
#     estimate_scope == "ensemble"
#   )


# Australia model-specific
#
# aus_access <- grouped_summary %>%
#   filter(
#     level == "Australia",
#     estimate_scope == "model_specific",
#     mod == "ACCESS-CM2"
#   )


# ============================================================
# Interpretation:
#
# an_total_est
#   cumulative attributable deaths in the complete period
#
# an_total_annual_mean_est
#   mean annual attributable deaths
#
# af_total_est
#   cumulative attributable deaths /
#   cumulative deaths during the period
#
# amr_total_est
#   cumulative attributable deaths /
#   population person-years * 100000
#
# CI:
#
# model-specific
#   2.5% and 97.5% across simulations of that model
#
# ensemble
#   point estimate = mean model-specific est
#
#   CI = pooled non-est
#        model x simulation values
#
# ============================================================

message(
  "\nFinished."
)

