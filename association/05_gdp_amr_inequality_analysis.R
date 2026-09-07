# ==============================================================================
# GDP-RELATED AMR INEQUALITY WITH MONTE CARLO UNCERTAINTY
#
# Raw burden:
#
#   one parquet file per SA2
#
# Columns:
#
#   sim
#   locationID
#   scenario
#   mod
#   year
#   pop
#   deaths_annual
#
#   an_firepm25l
#   an_fireo3l
#   an_firepm25s
#   an_fireo3s
#
#
# Total attributable mortality:
#
#   AN =
#     an_firepm25l +
#     an_fireo3l +
#     an_firepm25s +
#     an_fireo3s
#
#
# AMR:
#
#   AMR = AN / population * 100000
#
#
# GDP ridit:
#
#   0 = lowest GDP
#   1 = highest GDP
#
#
# SII:
#
#   predicted AMR at GDP ridit 0
#   -
#   predicted AMR at GDP ridit 1
#
#   > 0 = higher AMR in low-GDP population
#
#
# RII:
#
#   predicted AMR at GDP ridit 0
#   /
#   predicted AMR at GDP ridit 1
#
#   > 1 = higher AMR in low-GDP population
#
#
# UNCERTAINTY
# ------------------------------------------------------------------------------
#
# Point estimate:
#
#   sim == "est"
#   calculate SII/RII separately for each GCM
#   then average across GCMs
#
#
# 95% CI:
#
#   sim != "est"
#
#   DO NOT average GCMs within each simulation.
#
#   Directly pool:
#
#       all GCM × simulation inequality draws
#
#   e.g.
#
#       17 GCM × 1000 simulations = 17000 draws
#
#   CI = empirical 2.5% and 97.5% quantiles
#
#
# IMPORTANT MEMORY STRATEGY
# ------------------------------------------------------------------------------
#
# We DO NOT bind all raw SA2 files.
#
# For weighted regression:
#
#       AMR = a + b * ridit
#
# weights = population
#
# sufficient statistics are:
#
#       S0 = sum(pop)
#       S1 = sum(pop * ridit)
#       S2 = sum(pop * ridit^2)
#
#       T0 = sum(pop * AMR)
#       T1 = sum(pop * ridit * AMR)
#
# Because:
#
#       pop * AMR = AN * 100000
#
# therefore:
#
#       T0 = sum(AN * 100000)
#       T1 = sum(AN * ridit * 100000)
#
# This allows us to process one SA2 at a time.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(tidyverse)
library(arrow)
library(data.table)

rm(list = ls())
gc()


# ==============================================================================
# 1. WORKING DIRECTORY
# ==============================================================================

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
# 2. SETTINGS
# ==============================================================================

BURDEN_DIR <-
  "sa2_burden_by_demography/spatial_cdr"


OUT_DIR <-
  "GDP_AMR_inequality_simulation"


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


SCENARIOS <- c(
  "ssp126",
  "ssp245",
  "ssp370"
)


SCENARIO_LABELS <- c(
  
  "ssp126" =
    "SSP1-2.6",
  
  "ssp245" =
    "SSP2-4.5",
  
  "ssp370" =
    "SSP3-7.0"
)


YEAR_MIN <- 2000L

YEAR_MAX <- 2100L


YEARS <- YEAR_MIN:YEAR_MAX


N_YEAR <-
  length(
    YEARS
  )


N_SCENARIO <-
  length(
    SCENARIOS
  )


# ------------------------------------------------------------------------------
# How often to call gc() in main loop
# ------------------------------------------------------------------------------

GC_EVERY <- 10L


# ==============================================================================
# 3. BURDEN FILES
# ==============================================================================

burden_files <- list.files(
  
  BURDEN_DIR,
  
  pattern = "\\.parquet$",
  
  full.names = TRUE
)


if (
  length(
    burden_files
  ) == 0
) {
  
  stop(
    "No burden parquet files found."
  )
}


cat(
  "\nNumber of SA2 burden files:",
  length(
    burden_files
  ),
  "\n"
)


# ==============================================================================
# 4. CHECK REQUIRED COLUMNS
# ==============================================================================

first_ds <- open_dataset(
  
  burden_files[1],
  
  format = "parquet"
)


required_cols <- c(
  
  "sim",
  "locationID",
  "scenario",
  "mod",
  "year",
  "pop",
  
  "an_firepm25l",
  "an_fireo3l",
  "an_firepm25s",
  "an_fireo3s"
)


missing_cols <- setdiff(
  
  required_cols,
  
  names(
    first_ds
  )
)


if (
  length(
    missing_cols
  ) > 0
) {
  
  stop(
    
    paste0(
      
      "Missing required columns:\n",
      
      paste(
        missing_cols,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 5. IDENTIFY REFERENCE GCM
#
# Population should be identical across GCMs and simulations.
#
# We only need one GCM to construct population-based GDP rank.
# ==============================================================================

model_check <- first_ds %>%
  
  filter(
    sim ==
      "est"
  ) %>%
  
  select(
    mod
  ) %>%
  
  distinct() %>%
  
  collect() %>%
  
  as_tibble()


MODELS <-
  sort(
    unique(
      model_check$mod
    )
  )


cat(
  "\nNumber of GCMs:",
  length(
    MODELS
  ),
  "\n"
)


print(
  MODELS
)


REF_MOD <-
  MODELS[1]


cat(
  "\nReference model used for population ranking:",
  REF_MOD,
  "\n"
)


# ==============================================================================
# 6. POPULATION QC
#
# Check one SA2:
#
# population should not vary across GCM/simulation for the same
# scenario/year.
# ==============================================================================

pop_qc_sample <- first_ds %>%
  
  filter(
    
    scenario ==
      "ssp126",
    
    year ==
      2050
  ) %>%
  
  select(
    sim,
    mod,
    pop
  ) %>%
  
  collect() %>%
  
  as_tibble()


pop_qc <- pop_qc_sample %>%
  
  summarise(
    
    n_population_values =
      n_distinct(
        pop
      ),
    
    min_pop =
      min(
        pop,
        na.rm = TRUE
      ),
    
    max_pop =
      max(
        pop,
        na.rm = TRUE
      )
  )


cat(
  "\n============================================================\n"
)

cat(
  "Population QC\n"
)

cat(
  "============================================================\n"
)


print(
  pop_qc
)


if (
  pop_qc$n_population_values != 1
) {
  
  warning(
    paste0(
      "Population differs across simulations/models in the sample. ",
      "The current inequality code assumes population is invariant ",
      "across GCM and simulation."
    )
  )
}


# ==============================================================================
# 7. READ GDP PER CAPITA
# ==============================================================================

gdp <- read_feather(
  
  "sa2_GDPpc_2000_2100_ssp123.feather"
) %>%
  
  as_tibble() %>%
  
  transmute(
    
    locationID,
    
    scenario =
      recode(
        
        scenario,
        
        "ssp1" =
          "ssp126",
        
        "ssp2" =
          "ssp245",
        
        "ssp3" =
          "ssp370"
      ),
    
    year =
      as.integer(
        year
      ),
    
    gdp_pc =
      as.numeric(
        gdp_pc
      )
  ) %>%
  
  filter(
    
    scenario %in%
      SCENARIOS,
    
    year >=
      YEAR_MIN,
    
    year <=
      YEAR_MAX
  ) %>%
  
  distinct(
    
    locationID,
    scenario,
    year,
    
    .keep_all =
      TRUE
  )


# ==============================================================================
# 8. FUNCTION:
#
# READ POPULATION FOR ONE SA2
#
# Only:
#
#   sim == est
#   one reference GCM
#
# This reduces ~5.1 million rows to ~303 rows per SA2.
# ==============================================================================

read_population_one_sa2 <- function(f) {
  
  
  ds <- open_dataset(
    
    f,
    
    format =
      "parquet"
  )
  
  
  dat <- ds %>%
    
    filter(
      
      sim ==
        "est",
      
      mod ==
        REF_MOD,
      
      scenario %in%
        SCENARIOS,
      
      year >=
        YEAR_MIN,
      
      year <=
        YEAR_MAX
    ) %>%
    
    select(
      
      locationID,
      scenario,
      year,
      pop
    ) %>%
    
    collect() %>%
    
    as_tibble() %>%
    
    transmute(
      
      locationID =
        as.character(
          locationID
        ),
      
      scenario =
        as.character(
          scenario
        ),
      
      year =
        as.integer(
          year
        ),
      
      pop =
        as.numeric(
          pop
        )
    )
  
  
  if (
    nrow(
      dat
    ) == 0
  ) {
    
    stop(
      paste0(
        "No population data found in:\n",
        f
      )
    )
  }
  
  
  locs <- unique(
    dat$locationID
  )
  
  
  if (
    length(
      locs
    ) != 1
  ) {
    
    stop(
      paste0(
        "Expected exactly one locationID in:\n",
        f
      )
    )
  }
  
  
  dat
}


# ==============================================================================
# 9. BUILD POPULATION × GDP RANK DATA
#
# First lightweight pass through all SA2 files.
# ==============================================================================

population_list <-
  vector(
    "list",
    length(
      burden_files
    )
  )


file_location <-
  character(
    length(
      burden_files
    )
  )


cat(
  "\n============================================================\n"
)

cat(
  "PASS 1: Building population-based GDP ranks\n"
)

cat(
  "============================================================\n"
)


for (
  i in seq_along(
    burden_files
  )
) {
  
  
  if (
    i == 1 ||
    i %% 100 == 0 ||
    i ==
    length(
      burden_files
    )
  ) {
    
    cat(
      "Population pass:",
      i,
      "/",
      length(
        burden_files
      ),
      "\n"
    )
  }
  
  
  tmp <-
    read_population_one_sa2(
      burden_files[i]
    )
  
  
  file_location[i] <-
    unique(
      tmp$locationID
    )[1]
  
  
  population_list[[i]] <-
    tmp
}


rank_base <-
  bind_rows(
    population_list
  )


rm(
  population_list
)

gc()


# ==============================================================================
# 10. JOIN GDP
# ==============================================================================

rank_base <- rank_base %>%
  
  left_join(
    
    gdp,
    
    by =
      c(
        "locationID",
        "scenario",
        "year"
      )
  )


rank_qc <- rank_base %>%
  
  summarise(
    
    n =
      n(),
    
    missing_gdp =
      sum(
        !is.finite(
          gdp_pc
        )
      ),
    
    missing_pop =
      sum(
        !is.finite(
          pop
        )
      )
  )


cat(
  "\n============================================================\n"
)

cat(
  "GDP / population matching QC\n"
)

cat(
  "============================================================\n"
)


print(
  rank_qc
)


if (
  rank_qc$missing_gdp > 0
) {
  
  stop(
    "Missing GDP values found. Fix GDP matching before continuing."
  )
}


if (
  rank_qc$missing_pop > 0
) {
  
  stop(
    "Missing population values found."
  )
}


# ==============================================================================
# 11. CONSTRUCT EXACT POPULATION GDP RIDIT
#
# IMPORTANT:
#
# exact GDP ties are collapsed first.
#
# All SA2s with identical GDP receive exactly the same ridit.
# ==============================================================================

rank_dt <-
  as.data.table(
    rank_base
  )


rank_dt <- rank_dt[
  
  is.finite(
    gdp_pc
  ) &
    
    is.finite(
      pop
    ) &
    
    pop > 0
]


# ------------------------------------------------------------------------------
# Total population for each exact GDP tie
# ------------------------------------------------------------------------------

tie_dt <- rank_dt[
  ,
  .(
    tie_population =
      sum(
        pop
      )
  ),
  by =
    .(
      scenario,
      year,
      gdp_pc
    )
]


setorder(
  
  tie_dt,
  
  scenario,
  year,
  gdp_pc
)


# ------------------------------------------------------------------------------
# GDP population rank
# ------------------------------------------------------------------------------

tie_dt[
  ,
  previous_population :=
    cumsum(
      tie_population
    ) -
    tie_population,
  by =
    .(
      scenario,
      year
    )
]


tie_dt[
  ,
  total_population :=
    sum(
      tie_population
    ),
  by =
    .(
      scenario,
      year
    )
]


tie_dt[
  ,
  ridit :=
    (
      previous_population +
        0.5 *
        tie_population
    ) /
    total_population
]


# ------------------------------------------------------------------------------
# Return ridit to SA2
# ------------------------------------------------------------------------------

rank_dt <- merge(
  
  rank_dt,
  
  tie_dt[
    ,
    .(
      scenario,
      year,
      gdp_pc,
      ridit
    )
  ],
  
  by =
    c(
      "scenario",
      "year",
      "gdp_pc"
    ),
  
  all.x =
    TRUE,
  
  sort =
    FALSE
)


# ==============================================================================
# 12. RANK QC
# ==============================================================================

rank_dt[
  ,
  .(
    min_ridit =
      min(
        ridit
      ),
    
    max_ridit =
      max(
        ridit
      ),
    
    total_pop =
      sum(
        pop
      )
  ),
  by =
    .(
      scenario,
      year
    )
][
  1:10
] %>%
  print()


if (
  any(
    !is.finite(
      rank_dt$ridit
    )
  )
) {
  
  stop(
    "Non-finite GDP ridit found."
  )
}


# ==============================================================================
# 13. REGRESSION CONSTANTS
#
# Weighted regression:
#
#   AMR = intercept + beta * ridit
#
# weights = population
#
#
# S0 = sum(w)
# S1 = sum(w*r)
# S2 = sum(w*r^2)
#
# These do NOT depend on GCM/simulation.
# ==============================================================================

rank_stats <- rank_dt[
  ,
  .(
    
    S0 =
      sum(
        pop
      ),
    
    S1 =
      sum(
        pop *
          ridit
      ),
    
    S2 =
      sum(
        pop *
          ridit^2
      )
  ),
  by =
    .(
      scenario,
      year
    )
]


rank_stats[
  ,
  denominator :=
    S2 -
    S1^2 /
    S0
]


if (
  any(
    !is.finite(
      rank_stats$denominator
    ) |
    rank_stats$denominator <= 0
  )
) {
  
  stop(
    "Invalid weighted regression denominator."
  )
}


# ==============================================================================
# 14. CREATE INTEGER INDEX FOR SCENARIO × YEAR
#
# Avoid expensive joins of 5 million rows inside each SA2 loop.
# ==============================================================================

make_rank_index <- function(
    scenario,
    year
) {
  
  
  scenario_index <-
    match(
      scenario,
      SCENARIOS
    )
  
  
  year_index <-
    as.integer(
      year
    ) -
    YEAR_MIN +
    1L
  
  
  (
    scenario_index -
      1L
  ) *
    N_YEAR +
    year_index
}


rank_stats[
  ,
  rank_index :=
    make_rank_index(
      scenario,
      year
    )
]


N_RANK <-
  N_SCENARIO *
  N_YEAR


S0_lookup <-
  rep(
    NA_real_,
    N_RANK
  )


S1_lookup <-
  rep(
    NA_real_,
    N_RANK
  )


S2_lookup <-
  rep(
    NA_real_,
    N_RANK
  )


DEN_lookup <-
  rep(
    NA_real_,
    N_RANK
  )


S0_lookup[
  rank_stats$rank_index
] <-
  rank_stats$S0


S1_lookup[
  rank_stats$rank_index
] <-
  rank_stats$S1


S2_lookup[
  rank_stats$rank_index
] <-
  rank_stats$S2


DEN_lookup[
  rank_stats$rank_index
] <-
  rank_stats$denominator


# ==============================================================================
# 15. LOCATION-SPECIFIC RIDIT LOOKUP
# ==============================================================================

rank_by_location <- split(
  
  rank_dt[
    ,
    .(
      scenario,
      year,
      ridit
    )
  ],
  
  rank_dt$locationID
)


# ==============================================================================
# 16. FUNCTION:
#
# READ ONE RAW BURDEN FILE
#
# Only keeps:
#
#   sim
#   scenario
#   mod
#   year
#   total attributable deaths
# ==============================================================================

read_burden_draws <- function(f) {
  
  
  ds <- open_dataset(
    
    f,
    
    format =
      "parquet"
  )
  
  
  dat <- ds %>%
    
    filter(
      
      scenario %in%
        SCENARIOS,
      
      year >=
        YEAR_MIN,
      
      year <=
        YEAR_MAX
    ) %>%
    
    transmute(
      
      sim,
      
      scenario,
      
      mod,
      
      year,
      
      
      an_total =
        
        an_firepm25l +
        
        an_fireo3l +
        
        an_firepm25s +
        
        an_fireo3s
      
    ) %>%
    
    collect()
  
  
  dt <-
    as.data.table(
      dat
    )
  
  
  dt[
    ,
    year :=
      as.integer(
        year
      )
  ]
  
  
  dt[
    ,
    sim :=
      as.character(
        sim
      )
  ]
  
  
  dt[
    ,
    scenario :=
      as.character(
        scenario
      )
  ]
  
  
  dt[
    ,
    mod :=
      as.character(
        mod
      )
  ]
  
  
  if (
    any(
      !is.finite(
        dt$an_total
      )
    )
  ) {
    
    stop(
      paste0(
        "Non-finite total burden detected in:\n",
        f
      )
    )
  }
  
  
  dt
}


# ==============================================================================
# 17. BUILD MASTER DRAW KEY
#
# Using first SA2 file.
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "Building master GCM × simulation × year key\n"
)

cat(
  "============================================================\n"
)


master_dt <-
  read_burden_draws(
    burden_files[1]
  )


master_key <- master_dt[
  ,
  .(
    sim,
    scenario,
    mod,
    year
  )
]


N_DRAW <-
  nrow(
    master_key
  )


cat(
  "\nNumber of draw rows:",
  format(
    N_DRAW,
    big.mark = ","
  ),
  "\n"
)


cat(
  "Number of simulations:",
  uniqueN(
    master_key$sim
  ),
  "\n"
)


cat(
  "Number of models:",
  uniqueN(
    master_key$mod
  ),
  "\n"
)


# ==============================================================================
# 18. MASTER SCENARIO × YEAR INDEX
# ==============================================================================

master_rank_index <-
  make_rank_index(
    
    master_key$scenario,
    
    master_key$year
  )


if (
  any(
    is.na(
      master_rank_index
    )
  )
) {
  
  stop(
    "Master draw key contains unknown scenario/year."
  )
}


# ==============================================================================
# 19. REGRESSION CONSTANTS FOR EVERY DRAW ROW
# ==============================================================================

S0_draw <-
  S0_lookup[
    master_rank_index
  ]


S1_draw <-
  S1_lookup[
    master_rank_index
  ]


S2_draw <-
  S2_lookup[
    master_rank_index
  ]


DEN_draw <-
  DEN_lookup[
    master_rank_index
  ]


if (
  any(
    !is.finite(
      S0_draw
    )
  )
) {
  
  stop(
    "Missing rank statistics for master draw rows."
  )
}


# ==============================================================================
# 20. INITIALISE SUFFICIENT-STATISTIC ACCUMULATORS
#
# T0 = sum(AN * 100000)
#
# T1 = sum(AN * ridit * 100000)
# ==============================================================================

T0 <-
  numeric(
    N_DRAW
  )


T1 <-
  numeric(
    N_DRAW
  )


# ==============================================================================
# 21. REMOVE FIRST FULL BURDEN OBJECT
# ==============================================================================

rm(
  master_dt
)

gc()


# ==============================================================================
# 22. MAIN PASS THROUGH ALL SA2 FILES
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "PASS 2: Accumulating GCM × simulation inequality statistics\n"
)

cat(
  "============================================================\n"
)


for (
  i in seq_along(
    burden_files
  )
) {
  
  
  f <-
    burden_files[i]
  
  
  loc <-
    file_location[i]
  
  
  if (
    i == 1 ||
    i %% 10 == 0 ||
    i ==
    length(
      burden_files
    )
  ) {
    
    cat(
      
      "[",
      
      format(
        Sys.time(),
        "%H:%M:%S"
      ),
      
      "] ",
      
      i,
      
      "/",
      
      length(
        burden_files
      ),
      
      " | ",
      
      loc,
      
      "\n",
      
      sep = ""
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Location ridit
  # --------------------------------------------------------------------------
  
  loc_rank <-
    rank_by_location[[loc]]
  
  
  if (
    is.null(
      loc_rank
    )
  ) {
    
    stop(
      paste0(
        "No GDP rank found for location: ",
        loc
      )
    )
  }
  
  
  loc_rank_index <-
    make_rank_index(
      
      loc_rank$scenario,
      
      loc_rank$year
    )
  
  
  loc_ridit_lookup <-
    rep(
      NA_real_,
      N_RANK
    )
  
  
  loc_ridit_lookup[
    loc_rank_index
  ] <-
    loc_rank$ridit
  
  
  ridit_draw <-
    loc_ridit_lookup[
      master_rank_index
    ]
  
  
  if (
    any(
      !is.finite(
        ridit_draw
      )
    )
  ) {
    
    stop(
      paste0(
        "Missing ridit for location: ",
        loc
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Read burden draws
  # --------------------------------------------------------------------------
  
  dt <-
    read_burden_draws(
      f
    )
  
  
  # --------------------------------------------------------------------------
  # CRITICAL:
  #
  # All location files are expected to have the same
  # sim/scenario/model/year row order.
  #
  # We explicitly verify this.
  # --------------------------------------------------------------------------
  
  same_key <-
    
    nrow(
      dt
    ) ==
    N_DRAW &&
    
    identical(
      dt$sim,
      master_key$sim
    ) &&
    
    identical(
      dt$scenario,
      master_key$scenario
    ) &&
    
    identical(
      dt$mod,
      master_key$mod
    ) &&
    
    identical(
      dt$year,
      master_key$year
    )
  
  
  if (
    !same_key
  ) {
    
    stop(
      paste0(
        "\nDraw-key order differs in:\n",
        f,
        "\n\n",
        "Do NOT accumulate this file without re-alignment."
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Population-weighted AMR sufficient statistics
  #
  # pop × AMR = AN × 100000
  # --------------------------------------------------------------------------
  
  an_scaled <-
    dt$an_total *
    100000
  
  
  T0 <-
    T0 +
    an_scaled
  
  
  T1 <-
    T1 +
    an_scaled *
    ridit_draw
  
  
  # --------------------------------------------------------------------------
  # Clean memory
  # --------------------------------------------------------------------------
  
  rm(
    dt,
    an_scaled,
    ridit_draw,
    loc_ridit_lookup,
    loc_rank
  )
  
  
  if (
    i %% GC_EVERY ==
    0
  ) {
    
    gc(
      verbose = FALSE
    )
  }
}


gc()


# ==============================================================================
# 23. CALCULATE WEIGHTED REGRESSION FOR EVERY DRAW
#
# Weighted slope:
#
# beta =
#
#   [T1 - T0*S1/S0]
#   ----------------
#   [S2 - S1^2/S0]
#
#
# intercept =
#
#   [T0 - beta*S1] / S0
#
# ==============================================================================

beta <-
  
  (
    T1 -
      T0 *
      S1_draw /
      S0_draw
  ) /
  DEN_draw


intercept <-
  
  (
    T0 -
      beta *
      S1_draw
  ) /
  S0_draw


predicted_low <-
  intercept


predicted_high <-
  intercept +
  beta


# ==============================================================================
# 24. SII
#
# low - high
#
# Since:
#
# high = intercept + beta
#
# therefore:
#
# SII = -beta
# ==============================================================================

SII <-
  predicted_low -
  predicted_high


# ==============================================================================
# 25. RII
#
# low / high
# ==============================================================================

RII <-
  
  fifelse(
    
    predicted_low > 0 &
      predicted_high > 0,
    
    predicted_low /
      predicted_high,
    
    NA_real_
  )


# ==============================================================================
# 26. NATIONAL AMR FOR QC
# ==============================================================================

mean_amr <-
  T0 /
  S0_draw


# ==============================================================================
# 27. BUILD DRAW-LEVEL RESULT
# ==============================================================================

draw_result <- copy(
  master_key
)


draw_result[
  ,
  `:=`(
    
    mean_amr =
      mean_amr,
    
    predicted_low =
      predicted_low,
    
    predicted_high =
      predicted_high,
    
    sii =
      SII,
    
    rii =
      RII
  )
]


# ==============================================================================
# 28. QC DRAW COUNTS
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "Draw-level QC\n"
)

cat(
  "============================================================\n"
)


draw_result[
  ,
  .(
    
    n =
      .N,
    
    n_models =
      uniqueN(
        mod
      ),
    
    n_sim =
      uniqueN(
        sim
      ),
    
    missing_sii =
      sum(
        !is.finite(
          sii
        )
      ),
    
    missing_rii =
      sum(
        !is.finite(
          rii
        )
      )
  ),
  by =
    scenario
] %>%
  print()


# ==============================================================================
# 29. SAVE ALL DRAW-LEVEL INEQUALITY RESULTS
# ==============================================================================

write_parquet(
  
  as.data.frame(
    draw_result
  ),
  
  file.path(
    
    OUT_DIR,
    
    "GDP_AMR_SII_RII_all_model_sim_draws.parquet"
  )
)


# ==============================================================================
# 30. POINT ESTIMATE
#
# sim == est
#
# Mean across GCMs.
# ==============================================================================

point <- draw_result[
  
  sim ==
    "est",
  
  .(
    
    amr_est =
      mean(
        mean_amr,
        na.rm = TRUE
      ),
    
    predicted_low_est =
      mean(
        predicted_low,
        na.rm = TRUE
      ),
    
    predicted_high_est =
      mean(
        predicted_high,
        na.rm = TRUE
      ),
    
    sii_est =
      mean(
        sii,
        na.rm = TRUE
      ),
    
    rii_est =
      mean(
        rii,
        na.rm = TRUE
      ),
    
    n_models_est =
      uniqueN(
        mod
      )
  ),
  
  by =
    .(
      scenario,
      year
    )
]


# ==============================================================================
# 31. 95% CI
#
# IMPORTANT:
#
# DIRECT TOTAL UNCERTAINTY
#
# Do NOT first average GCMs inside each simulation.
#
# Pool all:
#
#   simulation × GCM
#
# then calculate quantiles.
# ==============================================================================

ci <- draw_result[
  
  sim !=
    "est",
  
  .(
    
    amr_low =
      quantile(
        mean_amr,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    
    amr_high =
      quantile(
        mean_amr,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
    
    
    predicted_low_low =
      quantile(
        predicted_low,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    
    predicted_low_high =
      quantile(
        predicted_low,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
    
    
    predicted_high_low =
      quantile(
        predicted_high,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    
    predicted_high_high =
      quantile(
        predicted_high,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
    
    
    sii_low =
      quantile(
        sii,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    
    sii_high =
      quantile(
        sii,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
    
    
    rii_low =
      quantile(
        rii,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    
    rii_high =
      quantile(
        rii,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
    
    
    n_draw =
      .N,
    
    n_sim =
      uniqueN(
        sim
      ),
    
    n_models =
      uniqueN(
        mod
      )
  ),
  
  by =
    .(
      scenario,
      year
    )
]


# ==============================================================================
# 32. FINAL RESULT
# ==============================================================================

final_result <- merge(
  
  point,
  
  ci,
  
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


final_result[
  ,
  scenario_label :=
    factor(
      
      SCENARIO_LABELS[
        scenario
      ],
      
      levels =
        c(
          "SSP1-2.6",
          "SSP2-4.5",
          "SSP3-7.0"
        )
    )
]


setorder(
  
  final_result,
  
  scenario,
  year
)


# ==============================================================================
# 33. EXPECTED UNCERTAINTY DRAW COUNT
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "Final uncertainty draw counts\n"
)

cat(
  "============================================================\n"
)


final_result[
  ,
  unique(
    .(
      scenario,
      n_draw,
      n_sim,
      n_models
    )
  )
] %>%
  print()


# Expected approximately:
#
# n_draw   = 17000
# n_sim    = 1000
# n_models = 17


# ==============================================================================
# 34. FORMAT RESULT FOR TABLE
# ==============================================================================

final_result[
  ,
  `:=`(
    
    sii_text =
      sprintf(
        "%.3f (%.3f, %.3f)",
        sii_est,
        sii_low,
        sii_high
      ),
    
    rii_text =
      sprintf(
        "%.3f (%.3f, %.3f)",
        rii_est,
        rii_low,
        rii_high
      )
  )
]


# ==============================================================================
# 35. SAVE FINAL TABLE
# ==============================================================================

write_parquet(
  
  as.data.frame(
    final_result
  ),
  
  file.path(
    
    OUT_DIR,
    
    "GDP_AMR_SII_RII_annual_est_95CI.parquet"
  )
)


write_csv(
  
  as.data.frame(
    final_result
  ),
  
  file.path(
    
    OUT_DIR,
    
    "GDP_AMR_SII_RII_annual_est_95CI.csv"
  )
)


# ==============================================================================
# 36. SIMPLE MANUSCRIPT TABLE
# ==============================================================================

manuscript_table <- final_result[
  ,
  .(
    
    scenario,
    
    scenario_label,
    
    year,
    
    SII =
      sii_text,
    
    RII =
      rii_text,
    
    n_draw,
    
    n_sim,
    
    n_models
  )
]


write_csv(
  
  as.data.frame(
    manuscript_table
  ),
  
  file.path(
    
    OUT_DIR,
    
    "GDP_AMR_SII_RII_manuscript_table.csv"
  )
)
