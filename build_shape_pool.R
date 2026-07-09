# ============================================================================
#  build_shape_pool.R  - bakes everything the shape model needs at RUNTIME into
#  one artifact: data/shape_pool.rds
#
#  The kNN response pooling needs a pitch pool to draw neighbors from. Pulling
#  that from Drive at app start is slow and fragile (and won't fly on shinyapps).
#  So we trim the pitch table to ONLY the columns the query touches, attach the
#  xwOBACON outcome distribution, and save it. The app loads this file — no
#  Drive, no model repos, no network.
#
#  RUN (in the build): Rscript build_shape_pool.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
OUT_DIR <- Sys.getenv("OUT_DIR","data"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

DIMS <- c("RelSpeed","InducedVertBreak","HorzBreak","PlateLocHeight","PlateLocSide")
# optional cap so the artifact stays small enough to deploy
MAX_ROWS <- as.numeric(Sys.getenv("SHAPE_POOL_MAX", "300000"))

px <- load_pitches()
px <- attach_xwobacon(px)

keep <- c(DIMS, "bucket","PitcherThrows","BatterSide","BatterId","Batter",
          "PitcherId","Pitcher","season","InZone","IsSwing","IsWhiff",
          "IsContact","IsFoul","IsBIP","xwobacon",
          "xp_out","xp_1b","xp_2b","xp_3b","xp_hr")
miss <- setdiff(keep, names(px))
if (length(miss)) stop("shape pool missing columns: ", paste(miss, collapse=", "))

pool <- px %>%
  mutate(across(all_of(DIMS), ~suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(DIMS), ~!is.na(.)),
         BatterSide %in% c("Left","Right"), PitcherThrows %in% c("Left","Right")) %>%
  select(all_of(keep))

# if it's too big, subsample — but NEVER drop the current season, and keep the
# neighbor density balanced across (bucket, pitcher hand, batter side)
if (nrow(pool) > MAX_ROWS) {
  cur <- Sys.getenv("CURRENT_SEASON", "2026")
  keep_all <- pool %>% filter(season == cur)
  rest     <- pool %>% filter(season != cur)
  room <- max(0, MAX_ROWS - nrow(keep_all))
  if (room > 0 && nrow(rest) > room) {
    rest <- rest %>% group_by(bucket, PitcherThrows, BatterSide) %>%
      slice_sample(prop = room/nrow(rest)) %>% ungroup()
  }
  pool <- bind_rows(keep_all, rest)
  message("shape pool subsampled to ", nrow(pool), " rows (cap ", MAX_ROWS, ")")
}

# shrink storage: logicals stay logical, ids to character once
pool <- pool %>% mutate(across(c(InZone,IsSwing,IsWhiff,IsContact,IsFoul,IsBIP), as.logical),
                        BatterId = as.character(BatterId),
                        PitcherId = as.character(PitcherId))

# per-hitter profile for the similarity kernel (small, saved alongside)
hp <- pool %>% group_by(BatterId, Batter, BatterSide) %>% summarise(
  chase = mean(IsSwing & !InZone), whiff = sum(IsWhiff)/pmax(sum(IsSwing),1),
  xw = mean(xwobacon[IsBIP], na.rm=TRUE), pit = n(), .groups="drop")

saveRDS(list(pool = pool, hp = hp, DIMS = DIMS),
        file.path(OUT_DIR, "shape_pool.rds"))
message("== shape_pool.rds: ", nrow(pool), " pitches, ", nrow(hp), " batters (",
        round(file.size(file.path(OUT_DIR,"shape_pool.rds"))/1e6,1), " MB) ==")
