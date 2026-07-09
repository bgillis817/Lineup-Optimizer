# ============================================================================
#  build_shape_discipline.R  - STAGE 3 (rebuilt): full BY-ZONE plate discipline
#  for BOTH batter and pitcher, per shape (AutoPitchType) per leverage bucket.
#
#  Zone = InZone flag (in-zone vs out-of-zone).
#  Metrics (both sides, so they cross cleanly):
#    z_swing   = swings at in-zone / in-zone pitches
#    o_swing   = swings at out-of-zone / out-of-zone pitches   (chase)
#    z_contact = contact on in-zone swings / in-zone swings
#    o_contact = contact on out-of-zone swings / out-of-zone swings
#    whiff     = whiffs / swings
#    zone_rate = in-zone / all   (pitcher: command; batter: what he's seen)
#    xwobacon  = contact quality on balls in play (batter side meaningful)
#
#  Thin cells shrink toward the league-vs-shape baseline; large samples kept raw.
#
#  RUN: PITCH_CACHE=data/pitches_cache.rds XS_DIR=../xStatsNECBL Rscript build_shape_discipline.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
OUT_DIR <- Sys.getenv("OUT_DIR","data"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
FULL_TRUST <- as.numeric(Sys.getenv("FULL_TRUST","50"))

px <- load_pitches()
px <- px %>% mutate(shape = AutoPitchType) %>%
  filter(!is.na(shape), shape != "", shape != "Undefined")

# xwOBACON per batted ball (optional)
xw <- tryCatch(load_xwobacon(), error=function(e){ message("xwOBACON model not found: contact quality NA"); NULL })
px$xwobacon <- NA_real_
if (!is.null(xw)) {
  bip <- px %>% filter(IsBIP, !is.na(ExitSpeed), !is.na(Angle))
  if (nrow(bip) && "PitchUID" %in% names(px)) {
    pr <- xw(bip$ExitSpeed, bip$Angle, if ("Bearing" %in% names(bip)) bip$Bearing else 0)
    bip$xwobacon <- as.numeric(pr %*% WOBA_VEC)
    px <- px %>% left_join(bip %>% select(PitchUID, xwobacon), by="PitchUID", suffix=c("",".y")) %>%
      mutate(xwobacon = coalesce(xwobacon.y, xwobacon)) %>% select(-xwobacon.y)
  }
}

# full by-zone discipline for any pitch set
disc <- function(d) {
  inz <- d$InZone; sw <- d$IsSwing
  d %>% summarise(
    n         = n(),
    zone_rate = sum(inz)/n(),
    z_swing   = sum(sw & inz)/pmax(sum(inz),1),
    o_swing   = sum(sw & !inz)/pmax(sum(!inz),1),
    z_contact = sum(IsContact & inz)/pmax(sum(sw & inz),1),
    o_contact = sum(IsContact & !inz)/pmax(sum(sw & !inz),1),
    whiff     = sum(IsWhiff)/pmax(sum(sw),1),
    chase     = sum(sw & !inz)/pmax(sum(!inz),1),
    bip       = sum(IsBIP),
    xwobacon  = mean(xwobacon[IsBIP], na.rm=TRUE),
    .groups="drop")
}
RATES <- c("zone_rate","z_swing","o_swing","z_contact","o_contact","whiff","chase","xwobacon")

shrink_to <- function(val, cell_n, lg_val, K=40) {
  shr <- (val*cell_n + lg_val*K) / (cell_n + K)
  w <- pmin(cell_n / FULL_TRUST, 1)
  ifelse(is.na(val), lg_val, w*val + (1-w)*shr)
}

build_side <- function(id_cols) {
  lg <- px %>% group_by(shape) %>% group_modify(~disc(.x)) %>% ungroup()
  cells <- px %>% group_by(across(all_of(id_cols)), shape, bucket) %>%
    group_modify(~disc(.x)) %>% ungroup()
  lg_map <- lg %>% select(shape, all_of(RATES)) %>% rename_with(~paste0("lg_",.), all_of(RATES))
  cells <- cells %>% left_join(lg_map, by="shape")
  for (r in RATES) { cells[[paste0("raw_",r)]] <- cells[[r]]
    cells[[r]] <- shrink_to(cells[[r]], cells$n, cells[[paste0("lg_",r)]]) }
  list(cells = cells %>% select(-starts_with("lg_")), league = lg)
}

batter  <- build_side(c("BatterId","Batter","BatterSide","season"))
pitcher <- build_side(c("PitcherId","Pitcher","PitcherThrows","season"))

saveRDS(list(batter=batter, pitcher=pitcher, shapes=sort(unique(px$shape))),
        file.path(OUT_DIR, "shape_discipline_object.rds"))
message("== shape_discipline_object.rds: batters ", n_distinct(batter$cells$BatterId),
        ", pitchers ", n_distinct(pitcher$cells$PitcherId), " ==")

# ---- verification: a 2026 hitter's full by-zone discipline vs shapes --------
vol <- px %>% filter(season=="2026") %>% count(BatterId, Batter, name="tot")
who <- vol %>% filter(tot>=150) %>% arrange(desc(tot)) %>% slice(1)
cat("\nHitter:", who$Batter, "(2026 pitches:", who$tot, ") \u2014 by-zone discipline vs shape:\n")
batter$cells %>% filter(BatterId==who$BatterId, season=="2026") %>%
  group_by(shape) %>% summarise(pitches=sum(n),
    `Z-sw%`=round(100*weighted.mean(z_swing,n),0),
    `O-sw%`=round(100*weighted.mean(o_swing,n),0),
    `Z-ct%`=round(100*weighted.mean(z_contact,n),0),
    `O-ct%`=round(100*weighted.mean(o_contact,n),0),
    `whiff%`=round(100*weighted.mean(whiff,n),0),
    xwOBAcon=round(weighted.mean(xwobacon,pmax(bip,1)),3), .groups="drop") %>%
  arrange(desc(pitches)) %>% print()
