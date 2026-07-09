# ============================================================================
#  build_batter_shape.R  - STAGE 3: how each BATTER performs against each SHAPE.
#  Shape unit = AutoPitchType (movement-based, consistent with Stage 2).
#  Per (batter, shape, leverage bucket): plate discipline (whiff/chase/swing) +
#  xwOBACON (contact quality), with thin cells shrunk toward the LEAGUE-vs-that-
#  shape baseline by sample (so the number isolates this hitter's skill vs that
#  shape relative to an average hitter). Large samples keep their raw rate.
#
#  Verify: prints a known hitter's shape profile + the league baseline so you can
#  confirm his strengths/weaknesses vs shapes match what you know.
#
#  RUN: PITCH_CACHE=data/pitches_cache.rds XS_DIR=../xStatsNECBL Rscript build_batter_shape.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
OUT_DIR <- Sys.getenv("OUT_DIR","data"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
FULL_TRUST <- as.numeric(Sys.getenv("FULL_TRUST","50"))   # pitches to fully trust a shape cell

px <- load_pitches()
px <- px %>% mutate(shape = AutoPitchType) %>%
  filter(!is.na(shape), shape != "", shape != "Undefined",
         BatterSide %in% c("Left","Right"))

# xwOBACON per batted ball (optional; needs the saved model in XS_DIR)
xw <- tryCatch(load_xwobacon(), error=function(e){ message("xwOBACON model not found: contact quality skipped"); NULL })
if (!is.null(xw)) {
  bip <- px %>% filter(IsBIP, !is.na(ExitSpeed), !is.na(Angle))
  if (nrow(bip)) {
    pr <- xw(bip$ExitSpeed, bip$Angle, if ("Bearing" %in% names(bip)) bip$Bearing else 0)
    bip$xwobacon <- as.numeric(pr %*% WOBA_VEC)
    px <- px %>% left_join(bip %>% select(any_of("PitchUID"), xwobacon),
                           by = intersect("PitchUID", names(px)))
  } else px$xwobacon <- NA_real_
} else px$xwobacon <- NA_real_

# rate helper: full plate discipline split by in-zone / out-of-zone
brate <- function(d) d %>% summarise(
  n = n(),
  # zone splits
  z_pct   = sum(InZone)/n(),
  zswing  = sum(IsSwing &  InZone)/pmax(sum(InZone),1),      # swing at strikes
  oswing  = sum(IsSwing & !InZone)/pmax(sum(!InZone),1),     # chase
  zcontact= sum(IsContact & IsSwing &  InZone)/pmax(sum(IsSwing &  InZone),1),
  ocontact= sum(IsContact & IsSwing & !InZone)/pmax(sum(IsSwing & !InZone),1),
  zwhiff  = sum(IsWhiff &  InZone)/pmax(sum(IsSwing &  InZone),1),
  owhiff  = sum(IsWhiff & !InZone)/pmax(sum(IsSwing & !InZone),1),
  # overall
  swing   = sum(IsSwing)/n(),
  whiff   = sum(IsWhiff)/pmax(sum(IsSwing),1),
  chase   = oswing,
  contact = sum(IsContact)/pmax(sum(IsSwing),1),
  gb = sum(is_gb & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  fb = sum(is_fb & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  ld = sum(is_ld & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  bip     = sum(IsBIP),
  xwobacon= mean(xwobacon[IsBIP], na.rm=TRUE),
  .groups="drop")

RATES <- c("z_pct","zswing","oswing","zcontact","ocontact","zwhiff","owhiff",
           "swing","whiff","chase","contact","gb","fb","ld","xwobacon")

# league performance vs each shape (the baseline a thin batter regresses to)
lg_shape <- px %>% group_by(shape) %>% group_modify(~brate(.x)) %>% ungroup()
# batter x shape x bucket
b_shape <- px %>% group_by(BatterId, Batter, BatterSide, season, shape, bucket) %>%
  group_modify(~brate(.x)) %>% ungroup()

# shrink each cell toward league-vs-that-shape, but only for thin samples
lg_map <- lg_shape %>% select(shape, all_of(RATES)) %>% rename_with(~paste0("lg_",.), all_of(RATES))
bs <- b_shape %>% left_join(lg_map, by="shape")
shrink_to <- function(val, cell_n, lg_val, K=40) {
  raw <- val
  shr <- (val*cell_n + lg_val*K) / (cell_n + K)
  w <- pmin(cell_n / FULL_TRUST, 1)            # full raw at/above FULL_TRUST
  ifelse(is.na(raw), lg_val, w*raw + (1-w)*shr)
}
for (r in RATES) { bs[[paste0("raw_",r)]] <- bs[[r]]
  bs[[r]] <- shrink_to(bs[[r]], bs$n, bs[[paste0("lg_",r)]]) }
bs <- bs %>% select(-starts_with("lg_"))

batter_shape_object <- list(by_shape_bucket = bs, league_by_shape = lg_shape,
                            shapes = sort(unique(px$shape)))
saveRDS(batter_shape_object, file.path(OUT_DIR, "batter_shape_object.rds"))
message("== batter_shape_object.rds: ", n_distinct(bs$BatterId), " batters, ", nrow(bs), " rows ==")

# ---- verification: a high-volume 2026 hitter's shape profile vs league ------
vol <- px %>% filter(season=="2026") %>% count(BatterId, Batter, name="tot")
who <- vol %>% filter(tot >= 150) %>% arrange(desc(tot)) %>% slice(1)
cat("\nHitter:", who$Batter, " (2026 pitches:", who$tot, ")\n")
cat("(his xwOBACON & whiff by shape vs the league baseline)\n")
bs %>% filter(BatterId==who$BatterId, season=="2026") %>%
  group_by(shape) %>% summarise(pitches=sum(n),
    `Z-swing%`=round(100*weighted.mean(zswing, n),0),
    `O-swing%`=round(100*weighted.mean(oswing, n),0),
    `Z-con%`=round(100*weighted.mean(zcontact, n),0),
    `O-con%`=round(100*weighted.mean(ocontact, n),0),
    `xwOBAcon`=round(weighted.mean(xwobacon, pmax(bip,1)),3), .groups="drop") %>%
  arrange(desc(pitches)) %>% print()
cat("\nLeague baseline by shape:\n")
lg_shape %>% transmute(shape,
    `Z-swing%`=round(100*zswing,0), `O-swing%`=round(100*oswing,0),
    `Z-con%`=round(100*zcontact,0), `O-con%`=round(100*ocontact,0),
    `xwOBAcon`=round(xwobacon,3)) %>% print()
