# ============================================================================
#  build_pitcher_shape.R  - STAGE 2: attach each pitcher's QUALITY to his shapes.
#  Shape unit = AutoPitchType (movement-based, reliable; tags are ignored).
#  Per (pitcher, shape, leverage bucket): usage + whiff/chase/contact + GB/FB,
#  thin cells shrunk toward the pitcher's overall rate then league, with his
#  overall Pitching+ attached as the thin-sample backstop.
#
#  Verify: prints a nasty arm's arsenal so you can confirm his best shape grades
#  as high-whiff before Stage 3 rides on it.
#
#  RUN: PITCH_CACHE=data/pitches_cache.rds SP_DIR=../NECBLStuffPlus Rscript build_pitcher_shape.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
OUT_DIR <- Sys.getenv("OUT_DIR","data"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
SHRINK_SHAPE <- as.numeric(Sys.getenv("SHRINK_SHAPE","40"))   # pseudo-pitches toward pitcher overall
SHRINK_LG    <- as.numeric(Sys.getenv("SHRINK_LG","60"))      # then toward league

px <- load_pitches()
px <- attach_xwobacon(px)
px <- px %>% mutate(shape = AutoPitchType) %>%
  filter(!is.na(shape), shape != "", shape != "Undefined")

# complete rate set (identical to batter side): shape discipline (Z/O split),
# batted-ball (GB/FB/LD), and xwOBACON-allowed. For the pitcher these are the
# outcomes he INDUCES / ALLOWS.
rateset <- function(d) d %>% summarise(
  n = n(),
  z_pct    = sum(InZone)/n(),
  zswing   = sum(IsSwing &  InZone)/pmax(sum(InZone),1),
  oswing   = sum(IsSwing & !InZone)/pmax(sum(!InZone),1),
  zcontact = sum(IsContact & IsSwing &  InZone)/pmax(sum(IsSwing &  InZone),1),
  ocontact = sum(IsContact & IsSwing & !InZone)/pmax(sum(IsSwing & !InZone),1),
  zwhiff   = sum(IsWhiff &  InZone)/pmax(sum(IsSwing &  InZone),1),
  owhiff   = sum(IsWhiff & !InZone)/pmax(sum(IsSwing & !InZone),1),
  swing    = sum(IsSwing)/n(),
  whiff    = sum(IsWhiff)/pmax(sum(IsSwing),1),
  chase    = oswing,
  contact  = sum(IsContact)/pmax(sum(IsSwing),1),
  gb = sum(is_gb & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  fb = sum(is_fb & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  ld = sum(is_ld & IsBIP, na.rm=TRUE)/pmax(sum(IsBIP),1),
  bip = sum(IsBIP),
  xwobacon = mean(xwobacon[IsBIP], na.rm=TRUE),
  .groups="drop")

RATES <- c("z_pct","zswing","oswing","zcontact","ocontact","zwhiff","owhiff",
           "swing","whiff","chase","contact","gb","fb","ld","xwobacon")

league <- rateset(px %>% group_by(g=1)) %>% select(-g)
p_overall <- px %>% group_by(PitcherId) %>% group_modify(~rateset(.x)) %>% ungroup()
p_shape_bucket <- px %>% group_by(PitcherId, Pitcher, PitcherTeam, PitcherThrows, season, shape, bucket) %>%
  group_modify(~rateset(.x)) %>% ungroup()

# two-step shrink, but ONLY for thin samples. Shapes with >= FULL_TRUST pitches
# keep their raw rate untouched; below that, shrink toward pitcher-overall then
# league, with the pull fading in as the sample shrinks.
FULL_TRUST <- as.numeric(Sys.getenv("FULL_TRUST","60"))   # pitches to fully trust a shape
shrink <- function(cell_val, cell_n, pov_val, lg_val) {
  raw <- cell_val
  cv  <- ifelse(is.na(cell_val), lg_val, cell_val)          # NA (e.g. no BIP) -> league
  step1 <- (cv*cell_n + pov_val*SHRINK_SHAPE) / (cell_n + SHRINK_SHAPE)
  n_eff <- cell_n + SHRINK_SHAPE
  shr   <- (step1*n_eff + lg_val*SHRINK_LG) / (n_eff + SHRINK_LG)
  w <- pmin(cell_n / FULL_TRUST, 1)
  ifelse(is.na(raw), shr, w*raw + (1-w)*shr)
}

# Pitching+ (stuff) backstop: load FIRST so it can inform the thin-cell target.
sp_dir <- Sys.getenv("SP_DIR","../NECBLStuffPlus")
pp_files <- list.files(sp_dir, "necbl_pitching_plus_overall_.*\\.rds$", full.names=TRUE)
if (length(pp_files)) {
  pp <- readRDS(tail(sort(pp_files),1)) %>% mutate(PitcherId=as.character(PitcherId)) %>%
    select(PitcherId, pitching_plus)
  message("attached Pitching+ for ", nrow(pp), " pitchers")
} else { pp <- tibble(PitcherId=character(), pitching_plus=numeric())
  message("no Pitching+ file in ", sp_dir, " (stuff backstop = neutral)") }

# direction: +1 = higher Pitching+ means MORE of this rate, -1 = LESS.
PP_DIR <- c(z_pct=0, zswing=0, oswing=+1, zcontact=-1, ocontact=-1, zwhiff=+1,
            owhiff=+1, swing=0, whiff=+1, chase=+1, contact=-1, gb=+1, fb=-1,
            ld=-1, xwobacon=-1)

pov <- p_overall %>% select(PitcherId, all_of(RATES)) %>%
  rename_with(~paste0("pov_",.), all_of(RATES)) %>%
  left_join(pp, by="PitcherId") %>%
  mutate(ppf = ifelse(is.na(pitching_plus), 1, pitching_plus/100))
# adjust the shrink TARGET (his overall) toward his stuff, in the right direction
for (r in RATES) if (PP_DIR[[r]] != 0)
  pov[[paste0("pov_",r)]] <- pov[[paste0("pov_",r)]] * (pov$ppf ^ PP_DIR[[r]])
pov <- pov %>% select(PitcherId, starts_with("pov_"))

psb <- p_shape_bucket %>% left_join(pov, by="PitcherId")
for (r in RATES) psb[[paste0("raw_",r)]] <- psb[[r]]          # keep raw
for (r in RATES) psb[[r]] <- shrink(psb[[r]], psb$n, psb[[paste0("pov_",r)]], league[[r]])
psb <- psb %>% select(-starts_with("pov_")) %>% left_join(pp, by="PitcherId")

# usage per (pitcher, season, bucket): share of each shape
usage <- p_shape_bucket %>% group_by(PitcherId, season, bucket) %>%
  mutate(usage = n/sum(n)) %>% ungroup() %>%
  select(PitcherId, season, shape, bucket, usage)
psb <- psb %>% left_join(usage, by=c("PitcherId","season","shape","bucket"))

pitcher_shape_object <- list(
  by_shape_bucket = psb, overall = p_overall, league = league,
  shapes = sort(unique(px$shape)))
saveRDS(pitcher_shape_object, file.path(OUT_DIR, "pitcher_shape_object.rds"))
message("== pitcher_shape_object.rds: ",
        n_distinct(psb$PitcherId), " pitchers, ", nrow(psb), " shape-bucket rows ==")

# ---- verification: show a high-VOLUME 2026 arm's arsenal -------------------
vol26 <- px %>% filter(season=="2026") %>% count(PitcherId, Pitcher, name="tot")
p_over26 <- px %>% filter(season=="2026") %>% group_by(PitcherId) %>%
  group_modify(~rateset(.x)) %>% ungroup()
nasty <- p_over26 %>% inner_join(vol26, by="PitcherId") %>%
  filter(tot >= 150) %>% arrange(desc(whiff)) %>% slice(1)
cat("\nHigh-whiff 2026 arm (>=150 pitches):", nasty$Pitcher, " total:", nasty$tot, "\n")
cat("(usage/whiff by shape, 2026, ALL counts \u2014 out-pitch should grade highest)\n")
psb %>% filter(PitcherId==nasty$PitcherId, season=="2026") %>%
  group_by(shape) %>%
  summarise(pitches=sum(n),
            `raw whiff%`=round(100*weighted.mean(raw_whiff, n),1),
            `shrunk whiff%`=round(100*weighted.mean(whiff, n),1),
            `gb%`=round(100*weighted.mean(gb, n),1), .groups="drop") %>%
  arrange(desc(pitches)) %>% print()
