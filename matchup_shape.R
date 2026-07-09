# ============================================================================
#  matchup_shape.R  - STAGE 4: cross the pitcher's ACTUAL pitch cloud with each
#  batter's neighbor-pooled response (Stage 3 v2) to get per-hitter expected
#  outcomes vs that specific pitcher.
#
#  Logic: a pitcher IS his cloud of pitches (shape + location), sampled per
#  count bucket in proportion to how often he throws there. For each sampled
#  pitch we ask the batter's response (swing/chase/whiff/xwOBACON), then average
#  weighted by usage. The pitcher's quality is already encoded in the pitches
#  themselves (his stuff = shape, his command = location), so no extra scaling.
#
#  Exposes:
#    pitcher_cloud(pid, bucket, n)   -> sampled pitches (his real ones)
#    matchup_vs(batter_id, pid, ...) -> expected swing/chase/whiff/xwOBACON + xRV
#    lineup_damage(batters, pid)     -> per-hitter table for lineup construction
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
if (!exists("batter_response")) source("build_batter_response.R")   # provides px, Z, ctr, scl, hp, batter_response()

BUCKET_W <- c("0-0"=0.26,"ahead"=0.16,"even"=0.16,"behind"=0.16,
              "two_strike"=0.20,"3-2"=0.06)   # rough PA time spent per bucket

# the pitcher's real pitches in a bucket (his cloud), optionally subsampled
.CLOUD_CACHE <- new.env(parent = emptyenv())
pitcher_cloud <- function(pid, bucket, n = 40, season = NULL) {
  key <- paste(pid, bucket, n, season, sep="\u0001")
  hit <- .CLOUD_CACHE[[key]]; if (!is.null(hit)) return(hit)
  idx <- which(px$PitcherId == pid & px$bucket == bucket)
  if (!is.null(season)) idx <- idx[px$season[idx] == season]
  if (!length(idx)) idx <- which(px$PitcherId == pid)          # fallback: any count
  if (!length(idx)) return(NULL)
  if (length(idx) > n) idx <- sample(idx, n)
  out <- px[idx, c(DIMS, "PitcherThrows")]
  assign(key, out, envir = .CLOUD_CACHE)
  out
}

# batter vs pitcher: average the batter's response over the pitcher's cloud
matchup_vs <- function(batter_id, pid, season = NULL, n_per_bucket = 30,
                       buckets = names(BUCKET_W)) {
  rows <- map_dfr(buckets, function(bk) {
    cl <- pitcher_cloud(pid, bk, n_per_bucket, season)
    if (is.null(cl) || !nrow(cl)) return(tibble())
    ph <- cl$PitcherThrows[1]
    r <- map_dfr(seq_len(nrow(cl)), function(i) {
      p <- as.list(cl[i, DIMS])
      as_tibble(batter_response(batter_id, p, bucket = bk, pthrows = ph))
    })
    tibble(bucket = bk, w = BUCKET_W[[bk]],
           swing = mean(r$swing, na.rm=TRUE), chase = mean(r$chase, na.rm=TRUE),
           whiff = mean(r$whiff, na.rm=TRUE), xwobacon = mean(r$xwobacon, na.rm=TRUE))
  })
  if (!nrow(rows)) return(NULL)
  rows$w <- rows$w / sum(rows$w)
  list(by_bucket = rows,
       swing = weighted.mean(rows$swing, rows$w),
       chase = weighted.mean(rows$chase, rows$w),
       whiff = weighted.mean(rows$whiff, rows$w),
       xwobacon = weighted.mean(rows$xwobacon, rows$w, na.rm=TRUE))
}

# per-hitter expected damage vs a pitcher -> the input to lineup construction
lineup_damage <- function(batter_ids, pid, season = NULL, n_per_bucket = 25) {
  map_dfr(batter_ids, function(b) {
    m <- matchup_vs(b, pid, season, n_per_bucket)
    if (is.null(m)) return(tibble(BatterId=b))
    tibble(BatterId = b,
           Batter = hp$Batter[match(b, hp$BatterId)],
           swing = m$swing, chase = m$chase, whiff = m$whiff, xwobacon = m$xwobacon)
  }) %>% arrange(desc(xwobacon))
}

# ---- verification: a real starter vs a real lineup --------------------------
if (sys.nframe() == 0) {
  arm <- px %>% filter(season=="2026") %>% count(PitcherId, Pitcher, sort=TRUE) %>% slice(1)
  bats <- px %>% filter(season=="2026") %>% count(BatterId, Batter, sort=TRUE) %>% slice(1:9)
  cat("\nStarter:", arm$Pitcher, "\n")
  cat("Expected damage by hitter (pooled over his ACTUAL pitch cloud):\n")
  lineup_damage(bats$BatterId, arm$PitcherId, season="2026") %>%
    transmute(Batter, `swing%`=round(100*swing), `chase%`=round(100*chase),
              `whiff%`=round(100*whiff), xwOBACON=round(xwobacon,3)) %>% print()
}
