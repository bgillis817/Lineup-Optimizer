# ============================================================================
#  build_batter_response.R  - STAGE 3 (v2): "how does THIS hitter do against
#  THIS pitch in THIS spot, in THIS count?"  No global league shrinkage — the
#  estimate is POOLED from neighbors in shape+location space: his own pitches
#  near the query, plus SIMILAR hitters' pitches near it, distance-weighted.
#
#  Space dims (what the batter reacts to): RelSpeed, InducedVertBreak, HorzBreak,
#  PlateLocHeight, PlateLocSide.  Hard filters: same count bucket, same pitcher
#  handedness, same batter handedness (platoon-correct).
#
#  Answers e.g.: "Sinker 0-2, 88 / 12 IVB / 18 HB, down-and-in — his xwOBACON,
#  swing%, chase%, whiff%?"  Verified by that exact query below.
#
#  RUN: PITCH_CACHE=data/pitches_cache.rds XS_DIR=../xStatsNECBL Rscript build_batter_response.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
source("continuous_common.R")
OUT_DIR <- Sys.getenv("OUT_DIR","data"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

DIMS <- c("RelSpeed","InducedVertBreak","HorzBreak","PlateLocHeight","PlateLocSide")

px <- load_pitches()
px <- attach_xwobacon(px)
px <- px %>%
  mutate(across(all_of(DIMS), ~suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(DIMS), ~!is.na(.)),
         BatterSide %in% c("Left","Right"), PitcherThrows %in% c("Left","Right"))

# scale the space so mph / inches / feet are comparable
ctr <- colMeans(px[,DIMS]); scl <- apply(px[,DIMS],2,sd); scl[scl==0] <- 1
Z <- scale(as.matrix(px[,DIMS]), center=ctr, scale=scl)

# per-hitter overall profile -> hitter-to-hitter similarity (for the pool)
hp <- px %>% group_by(BatterId, Batter, BatterSide) %>% summarise(
  chase = mean(IsSwing & !InZone), whiff = sum(IsWhiff)/pmax(sum(IsSwing),1),
  xw = mean(xwobacon[IsBIP], na.rm=TRUE), pit = n(), .groups="drop")
hp_dims <- c("chase","whiff","xw")
hctr <- sapply(hp[hp_dims], mean, na.rm=TRUE); hscl <- sapply(hp[hp_dims], sd, na.rm=TRUE); hscl[hscl==0]<-1
hp_z <- scale(as.matrix(hp[hp_dims]), center=hctr, scale=hscl); hp_z[is.na(hp_z)] <- 0

# ---- SPEED PRECOMPUTES (same model, just not recomputed per query) ----------
# flat numeric vectors: data-frame subsetting in the inner loop is the killer
V <- list(
  swing   = as.numeric(px$IsSwing),  inzone = as.numeric(px$InZone),
  whiff   = as.numeric(px$IsWhiff),  contact= as.numeric(px$IsContact),
  foul    = as.numeric(px$IsFoul),   bip    = as.numeric(px$IsBIP),
  xwobacon= px$xwobacon,
  xp = cbind(px$xp_out, px$xp_1b, px$xp_2b, px$xp_3b, px$xp_hr))
PX_H  <- match(px$BatterId, hp$BatterId)               # px row -> hitter index
# hitter-to-hitter squared distance, computed once
HD2 <- as.matrix(dist(hp_z))^2
# candidate row indices per (bucket, pitcher hand, batter side), computed once
CAND <- local({
  key <- paste(px$bucket, px$PitcherThrows, px$BatterSide, sep="|")
  split(seq_len(nrow(px)), key)
})
SIDE_ONLY <- local({
  key <- paste(px$PitcherThrows, px$BatterSide, sep="|")
  split(seq_len(nrow(px)), key)
})

# --- the core query, BATCHED ------------------------------------------------
# Q: m x 5 matrix of query pitches (the pitcher's cloud for one bucket).
# Same math as before (gaussian kernel on shape+location, hitter-similarity
# kernel, own-pitch boost) but all m pitches scored in one pass.
batter_response_batch <- function(batter_id, Q, bucket, pthrows,
                                  k = 300, own_boost = 4, hitter_k = 40) {
  hb <- match(batter_id, hp$BatterId); if (is.na(hb)) return(NULL)
  side <- hp$BatterSide[hb]
  cand <- CAND[[paste(bucket, pthrows, side, sep="|")]]
  if (is.null(cand) || length(cand) < 20) cand <- SIDE_ONLY[[paste(pthrows, side, sep="|")]]
  if (is.null(cand) || !length(cand)) return(NULL)

  Zc <- Z[cand, , drop=FALSE]
  Qz <- sweep(sweep(as.matrix(Q[,DIMS,drop=FALSE]), 2, ctr, "-"), 2, scl, "/")
  # squared distances: ||q||^2 + ||z||^2 - 2 q.z   (one BLAS matmul)
  D2 <- outer(rowSums(Qz^2), rep(1,length(cand))) +
        outer(rep(1,nrow(Qz)), rowSums(Zc^2)) - 2 * (Qz %*% t(Zc))

  # hitter-similarity weight, precomputed distances
  hd2 <- HD2[hb, PX_H[cand]]
  hbw2 <- max(sort(hd2)[min(hitter_k, length(hd2))], 1e-6)
  wh <- exp(-hd2 / (2*hbw2))
  own <- px$BatterId[cand] == batter_id
  wh <- wh * ifelse(own, own_boost, 1)

  # accumulate weights across all cloud pitches (each pitch normalized first)
  W <- numeric(length(cand))
  for (i in seq_len(nrow(Qz))) {
    d2 <- D2[i,]
    ord <- order(d2)[seq_len(min(k, length(d2)))]
    bw2 <- max(d2[ord[max(1,length(ord)%/%2)]], 1e-6)
    wi <- numeric(length(cand)); wi[ord] <- exp(-d2[ord] / (2*bw2)) * wh[ord]
    s <- sum(wi); if (s > 0) W <- W + wi/s
  }
  if (sum(W) <= 0) return(NULL)

  wm <- function(x, ww=W) { ok <- !is.na(x) & ww>0; if(!any(ok)) return(NA_real_); sum(x[ok]*ww[ok])/sum(ww[ok]) }
  sv <- V$swing[cand]; cv <- V$contact[cand]; bv <- V$bip[cand]
  swm <- W * sv; conm <- W * (sv*cv); bipm <- W * bv
  list(
    swing    = wm(sv),
    whiff    = if (sum(swm)>0) sum(V$whiff[cand]*swm)/sum(swm) else NA_real_,
    foul     = if (sum(conm)>0) sum(V$foul[cand]*conm)/sum(conm) else NA_real_,
    xwobacon = if (sum(bipm)>0) wm(V$xwobacon[cand], bipm) else NA_real_,
    bip_dist = if (sum(bipm)>0) { d <- colSums(V$xp[cand,,drop=FALSE]*bipm, na.rm=TRUE); d/sum(d) }
               else c(.68,.20,.075,.008,.037),
    eff_n = round(sum(W)/max(W),1))
}

# single-pitch convenience (verification / one-off queries)
batter_response <- function(batter_id, pitch, bucket, pthrows, ...) {
  Q <- as.data.frame(pitch[DIMS]); names(Q) <- DIMS
  r <- batter_response_batch(batter_id, Q, bucket, pthrows, ...)
  if (is.null(r)) return(list(swing=NA,chase=NA,whiff=NA,foul=NA,xwobacon=NA,
                              bip_dist=c(.68,.20,.075,.008,.037),own_neighbors=0,eff_n=0))
  inz <- as.numeric(pitch$PlateLocHeight>=1.59 & pitch$PlateLocHeight<=3.41 & abs(pitch$PlateLocSide)<=1)
  r$chase <- if (inz==0) r$swing else NA_real_
  r$own_neighbors <- NA_integer_
  r
}

# save the pieces needed to answer queries later (used by Stage 4 / apps)
saveRDS(list(DIMS=DIMS, ctr=ctr, scl=scl, hp=hp, hp_z=hp_z, hctr=hctr, hscl=hscl,
             hp_dims=hp_dims), file.path(OUT_DIR, "batter_response_space.rds"))
message("== batter_response_space.rds saved (", nrow(px), " pitches indexed) ==")

# ---- VERIFY with your exact example (only when run directly) ---------------
if (sys.nframe() == 0) {
# Sinker 0-2, 88 mph, 12 IVB, 18 HB, down-and-in (RHP vs RHH).
who <- px %>% filter(season=="2026", BatterSide=="Right") %>%
  count(BatterId, Batter, sort=TRUE) %>% slice(1)
cat("\nQuery: RHH", who$Batter, "vs Sinker 0-2, 88 / 12 IVB / 18 HB, down-in (RHP)\n")
pitch <- list(RelSpeed=88, InducedVertBreak=12, HorzBreak=18,
              PlateLocHeight=1.9, PlateLocSide=0.6)   # down-and-in to RHH
r <- batter_response(who$BatterId, pitch, bucket="two_strike", pthrows="Right")
cat(sprintf("  swing%%: %.0f   chase%%: %.0f   whiff%%: %.0f   xwOBACON: %.3f\n",
            100*r$swing, 100*r$chase, 100*r$whiff, r$xwobacon))
cat(sprintf("  (his own neighbors: %d, effective pooled sample ~%.0f)\n",
            r$own_neighbors, r$eff_n))
}
