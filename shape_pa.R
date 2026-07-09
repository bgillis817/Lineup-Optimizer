# ============================================================================
#  shape_pa.R  - BRIDGE (Stage 4 -> the sims).
#  Converts the shape-aware pitch-level matchup (pitcher's real cloud x batter's
#  neighbor-pooled response) into a PA-level line: P(K, BB, HBP, out, 1B, 2B,
#  3B, HR) + xRV. Same shape as pa_eval() so the existing MC sims consume it
#  with no changes.
#
#  Per count (b,s): sample the pitcher's actual pitches for that bucket, ask the
#  batter's response at each, average -> per-pitch primitives:
#     no swing  -> ball (if out of zone) or called strike (if in zone)
#     swing     -> whiff (strike) | contact -> foul | in play -> BIP outcome
#  Then walk the count Markov to absorption.
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })

# map a (balls,strikes) count to the leverage bucket used everywhere else
count_bucket <- function(b, s) {
  if (b==0 && s==0) return("0-0")
  if (b==3 && s==2) return("3-2")
  if (s==2)         return("two_strike")
  if (b >  s)       return("behind")     # pitcher behind
  if (s >  b)       return("ahead")      # pitcher ahead
  "even"
}

HBP_RATE <- 0.010   # league-ish per-PA hit-by-pitch, applied as a small constant
.SHAPE_CACHE <- new.env(parent = emptyenv())

# per-pitch primitives for one bucket: ONE batched query over the pitcher's
# actual cloud (his cloud carries his stuff + command)
bucket_primitives <- function(batter_id, pid, bucket, season=NULL, n_pitch=20) {
  cl <- pitcher_cloud(pid, bucket, n_pitch, season)
  if (is.null(cl) || !nrow(cl)) return(NULL)
  ph <- cl$PitcherThrows[1]
  r <- batter_response_batch(batter_id, cl, bucket, ph)
  if (is.null(r) || is.na(r$swing)) return(NULL)
  # zone rate comes from the pitcher's OWN pitches (where he actually threw it)
  p_inz <- mean(cl$PlateLocHeight>=1.59 & cl$PlateLocHeight<=3.41 & abs(cl$PlateLocSide)<=1)
  p_swing <- r$swing
  p_whiff <- ifelse(is.na(r$whiff), 0.25, r$whiff)
  p_foul  <- ifelse(is.na(r$foul),  0.50, r$foul)
  bip <- r$bip_dist; bip <- bip/sum(bip)
  p_take  <- 1 - p_swing
  p_cont  <- p_swing * (1 - p_whiff)
  c(ball    = p_take * (1 - p_inz),
    cstrike = p_take * p_inz,
    sstrike = p_swing * p_whiff,
    foul    = p_cont * p_foul,
    bip     = p_cont * (1 - p_foul),
    out=bip[1], single=bip[2], double=bip[3], triple=bip[4], home_run=bip[5])
}

# full PA: walk the count Markov using bucket-specific primitives (memoized)
shape_pa_eval <- function(batter_id, pid, season=NULL, n_pitch=15, rv=NULL) {
  key <- paste(batter_id, pid, season, n_pitch, sep="\u0001")
  hit <- .SHAPE_CACHE[[key]]; if (!is.null(hit)) return(hit)
  prim <- list()
  for (bk in c("0-0","ahead","even","behind","two_strike","3-2")) {
    p <- bucket_primitives(batter_id, pid, bk, season, n_pitch)
    if (is.null(p)) p <- bucket_primitives(batter_id, pid, "0-0", season, n_pitch)
    prim[[bk]] <- p
  }
  if (all(vapply(prim, is.null, logical(1)))) return(NULL)
  fill <- prim[[which(!vapply(prim, is.null, logical(1)))[1]]]
  for (bk in names(prim)) if (is.null(prim[[bk]])) prim[[bk]] <- fill

  # absorbing outcome accumulators
  outc <- c(K=0, BB=0, out=0, single=0, double=0, triple=0, home_run=0)
  st <- matrix(0, 4, 3); st[1,1] <- 1        # [balls+1, strikes+1]
  for (iter in 1:30) {                        # plenty for a PA to absorb
    nxt <- matrix(0, 4, 3); moved <- FALSE
    for (b in 0:3) for (s in 0:2) {
      m <- st[b+1, s+1]; if (m <= 0) next
      moved <- TRUE
      p <- prim[[count_bucket(b, s)]]
      # ball
      if (b+1 >= 4) outc["BB"] <- outc["BB"] + m*p["ball"] else nxt[b+2, s+1] <- nxt[b+2, s+1] + m*p["ball"]
      # strikes (called + swinging)
      pk <- p["cstrike"] + p["sstrike"]
      if (s+1 >= 3) outc["K"] <- outc["K"] + m*pk else nxt[b+1, s+2] <- nxt[b+1, s+2] + m*pk
      # foul: adds a strike unless already 2
      if (s < 2) nxt[b+1, s+2] <- nxt[b+1, s+2] + m*p["foul"] else nxt[b+1, s+1] <- nxt[b+1, s+1] + m*p["foul"]
      # ball in play -> distribute over outcomes
      bipm <- m * p["bip"]
      for (o in c("out","single","double","triple","home_run")) outc[o] <- outc[o] + bipm*p[o]
    }
    st <- nxt
    if (!moved || sum(st) < 1e-9) break
  }
  # any residual mass (long fouls) -> proportional to current absorbing mix
  resid <- max(0, 1 - sum(outc))
  if (resid > 1e-9 && sum(outc) > 0) outc <- outc + resid * outc/sum(outc)

  line <- c(K=unname(outc["K"]), BB=unname(outc["BB"])*(1-HBP_RATE), HBP=HBP_RATE,
            out=unname(outc["out"]), single=unname(outc["single"]),
            double=unname(outc["double"]), triple=unname(outc["triple"]),
            home_run=unname(outc["home_run"]))
  line <- line/sum(line)
  RV <- if (is.null(rv)) c(K=-0.089, BB=0.056, HBP=0.31, out=-0.26, single=0.44,
                           double=0.75, triple=1.01, home_run=1.40) else rv
  res <- list(line=line, xrv=sum(line * RV[names(line)]))
  assign(key, res, envir = .SHAPE_CACHE)
  res
}
