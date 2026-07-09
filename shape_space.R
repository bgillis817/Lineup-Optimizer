# ============================================================================
#  shape_space.R  - STAGE 1 of the shape-aware matchup model.
#  Builds a pitch "shape space" from velo + IVB + HB + RelHeight + RelSide,
#  kept SEPARATE by pitcher handedness. Provides a nearest-neighbor lookup so
#  you can verify the geometry ("do similar pitches actually come back similar?")
#  BEFORE any batter/outcome math is layered on top.
#
#  Similarity = plain scaled Euclidean distance in shape space (per handedness).
#  Scaling puts velo(mph)/break(in)/release(ft) on comparable footing so no
#  single axis dominates the distance.
#
#  RUN (standalone check):
#    DRIVE_FOLDER="Navs CSVs" GDRIVE_KEY_PATH=sa_key.json Rscript shape_space.R
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse) })

SHAPE_DIMS <- c("RelSpeed","InducedVertBreak","HorzBreak","RelHeight","RelSide")

# build the shape space: cleaned pitch table -> per-hand scaled matrices + meta
build_shape_space <- function(px) {
  d <- px %>%
    mutate(across(all_of(SHAPE_DIMS), ~ suppressWarnings(as.numeric(.)))) %>%
    filter(if_all(all_of(SHAPE_DIMS), ~ !is.na(.)),
           PitcherThrows %in% c("Left","Right"))
  spaces <- list()
  for (hand in c("Left","Right")) {
    dh <- d %>% filter(PitcherThrows == hand)
    if (!nrow(dh)) next
    M <- as.matrix(dh[, SHAPE_DIMS])
    ctr <- colMeans(M); scl <- apply(M, 2, sd); scl[scl == 0] <- 1
    Z <- sweep(sweep(M, 2, ctr, "-"), 2, scl, "/")
    spaces[[hand]] <- list(
      Z = Z, center = ctr, scale = scl,
      meta = bind_cols(
        dh %>% select(any_of(c("Pitcher","PitcherId","PitcherTeam","AutoPitchType","season"))),
        as_tibble(M)))
  }
  spaces
}

# nearest neighbors of a single pitch (given as a named vector of SHAPE_DIMS)
nn_pitch <- function(spaces, hand, pitch_vec, k = 15) {
  sp <- spaces[[hand]]; if (is.null(sp)) stop("no pitches for hand ", hand)
  q <- (unlist(pitch_vec[SHAPE_DIMS]) - sp$center) / sp$scale
  d2 <- rowSums(sweep(sp$Z, 2, q, "-")^2)
  ord <- order(d2)[seq_len(min(k, nrow(sp$Z)))]
  sp$meta[ord, ] %>% mutate(dist = sqrt(d2[ord]))
}

# a pitcher's average shape (centroid per pitch type) -> useful for eyeballing
pitcher_arsenal <- function(spaces, pitcher_name) {
  bind_rows(lapply(names(spaces), function(h){
    m <- spaces[[h]]$meta
    m %>% filter(Pitcher == pitcher_name) %>%
      group_by(AutoPitchType) %>%
      summarise(n = n(), across(all_of(SHAPE_DIMS), ~ round(mean(.),1)), .groups="drop") %>%
      mutate(Hand = h)
  }))
}

# ---- standalone verification run -------------------------------------------
if (sys.nframe() == 0) {
  source("continuous_common.R")
  px <- load_pitches()
  spaces <- build_shape_space(px)
  saveRDS(spaces, "data/shape_space.rds")
  cat("== shape_space.rds built ==\n")
  for (h in names(spaces)) cat(sprintf("  %s-handed pitches: %d\n", h, nrow(spaces[[h]]$Z)))

  # DEMO: take a random slider (if present) and show its nearest neighbors,
  # so you can confirm similar pitches come back similar.
  demo_hand <- names(spaces)[1]; meta <- spaces[[demo_hand]]$meta
  sl <- which(meta$AutoPitchType %in% c("Slider","slider"))
  if (length(sl)) {
    i <- sl[1]; pv <- meta[i, SHAPE_DIMS]
    cat(sprintf("\nDemo: neighbors of a %s-handed %s (%.1f mph, IVB %.1f, HB %.1f):\n",
                demo_hand, meta$AutoPitchType[i], pv$RelSpeed, pv$InducedVertBreak, pv$HorzBreak))
    print(nn_pitch(spaces, demo_hand, pv, k=12) %>%
            select(Pitcher, AutoPitchType, all_of(SHAPE_DIMS), dist))
  }
}
