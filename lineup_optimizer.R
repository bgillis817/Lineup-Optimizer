# ============================================================================
#  lineup_optimizer.R  - OFFENSE: find the batting order that maximizes our
#  expected runs vs a given pitcher. Fast heuristic seed + local perturbations,
#  each candidate simmed for xRV headline + run distribution.
#
#  Shares the same kernel/objects as the matchup sim:
#    source("continuous_common.R"); source("matchup_kernel_continuous.R")
#    source("lineup_sim_continuous.R"); source("lineup_optimizer.R")
# ============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(gtools) })

# PA weight by slot (slot 1 bats most). Rough MLB-style decay, fine for ranking.
SLOT_PA_WEIGHT <- c(4.65,4.55,4.44,4.34,4.23,4.13,4.02,3.92,3.81)

# per-hitter offensive value vs this pitcher (xRV per PA + on-base tendency)
hitter_values <- function(km, batters, pid, season) {
  map_dfr(batters, function(b){
    e <- if (SHAPE_MODE()) {
      bid <- if (exists("hp")) hp$BatterId[match(b, hp$Batter)] else b
      tryCatch(shape_pa_eval(bid, pid, season), error=function(err) NULL) %||%
        pa_eval(km, b, season, pid, tto = 1)
    } else pa_eval(km, b, season, pid, tto = 1)
    line <- e$line
    obp_like <- sum(line[c("BB","HBP","single","double","triple","home_run")])
    power    <- sum(line[c("double","triple","home_run")] * c(2,3,4))
    tibble(Batter=b, pa_xrv=e$xrv, obp=obp_like, power=power)
  })
}
`%||%` <- function(a,b) if (is.null(a)) b else a

# heuristic seed order: value in the high-PA slots, on-base ahead of power.
# Returns a batter vector in slot order 1..9.
heuristic_order <- function(hv) {
  hv <- hv %>% mutate(score = pa_xrv)               # base offensive value (xRV is +good on offense)
  # slots 1-2: high OBP; 3-5: best overall value/power; 6-9: descending value
  top    <- hv %>% arrange(desc(score))
  ob     <- hv %>% arrange(desc(obp))
  # build: 1=high OBP, 2=next OBP, 3=best value, 4=best power, 5=next value, rest by value
  chosen <- character(); pick <- function(df) { for (b in df$Batter) if (!(b %in% chosen)) return(b); NA }
  o <- character(9)
  o[1] <- pick(ob); chosen <- c(chosen,o[1])
  o[2] <- pick(ob); chosen <- c(chosen,o[2])
  o[3] <- pick(top); chosen <- c(chosen,o[3])
  o[4] <- pick(hv %>% arrange(desc(power))); chosen <- c(chosen,o[4])
  o[5] <- pick(top); chosen <- c(chosen,o[5])
  rest <- top %>% filter(!(Batter %in% chosen)) %>% pull(Batter)
  o[6:9] <- rest[1:4]
  o
}

# expected outing xRV for an order (offense), TTO-aware, ~innings*4.3 PAs
order_xrv <- function(km, order, pid, season, innings = 9) {
  n_pa <- round(innings * 4.3)
  tot <- 0
  for (i in seq_len(n_pa)) {
    slot <- ((i-1) %% 9) + 1; tto <- min(((i-1) %/% 9) + 1, 3)
    tot <- tot + pa_eval(km, order[slot], season, pid, tto)$xrv
  }
  tot
}

# full sim (offense): reuse the base-out MC via outing_lines-style rows
order_lines <- function(km, order, pid, season, innings = 9) {
  n_pa <- round(innings * 4.3)
  map_dfr(seq_len(n_pa), function(i){
    slot <- ((i-1) %% 9) + 1; tto <- min(((i-1) %/% 9) + 1, 3)
    e <- pa_eval(km, order[slot], season, pid, tto)
    tibble(i=i, slot=slot, Batter=order[slot], tto=tto, xrv=e$xrv, line=list(e$line))
  })
}
sim_order <- function(km, order, pid, season, innings = 9, n_sims = 2500) {
  rows <- order_lines(km, order, pid, season, innings)
  runs <- outing_mc(rows, n_sims)                     # from lineup_sim_continuous.R
  list(order=order, xrv=sum(rows$xrv), mean_runs=mean(runs),
       p10=quantile(runs,.10), p90=quantile(runs,.90),
       big=mean(runs>=5), runs=runs, rows=rows)
}

# cache pa_eval per (batter, tto) once â€” makes order scoring a pure lookup
# If SHAPE_MODE is on, use the shape-aware PA line (pitcher's real cloud x the
# batter's neighbor-pooled response) instead of the aggregate log5 kernel.
SHAPE_MODE <- function() isTRUE(getOption("necbl.shape_mode", FALSE))
build_pa_cache <- function(km, pool, pid, season) {
  cache <- new.env()
  if (SHAPE_MODE()) {
    for (b in pool) {
      bid <- if (exists("hp")) hp$BatterId[match(b, hp$Batter)] else b
      e <- tryCatch(shape_pa_eval(bid, pid, season), error=function(err) NULL)
      if (is.null(e)) e <- pa_eval(km, b, season, pid, 1)     # fallback
      for (tto in 1:3) assign(paste(b, tto, sep="\u0001"), e, envir = cache)
    }
    return(cache)
  }
  for (b in pool) for (tto in 1:3) {
    e <- pa_eval(km, b, season, pid, tto)
    assign(paste(b, tto, sep="\u0001"), e, envir = cache)
  }
  cache
}
.get <- function(cache, b, tto) get(paste(b, min(tto,3), sep="\u0001"), envir = cache)

# cheap deterministic proxy: summed pitcher xRV over the outing (order-sensitive)
order_xrv_c <- function(cache, order, innings = 9) {
  n_pa <- round(innings * 4.3); tot <- 0
  for (i in seq_len(n_pa)) { slot <- ((i-1)%%9)+1; tto <- ((i-1)%/%9)+1
    tot <- tot + .get(cache, order[slot], tto)$xrv }
  tot
}
order_lines_c <- function(cache, order, innings = 9) {
  n_pa <- round(innings * 4.3)
  map_dfr(seq_len(n_pa), function(i){ slot <- ((i-1)%%9)+1; tto <- ((i-1)%/%9)+1
    e <- .get(cache, order[slot], tto)
    tibble(i=i, slot=slot, Batter=order[slot], tto=min(tto,3), xrv=e$xrv, line=list(e$line)) })
}
sim_order_c <- function(cache, order, innings = 9, n_sims = 1200) {
  rows <- order_lines_c(cache, order, innings)
  runs <- outing_mc(rows, n_sims)
  list(order=order, xrv=sum(rows$xrv), mean_runs=mean(runs),
       p10=quantile(runs,.10), p90=quantile(runs,.90),
       big=mean(runs>=5), runs=runs, rows=rows)
}

# ---- deterministic EXPECTED RUNS for an order (no MC noise) ----------------
# Propagates the base-out state distribution through the order using the same
# advancement constants as the sim, returning exact expected runs. This is the
# objective the search maximizes, so search and final ranking agree.
.adv_expected <- function(b, outs, line) {
  # b: length-3 logical bases (1st,2nd,3rd). returns list(states=list(b,outs,prob,runs))
  ADV <- list(s2_score=0.5, s1_to3=0.3, d1_score=0.6)
  out <- list()
  add <- function(nb, no, p, r) if (p>0) out[[length(out)+1]] <<- list(b=nb, outs=no, prob=p, runs=r)
  pout <- unname(line["K"] + line["out"]); add(b, outs+1L, pout, 0)
  pbb  <- unname(line["BB"] + line["HBP"])
  if (pbb>0){ nb<-b; r<-0
    if(nb[1]){ if(nb[2]){ if(nb[3]) r<-r+1; nb[3]<-TRUE }; nb[2]<-TRUE }; nb[1]<-TRUE
    add(nb, outs, pbb, r) }
  p1 <- unname(line["single"])
  if (p1>0){ # branch on runner advancement (expected split)
    for (s2 in list(c(TRUE,ADV$s2_score),c(FALSE,1-ADV$s2_score)))
    for (s1 in list(c(TRUE,ADV$s1_to3),c(FALSE,1-ADV$s1_to3))) {
      nb<-c(FALSE,FALSE,FALSE); r<-0
      if(b[3]) r<-r+1
      if(b[2]){ if(s2[[1]]) r<-r+1 else nb[3]<-TRUE }
      if(b[1]){ if(s1[[1]]) nb[3]<-TRUE else nb[2]<-TRUE }
      nb[1]<-TRUE
      add(nb, outs, p1*s2[[2]]*s1[[2]], r) } }
  p2 <- unname(line["double"])
  if (p2>0){ for (d1 in list(c(TRUE,ADV$d1_score),c(FALSE,1-ADV$d1_score))) {
      nb<-c(FALSE,FALSE,FALSE); r<-0
      if(b[3]) r<-r+1; if(b[2]) r<-r+1
      if(b[1]){ if(d1[[1]]) r<-r+1 else nb[3]<-TRUE }
      nb[2]<-TRUE
      add(nb, outs, p2*d1[[2]], r) } }
  p3 <- unname(line["triple"]); if (p3>0){ add(c(FALSE,FALSE,TRUE), outs, p3, sum(b)) }
  phr<- unname(line["home_run"]); if (phr>0){ add(c(FALSE,FALSE,FALSE), outs, phr, sum(b)+1) }
  out
}
expected_runs_c <- function(cache, order, innings = 9) {
  memo_key <- paste0("ER\u0001", paste(order, collapse="|"), "\u0001", innings)
  hit <- tryCatch(get(memo_key, envir = cache), error = function(e) NULL)
  if (!is.null(hit)) return(hit)
  n_pa <- round(innings * 4.3); L <- 9
  # state key: paste(b1,b2,b3,outs). start empty, 0 outs, mass 1.
  st <- new.env(); st[["0|0|0|0"]] <- 1; ER <- 0
  keyf <- function(b,o) paste(as.integer(b[1]),as.integer(b[2]),as.integer(b[3]),o,sep="|")
  for (i in seq_len(n_pa)) {
    slot <- ((i-1)%%L)+1; tto <- ((i-1)%/%L)+1
    line <- .get(cache, order[slot], tto)$line
    nx <- new.env()
    for (k in ls(st)) {
      m <- st[[k]]; if (m<=0) next
      parts <- strsplit(k,"|",fixed=TRUE)[[1]]; b <- as.logical(as.integer(parts[1:3])); o <- as.integer(parts[4])
      for (tr in .adv_expected(b, o, line)) {
        ER <- ER + m * tr$prob * tr$runs
        if (tr$outs >= 3L) { kk <- "0|0|0|0" } else kk <- keyf(tr$b, tr$outs)
        nx[[kk]] <- (if (is.null(nx[[kk]])) 0 else nx[[kk]]) + m * tr$prob
      }
    }
    st <- nx
  }
  assign(memo_key, ER, envir = cache)
  ER
}

# ---- search engines (objective = deterministic expected runs) --------------
.oscore <- function(cache, o, innings) expected_runs_c(cache, o, innings)

# A) local search: swap until no pairwise swap improves (respects locks)
local_search <- function(cache, seed, innings, locked_slots) {
  o <- seed; free <- setdiff(1:9, locked_slots); best <- .oscore(cache,o,innings)
  improved <- TRUE
  while (improved) { improved <- FALSE
    for (a in free) for (b in free) if (a < b) {
      o2 <- o; o2[c(a,b)] <- o2[c(b,a)]; s2 <- .oscore(cache,o2,innings)
      if (s2 > best) { o <- o2; best <- s2; improved <- TRUE } } }
  o
}

# B) simulated annealing: random swaps, accept worse early, cool down
anneal <- function(cache, seed, innings, locked_slots, iters = 4000, t0 = 0.5) {
  free <- setdiff(1:9, locked_slots); if (length(free) < 2) return(seed)
  o <- seed; s <- .oscore(cache,o,innings); best_o <- o; best_s <- s
  for (k in seq_len(iters)) {
    temp <- t0 * (1 - k/iters) + 1e-4
    ns <- sample(free, 2); o2 <- o; o2[ns] <- o2[rev(ns)]
    s2 <- .oscore(cache,o2,innings); d <- s2 - s
    if (d > 0 || runif(1) < exp(d/temp)) { o <- o2; s <- s2 }
    if (s > best_s) { best_o <- o; best_s <- s } }
  best_o
}

# MAIN: best batting order from a POOL (>=9), MAXIMIZING expected runs.
#   method: "heuristic" (fast seed+swaps), "local" (hill-climb), "anneal"
optimize_lineup <- function(km, pool, pid, season, innings = 9,
                            n_top = 10, n_sims = 1200, method = "heuristic",
                            battle = NULL, battle_max = NULL, locked = NULL) {
  pool <- unique(pool); stopifnot(length(pool) >= 9)
  cache <- build_pa_cache(km, pool, pid, season)
  hv_all <- hitter_values(km, pool, pid, season)
  val <- setNames(hv_all$pa_xrv, hv_all$Batter)
  locked <- locked[!vapply(locked, function(x) is.null(x)||is.na(x)||x=="", logical(1))]
  locked_batters <- unlist(locked, use.names = FALSE)
  locked_slots   <- as.integer(names(locked))

  apply_lock <- function(o) {
    if (!length(locked)) return(o)
    for (s in names(locked)) { s <- as.integer(s); want <- locked[[as.character(s)]]
      cur <- which(o == want)
      if (length(cur) && cur != s) { tmp <- o[s]; o[s] <- want; o[cur] <- tmp } }
    o
  }
  valid_nine <- function() {
    combos <- utils::combn(pool, 9, simplify = FALSE)
    if (length(locked_batters))
      combos <- Filter(function(c9) all(locked_batters %in% c9), combos)
    if (!is.null(battle) && !is.null(battle_max))
      combos <- Filter(function(c9) sum(battle %in% c9) <= battle_max, combos)
    combos
  }
  nine_sets <- if (length(pool)==9) list(pool) else valid_nine()
  if (length(nine_sets) > 20) {
    ss <- vapply(nine_sets, function(s) sum(val[s]), numeric(1))
    nine_sets <- nine_sets[order(ss, decreasing = TRUE)][1:20]
  }

  # 2-opt hill climb: from a seed, repeatedly try all pairwise swaps of unlocked
  # slots, keep improvements (max proxy xRV), until no swap helps. Converges to a
  # real local optimum instead of a random neighbor.
  free_slots <- setdiff(1:9, locked_slots)
  climb <- function(seed) {
    o <- apply_lock(seed); best <- expected_runs_c(cache, o, innings); improved <- TRUE
    while (improved) {
      improved <- FALSE
      for (a in free_slots) for (b in free_slots) if (a < b) {
        o2 <- o; o2[c(a,b)] <- o2[c(b,a)]
        v <- expected_runs_c(cache, o2, innings)
        if (v > best + 1e-9) { o <- o2; best <- v; improved <- TRUE }
      }
    }
    o
  }

  # per 9-set: seeds + search depending on method
  all_orders <- list()
  for (s9 in nine_sets) {
    hv <- hv_all %>% filter(Batter %in% s9)
    base_seeds <- list(
      heuristic_order(hv),
      hv %>% arrange(desc(pa_xrv)) %>% pull(Batter),
      hv %>% arrange(desc(power))  %>% pull(Batter))
    if (method == "heuristic") {
      # fast: just the heuristic seed, lightly improved
      all_orders[[length(all_orders)+1]] <- apply_lock(base_seeds[[1]])
    } else if (method == "local") {
      seeds <- c(base_seeds, list(sample(s9), sample(s9)))
      for (sd in seeds) all_orders[[length(all_orders)+1]] <- climb(sd)
    } else { # anneal: climb all sets, but anneal only the strongest few
      seeds <- c(base_seeds, list(sample(s9), sample(s9)))
      for (sd in seeds) all_orders[[length(all_orders)+1]] <- climb(sd)
      if (which(vapply(nine_sets, identical, logical(1), s9))[1] <= 5)
        all_orders[[length(all_orders)+1]] <- anneal(cache, apply_lock(base_seeds[[1]]), innings, locked_slots, iters=1500)
    }
  }
  key <- vapply(all_orders, paste, character(1), collapse="|"); all_orders <- all_orders[!duplicated(key)]

  # expand: add single-swap neighbors of the best orders so we can show a real
  # top-N of DISTINCT lineups (searches collapse to one optimum otherwise).
  free_all <- setdiff(1:9, locked_slots)
  seed_best <- all_orders[order(map_dbl(all_orders, ~ expected_runs_c(cache, .x, innings)),
                                decreasing = TRUE)]
  for (o in head(seed_best, 3)) {
    for (a in free_all) for (b in free_all) if (a < b) {
      o2 <- o; o2[c(a,b)] <- o2[c(b,a)]; all_orders[[length(all_orders)+1]] <- o2 }
  }
  key <- vapply(all_orders, paste, character(1), collapse="|"); all_orders <- all_orders[!duplicated(key)]

  # rank distinct orders by expected runs, sim the top n_top, final rank by MEAN RUNS
  scored <- tibble(idx = seq_along(all_orders),
                   px = map_dbl(all_orders, ~ expected_runs_c(cache, .x, innings))) %>%
    arrange(desc(px)) %>% head(n_top)
  results <- map(scored$idx, function(ix) sim_order_c(cache, all_orders[[ix]], innings, n_sims))
  results[order(map_dbl(results, ~ .x$mean_runs), decreasing = TRUE)]
}
