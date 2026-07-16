source("continuous_common.R")
px <- attach_xwobacon(load_pitches())
keep <- c("BatterId","Batter","BatterSide","BatterTeam","bucket","InZone","IsBIP",
          "IsContact","IsFoul","IsSwing","IsWhiff","PitcherId","Pitcher","PitcherThrows",
          "season","xp_1b","xp_2b","xp_3b","xp_hr","xp_out","xwobacon",
          "RelSpeed","InducedVertBreak","HorzBreak","PlateLocHeight","PlateLocSide",
          "GameID","PAofInning","Inning","Top/Bottom","ExitSpeed","Angle")
miss <- setdiff(keep, names(px)); if (length(miss)) cat("MISSING:", paste(miss, collapse=", "), "\n")
px <- px[, intersect(keep, names(px)), drop=FALSE]
cat("cols:", ncol(px), " rows:", nrow(px), "\n")
saveRDS(px, "data/px_cache.rds", compress="gzip")
cat("MB:", round(file.size("data/px_cache.rds")/1e6,1), "\n")
