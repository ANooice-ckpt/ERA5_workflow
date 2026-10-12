#!/usr/bin/env Rscript
# Small readable local summaries after the pack validator has passed.
args <- commandArgs(TRUE)
if (length(args)!=2L) stop("Usage: Rscript 08_summarize_pack.R config.R audit_dir")
cmd <- commandArgs(FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
source(file.path(dirname(normalizePath(script)),"00_core_loader.R"))
cfg <- stream_read_config(args[1L]); out <- args[2L]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
index <- readRDS(file.path(cfg$local$pack_dir,"pack_index.rds"))
annual <- do.call(rbind,lapply(index$year_files,function(y) data.frame(
  year=y$year,selected_cells=y$n_cells,days=y$n_days,bytes=y$bytes,
  GiB=y$bytes/1024^3,md5=y$md5)))
write.csv(annual,file.path(out,"weather_pack_by_year.csv"),row.names=FALSE)
write.csv(index$history_audit,file.path(out,"source_history_audit.csv"),row.names=FALSE)
write.csv(index$temperature_audit,file.path(out,"source_temperature_audit.csv"),row.names=FALSE)
status <- list(status="LOCAL_PREPARATION_COMPLETE",created=Sys.time(),
  pack_identity=index$identity,selected_cells=index$n_cells,tiles=index$n_tiles,
  tile_side=index$tile_side,days=sum(annual$days),years=annual$year,
  packed_bytes=sum(annual$bytes),uncompressed_weather_bytes=index$n_cells*sum(annual$days)*5*8,
  scientific_domain="0.25 degrees; -60 to +80 latitude; any-land overlap; 2015 population density >1 per full cell km2",
  server_prediction_status="NOT_RUN",cdd_scenarios_declared=length(cfg$cdd$scenarios)>0)
saveRDS(status,file.path(out,"local_preparation_status.rds"))
text <- c("Local preparation complete; no fitted-model projection has run.",
  paste("Selected grids:",index$n_cells),paste("Spatial tiles:",index$n_tiles),
  paste("Years:",paste(range(annual$year),collapse="-")),paste("UTC days:",sum(annual$days)),
  sprintf("Packed NetCDF total: %.3f GiB (%.3f GB)",sum(annual$bytes)/1024^3,sum(annual$bytes)/1e9),
  paste("Pack identity:",index$identity),
  "All annual packs passed MD5, coordinates/cell IDs, dimensions and source sample comparisons.",
  "CDD theta0/slope/kappa remain to be declared before server prediction.")
writeLines(text,file.path(out,"LOCAL_PREPARATION_STATUS.txt"))
cat(paste(text,collapse="\n"),"\n")
