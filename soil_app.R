
# app.R
# Soil Texture Comparison Explorer
#
# This Shiny app:
#   * accepts a CSV with columns site_id, lat, lon (or long)
#   * de-duplicates exact repeated rows
#   * runs three soil backends for every unique coordinate:
#       1) OpenLandMap / LandGIS v0.2 topsoil rasters (250 m, b0..0cm)
#       2) SoilGrids 2.0 topsoil rasters (250 m, 0-5 cm mean; g/kg -> %)
#       3) HWSD v2.0 top layer (30 arc-second raster + SQLite attribute DB; D1 0-20 cm)
#   * shows a map and selected-site comparison
#   * previews wide and long results
#   * exports wide and long CSV outputs
#
# Important HWSD note:
#   This version no longer depends on reading the official HWSD .mdb file through
#   ODBC or mdbtools. Instead it uses the official HWSD raster plus a SQLite copy
#   of the HWSD v2 attribute database, matching the documented cross-platform R
#   workflow and QBMS helper workflow.

required_pkgs <- c(
  "shiny", "bslib", "leaflet", "terra", "readr", "dplyr",
  "glue", "htmltools", "DT", "tibble", "tidyr", "DBI", "RSQLite"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install the missing packages before running this app:\n  ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(leaflet)
  library(terra)
  library(readr)
  library(dplyr)
  library(glue)
  library(htmltools)
  library(DT)
  library(tibble)
  library(tidyr)
  library(DBI)
  library(RSQLite)
})

options(shiny.maxRequestSize = 100 * 1024^2)
options(timeout = max(7200, getOption("timeout", 60)))
terra::terraOptions(progress = 0)
Sys.setenv(
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".vrt,.ovr,.tif,.tiff,.xml",
  GDAL_HTTP_MAX_RETRY = "3",
  GDAL_HTTP_RETRY_DELAY = "1"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (is.character(x) && length(x) == 1 && !nzchar(x)) return(y)
  x
}

# ------------------------------------------------------------------------------
# Source specifications
# ------------------------------------------------------------------------------

OPENLANDMAP_SPECS <- tibble::tribble(
  ~property,    ~filename,                                                          ~url,                                                                                                                                          ~md5,                                ~doi,
  "sand_pct",   "sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif", "https://zenodo.org/records/2525662/files/sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif?download=1", "2e00065c107b4ccb3064ebde255c18b2", "10.5281/zenodo.2525662",
  "silt_pct",   "sol_silt.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif", "https://zenodo.org/records/2525676/files/sol_silt.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif?download=1", "f69a36046485d87d0b9f05f18c58422c", "10.5281/zenodo.2525676",
  "clay_pct",   "sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif", "https://zenodo.org/records/2525663/files/sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif?download=1", "5c6ab29f9068a9fae746b5aa02d2d535", "10.5281/zenodo.2525663"
)

SOILGRIDS_VRTS <- c(
  sand_pct = "https://files.isric.org/soilgrids/latest/data/sand/sand_0-5cm_mean.vrt",
  silt_pct = "https://files.isric.org/soilgrids/latest/data/silt/silt_0-5cm_mean.vrt",
  clay_pct = "https://files.isric.org/soilgrids/latest/data/clay/clay_0-5cm_mean.vrt"
)

HWSD_URLS <- list(
  raster_zip = "https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/HWSD/HWSD2_RASTER.zip",
  sqlite = "https://www.isric.org/sites/default/files/HWSD2.sqlite"
)

BACKEND_ORDER <- c("openlandmap", "soilgrids", "hwsd_v2")
BACKEND_LABELS <- c(
  openlandmap = "OpenLandMap",
  soilgrids = "SoilGrids",
  hwsd_v2 = "HWSD v2.0"
)

# ------------------------------------------------------------------------------
# Cache directories
# ------------------------------------------------------------------------------

default_cache_dir <- tryCatch(
  {
    if ("R_user_dir" %in% getNamespaceExports("tools")) {
      tools::R_user_dir("soil_texture_comparison_explorer", which = "cache")
    } else {
      file.path(tempdir(), "soil_texture_comparison_cache")
    }
  },
  error = function(e) file.path(tempdir(), "soil_texture_comparison_cache")
)

CACHE_DIR <- Sys.getenv("SOIL_CACHE_DIR", unset = default_cache_dir)
CACHE_OPENLANDMAP <- file.path(CACHE_DIR, "openlandmap")
CACHE_HWSD <- file.path(CACHE_DIR, "hwsd")

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_OPENLANDMAP, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_HWSD, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

format_pct <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  format(round(x, 1), trim = TRUE, nsmall = 1, scientific = FALSE)
}

format_num <- function(x, digits = 7) {
  if (length(x) == 0 || is.na(x)) return("NA")
  format(x, digits = digits, trim = TRUE, scientific = FALSE)
}

safe_read_lines <- function(path) {
  if (!file.exists(path)) return(character())
  tryCatch(readLines(path, warn = FALSE), error = function(e) character())
}

md5_sidecar <- function(file_path) {
  paste0(file_path, ".md5")
}

safe_md5 <- function(file_path) {
  tryCatch(unname(tools::md5sum(file_path)), error = function(e) NA_character_)
}

is_md5_valid <- function(file_path, expected_md5) {
  if (!file.exists(file_path)) return(FALSE)

  sidecar_path <- md5_sidecar(file_path)
  if (file.exists(sidecar_path)) {
    cached_md5 <- trimws(paste(safe_read_lines(sidecar_path), collapse = ""))
    if (nzchar(cached_md5) && identical(tolower(cached_md5), tolower(expected_md5))) {
      return(TRUE)
    }
  }

  actual_md5 <- safe_md5(file_path)
  if (is.na(actual_md5)) return(FALSE)

  ok <- identical(tolower(actual_md5), tolower(expected_md5))
  if (ok) writeLines(expected_md5, sidecar_path)
  ok
}

download_file_robust <- function(url, destfile, mode = "wb") {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  tmpfile <- paste0(destfile, ".part")
  if (file.exists(tmpfile)) unlink(tmpfile, force = TRUE)

  ok <- FALSE
  err <- NULL

  if (requireNamespace("curl", quietly = TRUE)) {
    tryCatch({
      curl::curl_download(url = url, destfile = tmpfile, quiet = TRUE, mode = mode)
      ok <- TRUE
    }, error = function(e) {
      err <<- e$message
    })
  }

  if (!ok) {
    tryCatch({
      utils::download.file(
        url = url,
        destfile = tmpfile,
        mode = mode,
        method = "libcurl",
        quiet = TRUE
      )
      ok <- TRUE
    }, error = function(e) {
      err <<- e$message
    })
  }

  if (!ok || !file.exists(tmpfile) || is.na(file.info(tmpfile)$size) || file.info(tmpfile)$size <= 0) {
    if (file.exists(tmpfile)) unlink(tmpfile, force = TRUE)
    stop("Download failed for ", basename(destfile), if (!is.null(err)) paste0(": ", err) else ".", call. = FALSE)
  }

  if (file.exists(destfile)) unlink(destfile, force = TRUE)
  moved <- file.rename(tmpfile, destfile)
  if (!moved) {
    copied <- file.copy(tmpfile, destfile, overwrite = TRUE)
    unlink(tmpfile, force = TRUE)
    if (!copied) stop("Downloaded file could not be moved into cache: ", basename(destfile), call. = FALSE)
  }

  destfile
}

safe_unzip <- function(zipfile, exdir) {
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zipfile, exdir = exdir, overwrite = TRUE)
  invisible(exdir)
}

find_first_file <- function(root_dir, pattern, ignore.case = TRUE) {
  if (!dir.exists(root_dir)) return(NA_character_)
  hits <- list.files(
    root_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = ignore.case
  )
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

safe_crs <- function(x) {
  tryCatch(terra::crs(x), error = function(e) "") %||% ""
}

has_crs <- function(x) {
  nzchar(safe_crs(x))
}

stack_rasters <- function(rasters) {
  stopifnot(length(rasters) >= 1)
  out <- rasters[[1]]
  if (length(rasters) > 1) {
    for (i in 2:length(rasters)) {
      out <- c(out, rasters[[i]])
    }
  }
  out
}

weighted_mean_safe <- function(x, w) {
  idx <- is.finite(x) & is.finite(w) & w > 0
  if (!any(idx)) return(NA_real_)
  sum(x[idx] * w[idx]) / sum(w[idx])
}

normalize_texture_triplet <- function(sand, silt, clay) {
  vals <- c(sand, silt, clay)
  if (any(!is.finite(vals))) return(c(NA_real_, NA_real_, NA_real_))
  if (any(vals < 0)) return(c(NA_real_, NA_real_, NA_real_))
  total <- sum(vals)
  if (!is.finite(total) || total <= 0) return(c(NA_real_, NA_real_, NA_real_))
  vals / total * 100
}

classify_usda_texture <- function(sand, silt, clay) {
  vals <- normalize_texture_triplet(sand, silt, clay)
  sand <- vals[[1]]
  silt <- vals[[2]]
  clay <- vals[[3]]

  if (any(is.na(vals))) return(NA_character_)

  if ((silt + 1.5 * clay) < 15) return("Sand")
  if ((silt + 1.5 * clay) >= 15 && (silt + 2 * clay) < 30) return("Loamy sand")
  if (
    (clay >= 7 && clay < 20 && sand > 52 && (silt + 2 * clay) >= 30) ||
      (clay < 7 && silt < 50 && (silt + 2 * clay) >= 30)
  ) return("Sandy loam")
  if (clay >= 7 && clay < 27 && silt >= 28 && silt < 50 && sand <= 52) return("Loam")
  if ((silt >= 50 && clay >= 12 && clay < 27) || (silt >= 50 && silt < 80 && clay < 12)) return("Silt loam")
  if (silt >= 80 && clay < 12) return("Silt")
  if (clay >= 20 && clay < 35 && silt < 28 && sand > 45) return("Sandy clay loam")
  if (clay >= 27 && clay < 40 && sand > 20 && sand <= 45) return("Clay loam")
  if (clay >= 27 && clay < 40 && sand <= 20) return("Silty clay loam")
  if (clay >= 35 && sand > 45) return("Sandy clay")
  if (clay >= 40 && silt >= 40) return("Silty clay")
  if (clay >= 40 && sand <= 45 && silt < 40) return("Clay")

  NA_character_
}

format_fraction_sum_flag <- function(sand, silt, clay, tol = 0.75) {
  vals <- c(sand, silt, clay)
  if (any(!is.finite(vals))) return(NA_character_)
  total <- sum(vals)
  if (abs(total - 100) <= tol) return(NA_character_)
  "fraction_sum_not_100"
}

coalesce_flag <- function(...) {
  vals <- unlist(list(...))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(NA_character_)
  paste(unique(vals), collapse = "; ")
}

safe_rbind <- function(x) {
  if (length(x) == 0) return(tibble())
  dplyr::bind_rows(x)
}

# ------------------------------------------------------------------------------
# Input preparation
# ------------------------------------------------------------------------------

prepare_sites <- function(file_path) {
  df <- readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE, trim_ws = TRUE)
  original_n <- nrow(df)

  names(df) <- tolower(trimws(names(df)))
  if ("long" %in% names(df) && !"lon" %in% names(df)) {
    names(df)[names(df) == "long"] <- "lon"
  }

  required_cols <- c("site_id", "lat", "lon")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Uploaded CSV must contain columns site_id, lat, lon. Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- df %>%
    transmute(
      site_id = trimws(as.character(site_id)),
      lat = suppressWarnings(as.numeric(lat)),
      lon = suppressWarnings(as.numeric(lon))
    )

  dropped_missing <- sum(is.na(out$lat) | is.na(out$lon) | is.na(out$site_id) | out$site_id == "")
  out <- out %>% filter(!(is.na(lat) | is.na(lon) | is.na(site_id) | site_id == ""))

  dropped_range <- sum(out$lat < -90 | out$lat > 90 | out$lon < -180 | out$lon > 180)
  out <- out %>% filter(lat >= -90, lat <= 90, lon >= -180, lon <= 180)

  n_before_distinct <- nrow(out)
  out <- out %>% distinct(site_id, lat, lon)
  dropped_duplicate_rows <- n_before_distinct - nrow(out)

  if (nrow(out) == 0) {
    stop("No valid rows remain after cleaning the uploaded file.", call. = FALSE)
  }

  site_counts <- out %>% count(site_id, name = "site_id_n")
  out <- out %>%
    left_join(site_counts, by = "site_id") %>%
    mutate(
      site_label = ifelse(
        site_id_n > 1,
        glue("{site_id} [{format(lat, digits = 7)}, {format(lon, digits = 7)}]"),
        site_id
      )
    ) %>%
    select(site_id, lat, lon, site_label)

  coords <- out %>%
    distinct(lat, lon) %>%
    arrange(lat, lon) %>%
    mutate(coord_id = row_number()) %>%
    select(coord_id, lat, lon)

  out <- out %>% left_join(coords, by = c("lat", "lon"))

  notes <- c(
    glue("Rows read: {original_n}"),
    glue("Valid unique site rows retained: {nrow(out)}"),
    glue("Unique coordinates to query: {nrow(coords)}")
  )
  if (dropped_missing > 0) notes <- c(notes, glue("Dropped rows with missing or non-numeric site_id/lat/lon: {dropped_missing}"))
  if (dropped_range > 0) notes <- c(notes, glue("Dropped rows outside valid coordinate ranges: {dropped_range}"))
  if (dropped_duplicate_rows > 0) notes <- c(notes, glue("Dropped exact duplicate rows (same site_id, lat, lon): {dropped_duplicate_rows}"))
  if (any(site_counts$site_id_n > 1)) {
    notes <- c(notes, "Some site_id values were repeated with different coordinates; coordinates are shown in the selector to disambiguate them.")
  }

  list(
    sites = out,
    coords = coords,
    notes = notes
  )
}

# ------------------------------------------------------------------------------
# OpenLandMap backend
# ------------------------------------------------------------------------------

ensure_openlandmap_rasters <- function(progress_callback = NULL) {
  paths <- character(nrow(OPENLANDMAP_SPECS))

  for (i in seq_len(nrow(OPENLANDMAP_SPECS))) {
    spec_row <- OPENLANDMAP_SPECS[i, , drop = FALSE]
    if (is.function(progress_callback)) {
      progress_callback(
        label = glue("Preparing OpenLandMap: {spec_row$property[[1]]}"),
        detail = spec_row$filename[[1]]
      )
    }

    dest <- file.path(CACHE_OPENLANDMAP, spec_row$filename[[1]])
    if (file.exists(dest) && is_md5_valid(dest, spec_row$md5[[1]])) {
      paths[i] <- dest
      next
    }

    download_file_robust(spec_row$url[[1]], dest)
    if (!is_md5_valid(dest, spec_row$md5[[1]])) {
      unlink(dest, force = TRUE)
      unlink(md5_sidecar(dest), force = TRUE)
      stop("MD5 verification failed for downloaded file: ", basename(dest), call. = FALSE)
    }

    paths[i] <- dest
  }

  names(paths) <- OPENLANDMAP_SPECS$property
  paths
}

load_openlandmap_stack <- function(progress_callback = NULL) {
  raster_paths <- ensure_openlandmap_rasters(progress_callback = progress_callback)
  rasters <- lapply(raster_paths, terra::rast)

  ref <- rasters[[1]]
  same_geom <- all(vapply(
    rasters[-1],
    function(r) terra::compareGeom(ref, r, stopOnError = FALSE),
    logical(1)
  ))
  if (!same_geom) {
    stop("OpenLandMap sand, silt, and clay rasters do not share identical geometry.", call. = FALSE)
  }

  soil_stack <- stack_rasters(rasters)
  names(soil_stack) <- names(raster_paths)
  soil_stack
}

extract_openlandmap <- function(coords_df, progress_callback = NULL) {
  if (is.function(progress_callback)) {
    progress_callback(label = "OpenLandMap", detail = "Loading cached rasters")
  }

  soil_stack <- load_openlandmap_stack(progress_callback = progress_callback)
  pts <- terra::vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  vals <- terra::extract(soil_stack, pts, ID = FALSE) %>% tibble::as_tibble()

  out <- bind_cols(coords_df, vals) %>%
    mutate(
      source_id = "openlandmap",
      source_name = "OpenLandMap",
      native_depth = "b0..0cm",
      native_resolution_m = 250,
      extraction_method = "Point-in-raster-cell extraction from cached GeoTIFFs",
      data_source = "OpenLandMap / LandGIS v0.2 topsoil rasters (250 m, b0..0cm) from Zenodo DOI 10.5281/zenodo.2525662, 10.5281/zenodo.2525676, 10.5281/zenodo.2525663.",
      sand_pct = as.numeric(sand_pct),
      silt_pct = as.numeric(silt_pct),
      clay_pct = as.numeric(clay_pct),
      qa_flag = mapply(format_fraction_sum_flag, sand_pct, silt_pct, clay_pct),
      texture_class_usda = mapply(classify_usda_texture, sand_pct, silt_pct, clay_pct),
      error_message = NA_character_
    )

  out
}

# ------------------------------------------------------------------------------
# SoilGrids backend
# ------------------------------------------------------------------------------

open_soilgrids_layer <- function(url) {
  attempts <- list(
    function() terra::rast(url, vsi = TRUE),
    function() terra::rast(paste0("/vsicurl/", url))
  )

  errors <- character()
  for (attempt in attempts) {
    result <- tryCatch(attempt(), error = function(e) e)
    if (!inherits(result, "error")) {
      if (!has_crs(result)) {
        try(terra::crs(result) <- "ESRI:54052", silent = TRUE)
      }
      return(result)
    }
    errors <- c(errors, conditionMessage(result))
  }

  stop(
    "Could not open SoilGrids layer: ", url, "\n",
    paste(unique(errors), collapse = "\n"),
    call. = FALSE
  )
}

load_soilgrids_stack <- function(progress_callback = NULL) {
  rasters <- vector("list", length(SOILGRIDS_VRTS))
  for (i in seq_along(SOILGRIDS_VRTS)) {
    nm <- names(SOILGRIDS_VRTS)[[i]]
    if (is.function(progress_callback)) {
      progress_callback(label = "SoilGrids", detail = basename(SOILGRIDS_VRTS[[i]]))
    }
    rasters[[i]] <- open_soilgrids_layer(SOILGRIDS_VRTS[[i]])
  }

  ref <- rasters[[1]]
  same_geom <- all(vapply(
    rasters[-1],
    function(r) terra::compareGeom(ref, r, stopOnError = FALSE),
    logical(1)
  ))
  if (!same_geom) {
    stop("SoilGrids sand, silt, and clay rasters do not share identical geometry.", call. = FALSE)
  }

  soil_stack <- stack_rasters(rasters)
  names(soil_stack) <- names(SOILGRIDS_VRTS)
  if (!has_crs(soil_stack)) {
    try(terra::crs(soil_stack) <- "ESRI:54052", silent = TRUE)
  }
  soil_stack
}

extract_soilgrids <- function(coords_df, progress_callback = NULL) {
  soil_stack <- load_soilgrids_stack(progress_callback = progress_callback)
  pts <- terra::vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  pts_proj <- if (has_crs(soil_stack) && safe_crs(soil_stack) != "EPSG:4326") terra::project(pts, soil_stack) else pts

  vals <- terra::extract(soil_stack, pts_proj, ID = FALSE) %>% tibble::as_tibble()
  vals <- vals %>% mutate(across(everything(), ~ as.numeric(.) / 10))

  out <- bind_cols(coords_df, vals) %>%
    mutate(
      source_id = "soilgrids",
      source_name = "SoilGrids",
      native_depth = "0-5 cm",
      native_resolution_m = 250,
      extraction_method = "Point extraction from SoilGrids WebDAV/VRT layers; g/kg converted to percent by dividing by 10",
      data_source = "ISRIC SoilGrids 2.0 WebDAV/VRT layers sand_0-5cm_mean.vrt, silt_0-5cm_mean.vrt, clay_0-5cm_mean.vrt; values converted from g/kg to %.",
      sand_pct = as.numeric(sand_pct),
      silt_pct = as.numeric(silt_pct),
      clay_pct = as.numeric(clay_pct),
      qa_flag = mapply(format_fraction_sum_flag, sand_pct, silt_pct, clay_pct),
      texture_class_usda = mapply(classify_usda_texture, sand_pct, silt_pct, clay_pct),
      error_message = NA_character_
    )

  out
}

# ------------------------------------------------------------------------------
# HWSD v2 backend (cross-platform SQLite path; no MDB/ODBC required)
# ------------------------------------------------------------------------------

ensure_hwsd_assets <- function(progress_callback = NULL) {
  dir.create(CACHE_HWSD, recursive = TRUE, showWarnings = FALSE)

  bil_path <- find_first_file(CACHE_HWSD, "\\.bil$")
  if (is.na(bil_path)) {
    if (is.function(progress_callback)) {
      progress_callback(label = "HWSD v2.0", detail = "Downloading raster package")
    }
    zip_path <- file.path(CACHE_HWSD, "HWSD2_RASTER.zip")
    if (!file.exists(zip_path) || is.na(file.info(zip_path)$size) || file.info(zip_path)$size <= 0) {
      download_file_robust(HWSD_URLS$raster_zip, zip_path)
    }
    safe_unzip(zip_path, CACHE_HWSD)
    bil_path <- find_first_file(CACHE_HWSD, "\\.bil$")
  }

  if (is.na(bil_path)) {
    stop("HWSD raster package was downloaded but no .bil raster was found after unzip.", call. = FALSE)
  }

  sqlite_path <- find_first_file(CACHE_HWSD, "HWSD2\\.sqlite$")
  if (is.na(sqlite_path)) {
    if (is.function(progress_callback)) {
      progress_callback(label = "HWSD v2.0", detail = "Downloading SQLite attribute database")
    }
    sqlite_path <- file.path(CACHE_HWSD, "HWSD2.sqlite")
    download_file_robust(HWSD_URLS$sqlite, sqlite_path)
  }

  list(
    bil = bil_path,
    sqlite = sqlite_path
  )
}

query_hwsd_layers <- function(sqlite_path, smu_ids, layer = "D1") {
  smu_ids <- unique(as.integer(smu_ids[is.finite(smu_ids)]))
  if (length(smu_ids) == 0) return(tibble())

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqlite_path)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

  id_sql <- paste(smu_ids, collapse = ",")
  sql <- sprintf(
    paste(
      "select HWSD2_SMU_ID, SEQUENCE, SHARE, TOPDEP, BOTDEP,",
      "SAND, SILT, CLAY, TEXTURE_USDA",
      "from HWSD2_LAYERS",
      "where LAYER = '%s' and HWSD2_SMU_ID in (%s)",
      "order by HWSD2_SMU_ID, SEQUENCE"
    ),
    layer,
    id_sql
  )

  DBI::dbGetQuery(con, sql) %>% tibble::as_tibble()
}

aggregate_hwsd_components <- function(layer_df) {
  if (nrow(layer_df) == 0) {
    return(tibble(
      HWSD2_SMU_ID = integer(),
      sand_pct = numeric(),
      silt_pct = numeric(),
      clay_pct = numeric(),
      component_count = integer(),
      valid_component_count = integer(),
      share_sum_valid = numeric(),
      qa_flag = character()
    ))
  }

  layer_df %>%
    mutate(
      HWSD2_SMU_ID = as.integer(HWSD2_SMU_ID),
      SEQUENCE = as.integer(SEQUENCE),
      SHARE = as.numeric(SHARE),
      SAND = as.numeric(SAND),
      SILT = as.numeric(SILT),
      CLAY = as.numeric(CLAY),
      valid_triplet = is.finite(SHARE) & SHARE > 0 &
        is.finite(SAND) & SAND >= 0 &
        is.finite(SILT) & SILT >= 0 &
        is.finite(CLAY) & CLAY >= 0
    ) %>%
    group_by(HWSD2_SMU_ID) %>%
    summarise(
      component_count = n_distinct(SEQUENCE[is.finite(SEQUENCE)]),
      valid_component_count = sum(valid_triplet),
      share_sum_valid = sum(SHARE[valid_triplet], na.rm = TRUE),
      sand_pct = weighted_mean_safe(SAND[valid_triplet], SHARE[valid_triplet]),
      silt_pct = weighted_mean_safe(SILT[valid_triplet], SHARE[valid_triplet]),
      clay_pct = weighted_mean_safe(CLAY[valid_triplet], SHARE[valid_triplet]),
      qa_flag = case_when(
        valid_component_count == 0 ~ "no_valid_components",
        share_sum_valid < 99.5 ~ "partial_share_weight",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      qa_flag = mapply(coalesce_flag, qa_flag, mapply(format_fraction_sum_flag, sand_pct, silt_pct, clay_pct))
    )
}

extract_hwsd <- function(coords_df, progress_callback = NULL) {
  assets <- ensure_hwsd_assets(progress_callback = progress_callback)

  if (is.function(progress_callback)) {
    progress_callback(label = "HWSD v2.0", detail = "Reading raster and SQLite attributes")
  }

  hwsd_raster <- terra::rast(assets$bil)
  if (!has_crs(hwsd_raster)) {
    try(terra::crs(hwsd_raster) <- "EPSG:4326", silent = TRUE)
  }

  pts <- terra::vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  ids <- terra::extract(hwsd_raster, pts, ID = FALSE) %>% as.data.frame()
  id_col <- names(ids)[[1]]

  base <- bind_cols(coords_df, tibble(smu_id = as.integer(ids[[id_col]])))
  layer_df <- query_hwsd_layers(assets$sqlite, base$smu_id, layer = "D1")
  agg_df <- aggregate_hwsd_components(layer_df)

  out <- base %>%
    left_join(agg_df, by = c("smu_id" = "HWSD2_SMU_ID")) %>%
    mutate(
      source_id = "hwsd_v2",
      source_name = "HWSD v2.0",
      native_depth = "D1 (0-20 cm)",
      native_resolution_m = 1000,
      extraction_method = "Raster SMU lookup (30 arc-second) + SHARE-weighted aggregation across all valid D1 components in HWSD2_LAYERS",
      data_source = paste(
        "HWSD v2.0 raster from FAO/GAEZ plus SQLite attribute database hosted by ISRIC;",
        "D1 layer (0-20 cm); SHARE-weighted aggregation across components using HWSD2_LAYERS."
      ),
      texture_class_usda = mapply(classify_usda_texture, sand_pct, silt_pct, clay_pct),
      error_message = NA_character_
    ) %>%
    mutate(
      qa_flag = ifelse(is.na(smu_id), coalesce_flag(qa_flag, "no_smu_id"), qa_flag)
    )

  out
}

# ------------------------------------------------------------------------------
# Combined extraction
# ------------------------------------------------------------------------------

placeholder_backend_result <- function(coords_df, source_id, source_name, native_depth, native_resolution_m, extraction_method, data_source, error_message) {
  coords_df %>%
    mutate(
      source_id = source_id,
      source_name = source_name,
      native_depth = native_depth,
      native_resolution_m = native_resolution_m,
      extraction_method = extraction_method,
      data_source = data_source,
      sand_pct = NA_real_,
      silt_pct = NA_real_,
      clay_pct = NA_real_,
      texture_class_usda = NA_character_,
      qa_flag = "source_error",
      error_message = as.character(error_message)
    )
}

run_all_backends <- function(coords_df, progress_callback = NULL) {
  results <- list()

  results[["openlandmap"]] <- tryCatch(
    extract_openlandmap(coords_df, progress_callback = progress_callback),
    error = function(e) {
      placeholder_backend_result(
        coords_df = coords_df,
        source_id = "openlandmap",
        source_name = "OpenLandMap",
        native_depth = "b0..0cm",
        native_resolution_m = 250,
        extraction_method = "Point-in-raster-cell extraction from cached GeoTIFFs",
        data_source = "OpenLandMap / LandGIS v0.2 topsoil rasters from Zenodo.",
        error_message = conditionMessage(e)
      )
    }
  )

  results[["soilgrids"]] <- tryCatch(
    extract_soilgrids(coords_df, progress_callback = progress_callback),
    error = function(e) {
      placeholder_backend_result(
        coords_df = coords_df,
        source_id = "soilgrids",
        source_name = "SoilGrids",
        native_depth = "0-5 cm",
        native_resolution_m = 250,
        extraction_method = "Point extraction from SoilGrids WebDAV/VRT layers",
        data_source = "ISRIC SoilGrids 2.0 WebDAV/VRT layers.",
        error_message = conditionMessage(e)
      )
    }
  )

  results[["hwsd_v2"]] <- tryCatch(
    extract_hwsd(coords_df, progress_callback = progress_callback),
    error = function(e) {
      placeholder_backend_result(
        coords_df = coords_df,
        source_id = "hwsd_v2",
        source_name = "HWSD v2.0",
        native_depth = "D1 (0-20 cm)",
        native_resolution_m = 1000,
        extraction_method = "Raster SMU lookup + SHARE-weighted aggregation",
        data_source = "HWSD v2.0 raster + ISRIC-hosted SQLite attribute database.",
        error_message = conditionMessage(e)
      )
    }
  )

  long_df <- safe_rbind(results) %>%
    mutate(
      source_id = factor(source_id, levels = BACKEND_ORDER),
      source_name = factor(source_name, levels = BACKEND_LABELS[BACKEND_ORDER]),
      fraction_total = sand_pct + silt_pct + clay_pct
    ) %>%
    arrange(source_id, coord_id)

  long_df
}

join_results_to_sites <- function(sites_df, coords_long_df) {
  sites_df %>%
    left_join(coords_long_df, by = c("coord_id", "lat", "lon")) %>%
    mutate(
      source_id = as.character(source_id),
      source_name = as.character(source_name)
    ) %>%
    select(
      site_id, site_label, lat, lon, coord_id,
      source_id, source_name,
      sand_pct, silt_pct, clay_pct, texture_class_usda,
      native_depth, native_resolution_m,
      extraction_method, qa_flag, error_message, data_source, fraction_total
    ) %>%
    arrange(site_label, factor(source_id, levels = BACKEND_ORDER))
}

make_wide_export <- function(long_sites_df) {
  long_sites_df %>%
    select(
      site_id, lat, lon, source_id,
      sand_pct, silt_pct, clay_pct, texture_class_usda,
      native_depth, native_resolution_m, qa_flag, error_message
    ) %>%
    mutate(source_id = factor(source_id, levels = BACKEND_ORDER)) %>%
    pivot_wider(
      names_from = source_id,
      values_from = c(
        sand_pct, silt_pct, clay_pct, texture_class_usda,
        native_depth, native_resolution_m, qa_flag, error_message
      ),
      names_glue = "{source_id}_{.value}"
    ) %>%
    arrange(site_id, lat, lon)
}

make_long_export <- function(long_sites_df) {
  long_sites_df %>%
    transmute(
      site_id, lat, lon,
      source_id, source_name,
      sand_pct, silt_pct, clay_pct,
      texture_class_usda,
      native_depth, native_resolution_m,
      extraction_method, qa_flag, error_message,
      data_source
    ) %>%
    arrange(site_id, factor(source_id, levels = BACKEND_ORDER))
}

# ------------------------------------------------------------------------------
# UI helpers
# ------------------------------------------------------------------------------

metric_card <- function(title, sand, silt, clay, texture_class, depth, resolution, qa_flag = NA_character_, error_message = NA_character_, accent_class = "") {
  tags$div(
    class = paste("source-metric-card", accent_class),
    tags$div(class = "source-metric-title", title),
    if (!is.na(error_message) && nzchar(error_message)) {
      tags$div(
        class = "alert alert-warning py-2 px-3 mb-3",
        tags$strong("Source error: "),
        error_message
      )
    } else {
      tagList(
        tags$div(class = "metric-grid",
                 tags$div(class = "metric-item",
                          tags$div(class = "metric-label", "Sand"),
                          tags$div(class = "metric-value", format_pct(sand)),
                          tags$div(class = "metric-unit", "%")),
                 tags$div(class = "metric-item",
                          tags$div(class = "metric-label", "Silt"),
                          tags$div(class = "metric-value", format_pct(silt)),
                          tags$div(class = "metric-unit", "%")),
                 tags$div(class = "metric-item",
                          tags$div(class = "metric-label", "Clay"),
                          tags$div(class = "metric-value", format_pct(clay)),
                          tags$div(class = "metric-unit", "%"))
        ),
        tags$div(class = "source-subtable",
                 tags$div(tags$strong("Texture class: "), texture_class %||% "NA"),
                 tags$div(tags$strong("Native depth: "), depth %||% "NA"),
                 tags$div(tags$strong("Native resolution: "), ifelse(is.na(resolution), "NA", paste0(resolution, " m"))),
                 if (!is.na(qa_flag) && nzchar(qa_flag)) tags$div(tags$strong("QA flag: "), qa_flag)
        )
      )
    }
  )
}

selected_site_source_cards_ui <- function(site_df) {
  if (nrow(site_df) == 0) {
    return(tags$div(class = "alert alert-info", "No site selected yet."))
  }

  rows <- lapply(BACKEND_ORDER, function(src) {
    row <- site_df %>% filter(source_id == src)
    if (nrow(row) == 0) row <- tibble(
      source_name = BACKEND_LABELS[[src]],
      sand_pct = NA_real_, silt_pct = NA_real_, clay_pct = NA_real_,
      texture_class_usda = NA_character_, native_depth = NA_character_,
      native_resolution_m = NA_real_, qa_flag = NA_character_, error_message = "No record returned."
    )

    accent <- switch(
      src,
      openlandmap = "openlandmap-card",
      soilgrids = "soilgrids-card",
      hwsd_v2 = "hwsd-card",
      ""
    )

    tags$div(
      class = "col-lg-4 col-md-6 mb-3",
      metric_card(
        title = as.character(row$source_name[[1]] %||% BACKEND_LABELS[[src]]),
        sand = row$sand_pct[[1]],
        silt = row$silt_pct[[1]],
        clay = row$clay_pct[[1]],
        texture_class = row$texture_class_usda[[1]],
        depth = row$native_depth[[1]],
        resolution = row$native_resolution_m[[1]],
        qa_flag = row$qa_flag[[1]],
        error_message = row$error_message[[1]],
        accent_class = accent
      )
    )
  })

  tags$div(class = "row", rows)
}

source_summary_table <- function(site_df) {
  if (nrow(site_df) == 0) {
    return(placeholder_table("No site selected yet."))
  }

  out <- site_df %>%
    transmute(
      Source = source_name,
      `Sand (%)` = round(sand_pct, 1),
      `Silt (%)` = round(silt_pct, 1),
      `Clay (%)` = round(clay_pct, 1),
      `Texture class` = texture_class_usda,
      `Native depth` = native_depth,
      `Resolution (m)` = native_resolution_m,
      `QA flag` = qa_flag,
      `Error` = error_message
    )

  DT::datatable(
    out,
    rownames = FALSE,
    options = list(dom = "tip", paging = FALSE, autoWidth = TRUE, scrollX = TRUE)
  )
}

placeholder_table <- function(message) {
  DT::datatable(
    data.frame(Message = message),
    rownames = FALSE,
    options = list(dom = "t", paging = FALSE)
  )
}

source_metadata_table <- function() {
  meta <- tibble::tribble(
    ~Source, ~Native_top_layer, ~Nominal_resolution_m, ~Core_access_path,
    "OpenLandMap", "b0..0cm", 250, "Cached GeoTIFFs from Zenodo",
    "SoilGrids", "0-5 cm mean", 250, "Remote WebDAV/VRT",
    "HWSD v2.0", "D1 (0-20 cm)", 1000, "Raster SMU grid + SQLite attribute database"
  )

  DT::datatable(
    meta,
    rownames = FALSE,
    options = list(dom = "t", paging = FALSE, autoWidth = TRUE)
  )
}

# ------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------

ui <- page_sidebar(
  title = div("Soil Texture Comparison Explorer"),
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = TRUE,
  sidebar = sidebar(
    width = 360,
    open = "desktop",
    tags$div(
      class = "sidebar-section",
      tags$h5("Upload"),
      fileInput(
        inputId = "sites_file",
        label = "CSV with columns: site_id, lat, lon",
        accept = c(".csv", "text/csv", "text/comma-separated-values,text/plain")
      ),
      actionButton("run_lookup", "Run / refresh extraction", class = "btn btn-primary w-100 mt-2"),
      tags$div(
        class = "small text-muted mt-2",
        "Exact duplicate rows are removed automatically. Different site IDs may share the same coordinates; those are extracted once and joined back to all sites."
      )
    ),
    tags$hr(),
    tags$div(
      class = "sidebar-section",
      tags$h5("Selected site"),
      selectInput(
        inputId = "site_choice",
        label = NULL,
        choices = character(0),
        selected = NULL
      ),
      tags$div(class = "small text-muted", "The comparison cards and map update when you choose a site.")
    ),
    tags$hr(),
    tags$div(
      class = "sidebar-section",
      tags$h5("Downloads"),
      downloadButton("download_wide", "Download wide CSV", class = "btn btn-outline-primary w-100 mb-2"),
      downloadButton("download_long", "Download long CSV", class = "btn btn-outline-secondary w-100")
    ),
    tags$hr(),
    tags$div(
      class = "sidebar-section",
      tags$h5("Session summary"),
      uiOutput("session_notes")
    )
  ),
  tags$style(HTML("
    .bslib-sidebar-layout > .main { overflow: auto; }
    .app-hero { margin-bottom: 1rem; }
    .app-hero h4 { margin-bottom: .35rem; }
    .soft-card {
      background: #fff;
      border: 1px solid rgba(0,0,0,.08);
      border-radius: 14px;
      box-shadow: 0 4px 16px rgba(0,0,0,.04);
      padding: 1rem 1.1rem;
      margin-bottom: 1rem;
    }
    .source-metric-card {
      background: #fff;
      border: 1px solid rgba(0,0,0,.08);
      border-radius: 16px;
      box-shadow: 0 4px 16px rgba(0,0,0,.05);
      padding: 1rem 1rem 0.9rem 1rem;
      min-height: 100%;
    }
    .source-metric-title {
      font-weight: 700;
      font-size: 1.05rem;
      margin-bottom: 0.85rem;
    }
    .metric-grid {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: .7rem;
      margin-bottom: .85rem;
    }
    .metric-item {
      background: rgba(13,110,253,.05);
      border-radius: 12px;
      padding: .75rem .6rem;
      text-align: center;
    }
    .openlandmap-card .metric-item { background: rgba(25,135,84,.08); }
    .soilgrids-card .metric-item { background: rgba(13,110,253,.08); }
    .hwsd-card .metric-item { background: rgba(220,53,69,.08); }
    .metric-label {
      color: #52606d;
      font-size: .85rem;
      margin-bottom: .2rem;
    }
    .metric-value {
      font-size: 1.5rem;
      font-weight: 700;
      line-height: 1.1;
    }
    .metric-unit {
      color: #52606d;
      font-size: .8rem;
    }
    .source-subtable div { margin-bottom: .22rem; }
    .sidebar-section h5 { margin-bottom: .7rem; }
    .method-list li { margin-bottom: .35rem; }
    .status-pill {
      display: inline-block;
      padding: .25rem .55rem;
      border-radius: 999px;
      background: rgba(13,110,253,.08);
      color: #0d6efd;
      font-size: .84rem;
      margin-right: .4rem;
      margin-bottom: .4rem;
    }
    .map-wrap { min-height: 540px; }
  ")),
  navset_card_tab(
    height = NULL,

    nav_panel(
      "Selected site",
      tags$div(
        class = "app-hero soft-card",
        tags$h4("Three-source topsoil comparison"),
        tags$p(
          class = "mb-0",
          "Every query runs OpenLandMap, SoilGrids, and HWSD v2.0. Results are shown using each source's native top layer and native spatial support."
        )
      ),
      uiOutput("backend_alerts"),
      tags$div(class = "soft-card", uiOutput("selected_site_heading")),
      uiOutput("selected_site_cards"),
      tags$div(class = "soft-card map-wrap", leafletOutput("site_map", height = 520)),
      tags$div(class = "soft-card", DTOutput("selected_site_table"))
    ),

    nav_panel(
      "Results preview",
      tags$div(class = "soft-card", DTOutput("preview_wide_table"))
    ),

    nav_panel(
      "Comparison",
      tags$div(class = "soft-card", tags$h5("Per-source long table"), DTOutput("comparison_long_table")),
      tags$div(class = "soft-card", tags$h5("Source metadata"), source_metadata_table())
    ),

    nav_panel(
      "Methods & downloads",
      tags$div(
        class = "soft-card",
        tags$h5("Workflow summary"),
        tags$ul(
          class = "method-list",
          tags$li("The app cleans the uploaded CSV and removes exact duplicate rows."),
          tags$li("Unique coordinates are extracted once from each backend and then joined back to all matching site IDs."),
          tags$li("OpenLandMap uses cached 250 m topsoil GeoTIFFs (b0..0cm)."),
          tags$li("SoilGrids uses the 0-5 cm WebDAV/VRT layers and converts g/kg to percent by dividing by 10."),
          tags$li("HWSD v2.0 uses the 30 arc-second SMU raster plus the SQLite attribute database and aggregates all valid D1 components with SHARE weights."),
          tags$li("Texture class is derived consistently across all three backends from the extracted sand, silt, and clay fractions.")
        )
      ),
      tags$div(
        class = "soft-card",
        tags$h5("Assumptions and cautions"),
        tags$ul(
          class = "method-list",
          tags$li("The three products do not use identical topsoil depth support: OpenLandMap b0..0cm, SoilGrids 0-5 cm, HWSD D1 0-20 cm."),
          tags$li("Raster extraction uses the cell intersecting the point. On regular grids this is equivalent to the nearest cell-center lookup except on exact cell boundaries."),
          tags$li("Negative HWSD fraction values are treated as missing and excluded from the SHARE-weighted aggregation."),
          tags$li("When a backend fails, the app keeps the other sources and records the error in the outputs."),
          tags$li("QA flags identify source failures, missing SMU IDs, partial valid SHARE coverage, and fraction sums that do not total approximately 100%.")
        )
      ),
      tags$div(
        class = "soft-card",
        tags$h5("Download outputs"),
        tags$p("Use the buttons in the sidebar to export either a wide comparison CSV or a long provenance-rich CSV.")
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Server
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  prepared_sites <- reactiveVal(NULL)
  results_long <- reactiveVal(NULL)
  results_wide <- reactiveVal(NULL)

  output$session_notes <- renderUI({
    notes <- c(
      paste("Cache directory:", CACHE_DIR),
      "OpenLandMap and HWSD assets are cached locally after first download.",
      "SoilGrids is read remotely through its published VRT layers."
    )
    tags$div(lapply(notes, function(x) tags$div(class = "small mb-1", x)))
  })

  observeEvent(input$run_lookup, {
    req(input$sites_file)

    withProgress(message = "Preparing input and running soil queries", value = 0, {
      incProgress(0.05, detail = "Reading and validating uploaded CSV")
      prep <- prepare_sites(input$sites_file$datapath)
      prepared_sites(prep)

      updateSelectInput(
        session,
        "site_choice",
        choices = setNames(prep$sites$site_label, prep$sites$site_label),
        selected = prep$sites$site_label[[1]]
      )

      coords_df <- prep$coords

      progress_callback <- function(label = NULL, detail = NULL) {
        shiny::setProgress(
          message = label %||% "Working",
          detail = detail %||% NULL
        )
      }

      long_coords <- run_all_backends(coords_df, progress_callback = progress_callback)
      long_sites <- join_results_to_sites(prep$sites, long_coords)
      wide_sites <- make_wide_export(long_sites)

      results_long(long_sites)
      results_wide(wide_sites)

      incProgress(1, detail = "Finished")
    })
  })

  selected_site_df <- reactive({
    req(prepared_sites(), results_long(), input$site_choice)
    results_long() %>% filter(site_label == input$site_choice)
  })

  output$backend_alerts <- renderUI({
    req(results_long())
    errs <- results_long() %>%
      filter(!is.na(error_message) & nzchar(error_message)) %>%
      distinct(source_name, error_message)

    if (nrow(errs) == 0) return(NULL)

    tags$div(
      class = "soft-card",
      tags$h5("Backend notices"),
      lapply(seq_len(nrow(errs)), function(i) {
        tags$div(
          class = "alert alert-warning mb-2",
          tags$strong(as.character(errs$source_name[[i]])), ": ",
          errs$error_message[[i]]
        )
      })
    )
  })

  output$selected_site_heading <- renderUI({
    req(selected_site_df())
    site_df <- selected_site_df()
    first_row <- site_df[1, , drop = FALSE]

    tags$div(
      tags$h5(first_row$site_id[[1]], class = "mb-2"),
      tags$div(
        class = "status-pill",
        paste("Latitude", format_num(first_row$lat[[1]], digits = 7))
      ),
      tags$div(
        class = "status-pill",
        paste("Longitude", format_num(first_row$lon[[1]], digits = 7))
      ),
      tags$div(
        class = "status-pill",
        paste("Compared sources", nrow(site_df))
      )
    )
  })

  output$selected_site_cards <- renderUI({
    req(selected_site_df())
    selected_site_source_cards_ui(selected_site_df())
  })

  output$selected_site_table <- renderDT({
    req(selected_site_df())
    source_summary_table(selected_site_df())
  })

  output$preview_wide_table <- renderDT({
    req(results_wide())
    DT::datatable(
      results_wide(),
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })

  output$comparison_long_table <- renderDT({
    req(results_long())
    DT::datatable(
      make_long_export(results_long()),
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })

  output$site_map <- renderLeaflet({
    req(selected_site_df())
    site_df <- selected_site_df()
    lat <- site_df$lat[[1]]
    lon <- site_df$lon[[1]]
    site_name <- site_df$site_id[[1]]

    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addCircleMarkers(
        lng = lon,
        lat = lat,
        radius = 7,
        stroke = TRUE,
        weight = 2,
        fillOpacity = 0.9,
        popup = paste0(
          "<strong>", htmlEscape(site_name), "</strong><br/>",
          "Lat: ", format_num(lat, digits = 7), "<br/>",
          "Lon: ", format_num(lon, digits = 7)
        )
      ) %>%
      setView(lng = lon, lat = lat, zoom = 8)
  })

  output$download_wide <- downloadHandler(
    filename = function() {
      paste0("soil_texture_comparison_wide_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results_wide())
      readr::write_csv(results_wide(), file, na = "")
    }
  )

  output$download_long <- downloadHandler(
    filename = function() {
      paste0("soil_texture_comparison_long_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results_long())
      readr::write_csv(make_long_export(results_long()), file, na = "")
    }
  )
}

shinyApp(ui, server)
