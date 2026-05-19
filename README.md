# Soil Texture Comparison Explorer

This Shiny app reads a CSV with columns `site_id`, `lat`, and `lon`, then runs **all three backends automatically** for every site:

1. **OpenLandMap / LandGIS v0.2**
   - 250 m topsoil texture layers
   - App uses the `b0..0cm` standard-depth band
2. **SoilGrids 2.0 / ISRIC**
   - 250 m texture predictions
   - App uses the `0-5 cm` mean layers
   - SoilGrids stores texture fractions in `g/kg`; the app converts to `%` by dividing by `10`
3. **HWSD v2.0 / FAO-IIASA**
   - ~1 km raster of soil mapping units plus linked attribute database
   - App uses **D1 (0-20 cm)**
   - App aggregates **all soil components in the mapping unit using `SHARE` as weights**

## Main features

- Upload CSV with `site_id`, `lat`, `lon`
- Exact duplicate rows removed automatically
- Repeated coordinates extracted once and joined back to all matching site IDs
- Automatic three-source comparison for every run
- Selected-site map view
- Results preview tab (wide comparison table)
- Comparison tab (selected-site table + long provenance table)
- Download buttons for:
  - wide comparison CSV
  - long provenance CSV

## Expected input

```csv
site_id,lat,lon
SITE_A,40.123,-105.456
SITE_B,39.987,-104.321
```

A column named `long` is also accepted and renamed internally to `lon`.

## Installation

Required R packages:

```r
install.packages(c(
  "shiny", "bslib", "leaflet", "terra", "readr", "dplyr",
  "glue", "htmltools", "DT", "tibble", "tidyr"
))
```

Optional helpers:

```r
install.packages(c("curl", "DBI", "odbc"))
```

## HWSD requirement

HWSD is distributed as a raster package plus a Microsoft Access `.mdb` database.
The app can currently read that database through **either**:

- an ODBC Access-compatible driver (common on Windows), **or**
- `mdbtools` command-line utilities (`mdb-tables`, `mdb-export`), often easiest on Linux/macOS.

If neither is available, the OpenLandMap and SoilGrids backends can still succeed, but HWSD will report a source-specific error.

## Running the app

```r
shiny::runApp("app.R")
```

## Outputs

### Wide comparison CSV
One row per site, with side-by-side columns for all three sources.

### Long provenance CSV
One row per `site x source`, including:
- fractions
- derived USDA texture class
- source status and message
- native depth
- native resolution
- extraction method
- data source and license fields

## Important interpretation notes

- The three sources are compared at their **native top layers**; they are **not depth-harmonized**.
- OpenLandMap, SoilGrids, and HWSD differ in spatial support, modeling method, and provenance.
- HWSD is a **mapping-unit** product, not a continuous raster prediction like the two 250 m products.
- Fractions are exported as extracted; the USDA texture class is derived from a normalized triplet only for classification purposes.

## Cache behavior

The app caches downloaded files under an R cache directory and reuses them across sessions.
You can override the cache root with:

```r
Sys.setenv(SOIL_CACHE_DIR = "/path/to/cache")
```

## Notes for publication-oriented use

Before public release, it is worth validating a handful of fixed benchmark coordinates manually against each source and reviewing downstream licensing constraints for any redistributed outputs.
