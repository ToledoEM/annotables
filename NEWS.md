# v0.3.0

- Recipes pinned to `jun2026.archive.ensembl.org` (Ensembl 116). Bundled tables
  are built from Ensembl 114.
- Fixed `HTTP 405 Method Not Allowed` from `getBM()`: `www.ensembl.org`
  308-redirects BioMart to a `GET`-only endpoint that rejects biomaRt's POST.
- `get_data()` re-pins `mart@host` after `useMart()`, which overwrites it with
  `www.ensembl.org` from the BioMart registry.
- Refactored `build.Rmd` for resilience against dropouts and cut connections
  mid-build. Each table is now downloaded, saved and logged individually instead
  of saving everything at the end, so an interrupted run keeps what finished.
  Completed tables are recorded in `data-raw/build_progress.log` and skipped on
  rerun; the log is deleted once all tables are present.
- Added retry with backoff on 502/503/504/429, a 3s inter-request delay, and a
  per-host mart cache.
- Retries now also cover `"is not valid"` and `"unavailable"`: an overloaded
  archive serves an HTML error page with HTTP 200, which biomaRt misreports as
  an invalid dataset name.
- Set `options(timeout = 600)`; biomaRt uses R's 60s default, and a cold archive
  host was measured taking 51.7s to return the registry alone.
- `tx2gene` no longer overwrites `recipes` in place.
- `load_recipes` globs `*.yml` only.
- Added `data-raw/switch-host.sh` and `data-raw/recipes/README.md`.

# v0.2.0

- Added host to data retrieval function
- Updated all recipes to use https
