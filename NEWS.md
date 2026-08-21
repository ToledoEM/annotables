# v0.3.0

## Data

All tables rebuilt from **Ensembl 116** (previously Ensembl 114). Recipes are
pinned to the June 2026 archive, `jun2026.archive.ensembl.org`.

Recipes select data by species, and BioMart always returns whatever assembly is
current for that species, so several datasets moved to a new assembly with this
release:

| Dataset | Species | Assembly now served |
|---|---|---|
| `grch38` | Human | GRCh38 |
| `grch37` | Human | GRCh37 |
| `grcm39` | Mouse | GRCm39 |
| `grcr8` | Rat | GRCr8 |
| `grcg7b` | Chicken | bGalGal1.mat.broiler.GRCg7b |
| `wbcel235` | Roundworm | WBcel235 |
| `bdgp654` | Fruitfly | BDGP6.54 |
| `grcz11` | Zebrafish | GRCz11 |
| `mmul10` | Macaque | Mmul_10 |
| `roscfam1` | Dog | ROS_Cfam_1.0 |
| `sscrofa111` | Pig | Sscrofa11.1 |

## Breaking: datasets renamed after their assemblies

Datasets are named after the **genome assembly** from now on, not the species.

A species-based name stops describing its own contents as soon as Ensembl
switches assembly, and nothing fails to make that visible — the data simply
changes underneath the name. Naming after the assembly means a stale name is
obvious.

Eight datasets were renamed. The old names were removed rather than deprecated,
so existing code errors instead of quietly returning a different assembly:

| Old | New | Species | Assembly |
|---|---|---|---|
| `rnor6` | `grcr8` | Rat | GRCr8 |
| `grcm38` | `grcm39` | Mouse | GRCm39 |
| `galgal5` | `grcg7b` | Chicken | bGalGal1.mat.broiler.GRCg7b |
| `mmul801` | `mmul10` | Macaque | Mmul_10 |
| `cfamiliaris` | `roscfam1` | Dog | ROS_Cfam_1.0 |
| `bdgp6` | `bdgp654` | Fruitfly | BDGP6.54 |
| `drerio` | `grcz11` | Zebrafish | GRCz11 |
| `sscrofa` | `sscrofa111` | Pig | Sscrofa11.1 |

This applies to the `_tx2gene` tables too. `grch37`, `grch38` and `wbcel235`
were already assembly-named and are unchanged.

## Build

- Fixed `HTTP 405 Method Not Allowed` from `getBM()`. `www.ensembl.org`
  308-redirects BioMart to a `GET`-only endpoint, which rejects the POST biomaRt
  uses to submit query XML.
- `get_mart()` constructs the `Mart` object directly and calls `useDataset()` on
  it, instead of `useMart()`. The BioMart registry advertises
  `host="www.ensembl.org"` for every mart, so `useMart()` would send its
  dataset-list request to the broken host before returning; patching
  `mart@host` afterwards was too late.
- Rebuilt for resilience against dropped connections. Each table is downloaded,
  saved and logged individually rather than saving everything at the end, so an
  interrupted run keeps what finished. Completed tables are recorded in
  `data-raw/build_progress.log` and skipped on rerun; the log is removed once
  all tables are present.
- Retries with backoff on 502/503/504/429, plus a 3s delay between requests and
  a per-host mart cache.
- Retries also cover `"is not valid"`, `"unavailable"` and
  `"incompatible version"`. An overloaded archive serves an HTML error page with
  HTTP 200, which biomaRt parses as data and misreports as a bad dataset name or
  a version mismatch.
- Set `options(timeout = 600)`. biomaRt uses R's 60s default, and a cold archive
  host was measured taking 51.7s to return the registry alone.
- `tx2gene` no longer overwrites `recipes` in place.
- `load_recipes` globs `*.yml` only.
- Added `data-raw/switch-host.sh` and `data-raw/recipes/README.md`.

# v0.2.0

- Added host to data retrieval function
- Updated all recipes to use https
