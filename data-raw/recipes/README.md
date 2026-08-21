# Recipe hosts — why they point at an Ensembl archive

**Do not change `host:` back to `https://www.ensembl.org`.** That host is what
breaks the build. Read this first.

## The symptom

`data-raw/build.Rmd` dies in the `[gene_annotation_tables]` chunk:

```
Error in `httr2::req_perform()`:
! HTTP 405 Method Not Allowed.
  4. biomaRt:::.submitQueryXML(...)
```

## The cause

Ensembl redirects BioMart on `www.ensembl.org` to a release archive, and that
archive endpoint accepts `GET` only. biomaRt always POSTs its query XML, so
every `getBM()` call is rejected:

```
$ curl -D - -X POST https://www.ensembl.org/biomart/martservice
HTTP/2 405
allow: GET

$ curl -D - "https://www.ensembl.org/biomart/martservice?type=registry"
HTTP/2 308
location: https://jun2026.archive.ensembl.org/biomart/martservice?type=registry
```

Posting to the archive host directly works fine — HTTP 200, Ensembl Genes 116.

## The trap: `useMart()` overwrites the host

Setting `host:` in these YAML files is **necessary but not sufficient**.
`useMart()` fetches the BioMart registry, and the registry advertises

```xml
<MartURLLocation ... host="www.ensembl.org" path="/biomart/martservice" port="80" ... />
```

biomaRt then silently rewrites the mart object's host to that value. Even when
called with `host = "https://jun2026.archive.ensembl.org"`, the returned mart
comes back as:

```
mart@host  #> https://www.ensembl.org:443/biomart/martservice
```

…and the next `getBM()` 405s again.

So the fix has two halves, and **both are required**:

1. `host:` in each recipe here points at the archive.
2. `get_data()` in `data-raw/build.Rmd` re-pins `mart@host` from `recipe$host`
   immediately after `useMart()` returns, undoing the registry override.

Change one without the other and the 405 comes straight back.

## Exception: `grch37.yml`

`grch37.yml` stays on `https://grch37.ensembl.org`. That is a separate,
permanently maintained archive for the GRCh37 assembly and it answers POST
normally. It is deliberately excluded from host switching.

## Archives are release-pinned

`jun2026.archive.ensembl.org` serves **Ensembl release 116**. Pinning an archive
means the build is reproducible but frozen at that release — moving to a newer
release is an explicit decision, not a silent drift.

If Ensembl finishes its transition and `www.ensembl.org` starts accepting POST
again, or you want a different release, use the helper rather than editing 10
files by hand:

```bash
data-raw/switch-host.sh --dry-run https://www.ensembl.org   # preview
data-raw/switch-host.sh https://www.ensembl.org             # apply
data-raw/switch-host.sh jun2026                             # re-pin to the archive
```

Because `get_data()` reads `recipe$host`, the script is the whole switch — no R
code needs editing to change hosts.

## Unrelated failure: HTTP 403

`useast.ensembl.org` and `asia.ensembl.org` return **403 Forbidden** when you
query them too quickly. That is IP rate limiting, not this bug — rate limits
never return 405.

Ensembl issues **no API keys**; BioMart is open and unauthenticated, so no
credential will clear a 403. Wait 15–30 minutes. Do not add automatic retries,
which prolong the block.

## Rebuilding locally

### Run it from `data-raw/`, not the repo root

`build.Rmd` sets its own working directory in the setup chunk:

```r
opts_knit$set(root.dir = normalizePath(".."))
```

That `".."` is resolved **relative to wherever R was started**, and because it
runs inside the document it **overrides** any `knit_root_dir` passed to
`rmarkdown::render()`. RStudio knits with the Rmd's own directory as the
starting point, so `".."` lands on the repo root and everything works. From a
terminal started at the repo root, `".."` resolves to the *parent of the repo*,
where `data-raw/recipes` does not exist.

This fails **silently and destructively**: `dir()` finds no `.yml` files,
`recipes` becomes an empty list, `lapply(recipes, get_data)` "succeeds" without
issuing a single query, and the empty result gets cached. The first visible
error appears two chunks later, in `[tx2gene]`:

```
Error in `names(recipes) <- paste(names(recipes), "tx2gene", sep = "_")`:
! 'names' attribute [1] must be the same length as the vector [0]
```

If you see that, the cause is the working directory — not the recipes and not
biomaRt.

### The correct invocation

```bash
cd data-raw
rm -rf build_cache
Rscript -e 'rmarkdown::render("build.Rmd")'
```

`cd data-raw` is what makes `".."` resolve to the repo root, matching RStudio.

### Always clear the cache first

`opts_chunk$set(cache = TRUE)` means a previously failed or empty chunk is
replayed instead of re-run, so a real fix can look like it did nothing. Worse,
`[tx2gene]` **mutates** `recipes` in place (renaming entries and swapping
attributes), so a partially cached run feeds it already-consumed input.

```bash
rm -rf data-raw/build_cache
```

### Transient server errors (502/503/504/429)

`build.Rmd` handles these itself. A **504 Gateway Timeout** means the Ensembl
host took too long to answer — a server-side hiccup, not a configuration
problem, and typically momentary. Three mechanisms in the `functions` chunk:

- `with_retry()` retries transient HTTP failures (502, 503, 504, 429, timeouts)
  up to 3 times with exponential backoff. Permanent errors such as **405 are
  not retried** — they abort immediately rather than wasting a long run.
- `ENSEMBL_DELAY` (default 3s) spaces out requests. Archive hosts run on smaller
  infrastructure than the main site and degrade under back-to-back queries.
- `get_mart()` caches one mart per host+dataset. `useMart()` refetches the whole
  BioMart registry on every call, so the build would otherwise download the same
  registry document 22 times; caching cuts that to 11.

If a run still dies on 504 after retries, the host is genuinely unwell — wait
and rerun rather than lowering `ENSEMBL_DELAY`. Raising it is the safe direction.

### "The given dataset: X, is not valid" — usually an outage, not a bad recipe

This error is misleading. When the archive is briefly overloaded it serves an
HTML **"Service unavailable"** page with **HTTP 200**, and biomaRt parses that
page as the dataset list. The dataset is naturally absent from the HTML, so it
reports the name as invalid.

Observed directly: the dataset-list request returned 2,218 bytes of HTML at one
moment and the correct 25,964-byte TSV ten seconds later.

Because the HTTP status is 200, nothing looks like a network failure — which is
why `with_retry()` matches on `"is not valid"` and `"unavailable"` as well as on
status codes. Only suspect the recipe if the name is still rejected once the
host is reliably serving the real list:

```bash
curl -s "https://jun2026.archive.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL" \
  | grep drerio
```

An empty result plus a small response body means the host is down, not that the
dataset name is wrong.

### Request timeouts

biomaRt passes `getOption("timeout")` to `req_timeout()`, and R's default is
**60 seconds**. A cold archive host was measured taking **51.7s** just to return
the registry, leaving almost no headroom before the request aborts and surfaces
as a gateway timeout. `build.Rmd` sets `options(timeout = 600)`.

### Resuming an interrupted build

The archive fails often enough that a full 22-table build rarely completes on
the first attempt, so the build records its progress and resumes.

Each table is downloaded, saved to `data/` and logged one at a time, rather than
downloading everything and saving at the end. Completed table names are appended
to `data-raw/build_progress.log`, and a rerun skips anything listed there:

```
skipping grch38 (already built)
skipping grcm38 (already built)
drerio_gene_ensembl
```

So an interrupted run costs only the tables that had not finished. Just rerun
the same command — no flags, no manual bookkeeping.

The log is **deleted automatically** once every gene and tx2gene table is
present, leaving no state behind after a clean build. If it still exists when a
run ends, the build is incomplete and the final message names what is missing.

To force a full rebuild from scratch, delete it yourself:

```bash
rm -f data-raw/build_progress.log
```

The two download chunks set `cache=FALSE`, since this log supersedes knitr's
cache for resuming and a stale chunk cache would only mask real progress.

### Confirming the run actually did something

A successful rebuild refetches from Ensembl and rewrites the data files. Check:

```bash
git status --short data/
```

Modified `data/*.rda` means real queries ran. **An empty result means nothing
was fetched** — the run was a no-op regardless of what the console printed.

### What a full rebuild costs

22 BioMart queries (11 gene tables + 11 tx2gene), several MB each, taking many
minutes. It rewrites every `data/*.rda`, and gene counts and descriptions will
shift to whatever Ensembl release the recipe hosts point at. Treat it as a
deliberate data update, not a routine check.

Querying that fast can trip Ensembl's rate limiter; a mid-run **403** is
throttling, not a bug (see above). Wait rather than retrying.

### Generated files

`README.md` in the repo root is knitted from `README.Rmd` — edit the `.Rmd`.
A rebuild also leaves `data-raw/build.html`, which is a byproduct, not output.
