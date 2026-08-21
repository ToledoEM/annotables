#!/usr/bin/env bash
#
# Switch the Ensembl BioMart `host:` field across all recipe YAMLs.
#
# Ensembl's www host currently 308-redirects BioMart to a GET-only archive and
# answers 405 to biomaRt's POST, so the recipes are pinned to an archive host.
# Use this script to move off that pin (once www is healthy again) or to re-pin
# to a different Ensembl release.
#
# See data-raw/recipes/README.md for the full diagnosis.
#
# Usage:
#   switch-host.sh [--dry-run] <host-url|archive-tag>
#
# Examples:
#   switch-host.sh jun2026                     # -> https://jun2026.archive.ensembl.org
#   switch-host.sh https://www.ensembl.org     # back to the normal host
#   switch-host.sh --dry-run https://www.ensembl.org
#
# grch37.yml is ALWAYS skipped: https://grch37.ensembl.org is a separate,
# permanently maintained archive that works and must not be switched.

set -euo pipefail

readonly SKIP_FILE="grch37.yml"

usage() {
    sed -n '3,21p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-1}"
}

dry_run=false
target=""

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run|-n) dry_run=true; shift ;;
        -h|--help)    usage 0 ;;
        -*)           echo "error: unknown option: $1" >&2; usage ;;
        *)
            [ -n "$target" ] && { echo "error: multiple hosts given: '$target' and '$1'" >&2; usage; }
            target="$1"; shift ;;
    esac
done

[ -z "$target" ] && { echo "error: no host given" >&2; usage; }

command -v yq >/dev/null 2>&1 || {
    echo "error: yq is required (brew install yq)" >&2
    exit 1
}

# Bare archive tag (e.g. "jun2026") expands to the full archive URL.
if [[ "$target" != http://* && "$target" != https://* ]]; then
    if [[ "$target" =~ ^[a-z]{3}[0-9]{4}$ ]]; then
        target="https://${target}.archive.ensembl.org"
    else
        echo "error: '$target' is neither a URL nor an archive tag like 'jun2026'" >&2
        exit 1
    fi
fi

target="${target%/}"

recipes_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/recipes"
[ -d "$recipes_dir" ] || { echo "error: recipes dir not found: $recipes_dir" >&2; exit 1; }

if [ "$dry_run" = true ]; then
    echo "DRY RUN — no files written"
fi
echo "target host: $target"
echo

changed=0
for f in "$recipes_dir"/*.yml; do
    name="$(basename "$f")"

    if [ "$name" = "$SKIP_FILE" ]; then
        printf '  %-20s SKIPPED (separate live archive)\n' "$name"
        continue
    fi

    current="$(yq '.host' "$f")"

    if [ "$current" = "$target" ]; then
        printf '  %-20s unchanged  %s\n' "$name" "$current"
        continue
    fi

    printf '  %-20s %s -> %s\n' "$name" "$current" "$target"
    if [ "$dry_run" = false ]; then
        yq -i ".host = \"$target\"" "$f"
    fi
    changed=$((changed + 1))
done

echo
if [ "$dry_run" = true ]; then
    echo "$changed file(s) would change. Re-run without --dry-run to apply."
    exit 0
fi

echo "$changed file(s) updated."
echo
echo "No R changes needed: get_data() in data-raw/build.Rmd re-pins mart@host"
echo "from recipe\$host, so this script is the whole switch."
echo "If rebuilding, clear the knitr cache first: rm -rf data-raw/build_cache"
