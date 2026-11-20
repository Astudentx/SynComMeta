#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
source "$(dirname -- "${BASH_SOURCE[0]}")/tables.bash"

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
[[ $# -ge 4 && $1 == --outdir ]] || fail 'Usage: syncom-metaflow merge --outdir NEW_DIR BATCH_DIR BATCH_DIR [...]'
out=$(realpath -m -- "$2")
shift 2

declare -a batches=()
declare -A sample_batch=()
expected_config=
for path in "$@"; do
    batch=$(realpath -e -- "$path")
    [[ -f $batch/.syncom.lock ]] || fail "Not a run directory: $batch"
    exec {fd}<"$batch/.syncom.lock"
    flock -sn "$fd" || fail "Batch is busy: $batch"
    IFS=$'\t' read -r state plan < "$batch/provenance/status.tsv"
    [[ $state == complete ]] || fail "Batch is incomplete: $batch"
    [[ -f $batch/provenance/summary_catalog.tsv && -f $plan/reads.tsv ]] ||
        fail "Batch is missing its summary catalog or sample plan: $batch"
    fingerprint=$(<"$batch/provenance/compatibility.sha256")
    [[ $fingerprint =~ ^[a-f0-9]{64}$ ]] || fail "Invalid compatibility fingerprint: $batch"
    if [[ -z $expected_config ]]; then
        expected_config=$fingerprint
        catalog="$batch/provenance/summary_catalog.tsv"
    fi
    [[ $fingerprint == "$expected_config" ]] && cmp -s "$catalog" "$batch/provenance/summary_catalog.tsv" ||
        fail "Incompatible database, settings, tools, or module selection: $batch"
    mapfile -t batch_ids < <(cut -f1 "$plan/reads.tsv" | sort -u)
    for id in "${batch_ids[@]}"; do
        [[ ! ${sample_batch[$id]+present} ]] || fail "Duplicate sample ID across batches: $id"
        sample_batch[$id]=$batch
    done
    batches+=("$batch")
done
mapfile -t ids < <(printf '%s\n' "${!sample_batch[@]}" | sort)

mkdir -- "$out" || fail "Merge needs a new output directory: $out"
mkdir "$out"/{summary,provenance}
tmp=
trap 'if [[ -n $tmp ]]; then rm -rf -- "$tmp"; fi' EXIT
trap 'printf "failed\n" > "$out/provenance/status.tsv"' ERR
printf 'running\n' > "$out/provenance/status.tsv"
cp -- "$catalog" "$out/provenance/summary_catalog.tsv"
printf '%s\n' "$expected_config" > "$out/provenance/compatibility.sha256"
{
    printf 'sample_id\tbatch_dir\n'
    for id in "${ids[@]}"; do printf '%s\t%s\n' "$id" "${sample_batch[$id]}"; done
} > "$out/provenance/samples.tsv"
while IFS=$'\t' read -r dir src file key value dest second_value second_dest; do
    tables=()
    for batch in "${batches[@]}"; do tables+=("$batch/summary/$file"); done
    awk 'FNR==1 { if (NR==1) first=$0; else if ($0!=first) exit 1 }
        FNR>1 || NR==1' "${tables[@]}" > "$out/summary/$file"
    catalog_matrices "$out/summary" "$file" "$key" "$value" "$dest" "$second_value" "$second_dest"
done < "$catalog"
printf 'complete\n' > "$out/provenance/status.tsv"
printf 'Merged %s sample(s) from %s batches.\n' "${#ids[@]}" "${#batches[@]}"
