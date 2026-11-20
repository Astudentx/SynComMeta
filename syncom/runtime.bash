#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

plan=$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")
source "$plan/settings.sh"
if [[ ${1:-} != --lock-held || ! /proc/$$/fd/9 -ef $out/.syncom.lock ]]; then
    exec 9>"$out/.syncom.lock"
fi
flock -n 9 || { printf 'Output directory is busy: %s\n' "$out" >&2; exit 1; }
[[ $out/shell -ef $plan ]] || { printf 'This plan has been superseded.\n' >&2; exit 1; }

source "$plan/commands.sh"
mkdir -p "$out"/{logs,work,provenance,summary}
declare -a ids=()
mapfile -t ids < <(cut -f1 "$plan/reads.tsv" | sort -u)
method_sha=$(cat "$plan/method.sha256")
if ((M02)); then export PATH="${bwa%/*}:$PATH"; fi
tmp=
trap 'if [[ -n $tmp ]]; then rm -rf -- "$tmp"; fi' EXIT
trap 'printf "interrupted\t%s\n" "$plan" > "$out/provenance/status.tsv"; exit 130' INT
trap 'printf "interrupted\t%s\n" "$plan" > "$out/provenance/status.tsv"; exit 143' TERM
trap 'printf "failed\t%s\n" "$plan" > "$out/provenance/status.tsv"; printf "Run failed. See %s/logs\n" "$out" >&2' ERR
printf 'running\t%s\n' "$plan" > "$out/provenance/status.tsv"

# The immutable plan plus file identities are used to invalidate affected tasks.
identity() { stat -Lc '%n|%d:%i|%s|%y|%z' -- "$@"; }
digest() { sha256sum | cut -d ' ' -f1; }
task() {
    local name=$1 dir=$2 fingerprint=$3
    shift 3
    local stamp="$dir/.complete" log="$out/logs/$sample.$name.log"
    local product current
    mkdir -p "$dir"
    if [[ -f $stamp ]] && current=$(identity "$@" 2>/dev/null) &&
        [[ $(cat -- "$stamp") == "$fingerprint"$'\n'"$current" ]]; then
        printf '%s %s cached\n' "$sample" "$name"
        return
    fi
    rm -f -- "$stamp"
    tmp=$(mktemp -d "$dir/.task.XXXXXXXX")
    printf '%s %s started (%s)\n' "$sample" "$name" "$log"
    # A simple command preserves errexit inside the module function.
    { set -x; "$name"; set +x; } >"$log" 2>&1
    for product in "$@"; do
        [[ -s "$tmp/${product##*/}" ]]
        mv -f -- "$tmp/${product##*/}" "$product"
    done
    { printf '%s\n' "$fingerprint"; identity "$@"; } > "$tmp/completed"
    mv -f -- "$tmp/completed" "$stamp"
    rm -rf -- "$tmp"
    tmp=
}

if ((M01)); then fastp_version=$("$fastp" --version 2>&1); fi
if ((M02 || M03 || M04)); then bwa_version=$("$bwa" 2>&1 || true); fi
if ((M03 || M04)); then
    samtools_version=$("$samtools" --version)
    ref_signature=$(identity "$reference" "$reference".{amb,ann,bwt,pac,sa,fai} "$contig_map")
fi
if ((M02)); then
    motus_version=$("$motus" profile -db "$motus_db" --version 2>&1)
    motus_signature=$(identity "$motus_db"/db_mOTU_*)
fi
if ((M04)); then
    featurecounts_version=$("$featureCounts" -v 2>&1)
    gene_signature=$(identity "$genes_saf" "$annotations")
fi
if ((M05)); then bcftools_version=$("$bcftools" --version); fi
if ((M06)); then bed_signature=$(identity "$genes_bed"); fi
metadata_signature=none
if [[ -f $plan/strains.tsv ]]; then metadata_signature=$(digest < "$plan/strains.tsv"); fi
{
    printf 'db_version\t%s\nthreads\t%s\n' "$db_version" "$threads"
    if ((M01)); then printf 'fastp\t%s\n' "${fastp_version%%$'\n'*}"; fi
    if ((M02)); then printf 'motus\t%s\n' "${motus_version%%$'\n'*}"; fi
    if ((M02 || M03 || M04)); then printf '%s\n' "$bwa_version" | awk '/^Version:/ { print "bwa\t" $2 }'; fi
    if ((M03 || M04)); then printf 'samtools\t%s\n' "${samtools_version%%$'\n'*}"; fi
    if ((M04)); then printf '%s\n' "$featurecounts_version" | awk '/featureCounts/ { print "featureCounts\t" $2 }'; fi
    if ((M05)); then printf 'bcftools\t%s\n' "${bcftools_version%%$'\n'*}"; fi
} > "$out/provenance/versions.tsv"
{
    printf '%s\n' "$method_sha" "$clean" "$M01,$M02,$M03,$M04,$M05,$M06"
    printf '%s\n' "${ref_signature:-}" "${motus_signature:-}" "${gene_signature:-}" "${bed_signature:-}" "$metadata_signature"
    if ((M05 || M06)); then printf '%s\n' "$min_genome_depth" "$min_genome_breadth"; fi
    cat "$out/provenance/versions.tsv"
} | digest > "$out/provenance/compatibility.sha256"
cp -- "$plan/samples.input.tsv" "$out/provenance/samples.tsv"
cp -- "$plan/database.tsv" "$out/provenance/database.tsv"
if [[ -f $plan/strains.tsv ]]; then cp -- "$plan/strains.tsv" "$out/provenance/strains.tsv"; fi

if ((M03)); then
    groups="$out/provenance/reference_groups.tsv"
    prepare_reference_groups
fi

