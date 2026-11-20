#!/usr/bin/env bash
set -Eeuo pipefail
root=$(dirname "$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")")
cd "$root"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/syncom smoke.XXXXXXXX")
trap 'rm -rf -- "$scratch"' EXIT

# Python only generates the synthetic FASTQ fixture.
python3 - "$scratch" <<'PY'
import gzip
import random
import sys
from pathlib import Path

root = Path(sys.argv[1])
randomizer = random.Random(1234)
reference = "".join(randomizer.choices("ACGT", k=12000))
(root / "tiny.fa").write_text(">reference_1\n" + reference + "\n")
(root / "map.tsv").write_text(
    "contig_id\tgenome_id\toriginal_contig_id\tlength_bp\n"
    "reference_1\tTEST01\tref\t12000\n"
)
(root / "manifest.tsv").write_text(
    "db_version\ttest\nreference\ttiny.fa\ncontig_map\tmap.tsv\n"
)
samples = ["sample_id\tr1\tr2\tlane_id"]
for sample, lane, offsets in [("s1", "L1", range(0, 20)), ("s1", "L2", range(20, 40)), ("s2", "L1", range(40))]:
    stem = f"{sample} {lane}"
    samples.append(f"{sample}\t{stem} R1.fq.gz\t{stem} R2.fq.gz\t{lane}")
    for mate in (1, 2):
        with gzip.open(root / f"{stem} R{mate}.fq.gz", "wt") as out:
            for offset in offsets:
                start = 100 * offset
                seq = reference[start : start + 150]
                if mate == 2:
                    seq = reference[start + 180 : start + 330].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                out.write(f"@{sample}-{offset}/{mate}\n{seq}\n+\n" + "I" * len(seq) + "\n")
(root / "samples.tsv").write_text("\n".join(samples) + "\n")
PY

tools=${1:-"$root/config/tools.tsv"}
tools=$(realpath -- "$tools")
configured() {
    local exe
    exe=$(awk -F '\t' -v name="$1" '$1==name { print $2 }' "$tools")
    if [[ $exe == */* && $exe != /* ]]; then exe=${tools%/*}/$exe; fi
    command -v -- "${exe:-$1}"
}
"$(configured bwa)" index "$scratch/tiny.fa"
"$(configured samtools)" faidx "$scratch/tiny.fa"
args=(--samples "$scratch/samples.tsv" --db "$scratch" --modules M03 --tools "$tools")
./syncom-metaflow run "${args[@]}" --outdir "$scratch/direct"
./syncom-metaflow plan "${args[@]}" --outdir "$scratch/exported"
(
    cd /
    nohup bash "$scratch/exported/scripts/run_all.sh" > "$scratch/exported/run.log" 2>&1 &
    wait "$!"
)
for file in qc.tsv strains.tsv strains_matrix.tsv; do
    cmp "$scratch/direct/summary/$file" "$scratch/exported/summary/$file"
done
awk -F '\t' 'NR>1 { if ($2!=40 || $3!=40) exit 1; n++ } END { if (n!=2) exit 1 }' "$scratch/direct/summary/qc.tsv"
awk -F '\t' 'NR>1 { if ($3!="TEST01" || $5!=80 || $9!=1) exit 1; n++ } END { if (n!=2) exit 1 }' "$scratch/direct/summary/strains.tsv"

bash "$scratch/exported/scripts/run_all.sh" > "$scratch/resume.log"
[[ $(grep -c ' cached$' "$scratch/resume.log") == 6 ]]
# Substitute a missing BAM index for an incomplete alignment from an interrupted run.
mv "$scratch/exported/work/s1/alignment/aligned.bam.bai" "$scratch/incomplete.bai"
bash "$scratch/exported/scripts/run_all.sh" > "$scratch/recovered.log"
grep -q '^s1 align_reads started' "$scratch/recovered.log"
grep -q '^s2 align_reads cached' "$scratch/recovered.log"
cmp "$scratch/direct/summary/strains.tsv" "$scratch/exported/summary/strains.tsv"
printf 'Smoke test passed.\n'
