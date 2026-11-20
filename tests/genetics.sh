#!/usr/bin/env bash
set -Eeuo pipefail
root=$(dirname "$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")")
cd "$root"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/syncom genetics.XXXXXXXX")
trap 'rm -rf -- "$scratch"' EXIT

python3 - "$scratch" <<'PY'
import gzip
import random
import sys
from pathlib import Path

root = Path(sys.argv[1])
reference = "".join(random.Random(1256).choices("ACGT", k=12000))
alt = next(base for base in "ACGT" if base != reference[1700])
mutant = reference[:1700] + alt + reference[1701:]
(root / "expected.tsv").write_text(f"1701\t{reference[1700]}\t{alt}\n")
(root / "ref.fa").write_text(">reference_1\n" + reference + "\n")
(root / "map.tsv").write_text(
    "contig_id\tgenome_id\toriginal_contig_id\tlength_bp\nreference_1\tTEST01\tref\t12000\n"
)
(root / "genes.saf").write_text(
    "GeneID\tChr\tStart\tEnd\tStrand\nTEST01:g1\treference_1\t1001\t2600\t+\n"
    "TEST01:g2\treference_1\t11001\t11500\t+\n"
)
(root / "genes.bed").write_text(
    "reference_1\t1000\t2600\tTEST01:g1\tTEST01\nreference_1\t11000\t11500\tTEST01:g2\tTEST01\n"
)
(root / "annotations.tsv").write_text(
    "gene_id\tsource\tfunction_id\nTEST01:g1\tKO\tK00001\nTEST01:g2\tKO\tK00002\n"
)
(root / "manifest.tsv").write_text(
    "db_version\ttest\nreference\tref.fa\ncontig_map\tmap.tsv\n"
    "genes_saf\tgenes.saf\ngenes_bed\tgenes.bed\nannotations\tannotations.tsv\n"
)
header = "sample_id\tr1\tr2\n"
samples = []
for sample, step in [("high", 12), ("low", 200)]:
    row = f"{sample}\t{sample} R1.fq.gz\t{sample} R2.fq.gz\n"
    samples.append(row)
    (root / f"{sample}.tsv").write_text(header + row)
    for mate in (1, 2):
        with gzip.open(root / f"{sample} R{mate}.fq.gz", "wt") as stream:
            for i, start in enumerate(range(0, 10000, step)):
                seq = mutant[start : start + 150]
                if mate == 2:
                    seq = mutant[start + 180 : start + 330].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                stream.write(f"@{sample}-{i}/{mate}\n{seq}\n+\n" + "I" * len(seq) + "\n")
(root / "combined.tsv").write_text(header + "".join(samples))
PY

tools=$(realpath -- "${1:-"$root/config/tools.tsv"}")
configured() {
    local path
    path=$(awk -F '\t' -v name="$1" '$1==name { print $2 }' "$tools")
    if [[ $path == */* && $path != /* ]]; then path=${tools%/*}/$path; fi
    command -v -- "${path:-$1}"
}
"$(configured bwa)" index "$scratch/ref.fa"
"$(configured samtools)" faidx "$scratch/ref.fa"
for run in combined high low; do
    ./syncom-metaflow run --db "$scratch" --tools "$tools" --samples "$scratch/$run.tsv" \
        --modules M04,M05,M06 --outdir "$scratch/$run"
done
summary=$scratch/combined/summary
tail -n +2 "$summary/snv.tsv" | cut -f5-7 | diff - "$scratch/expected.tsv"
awk -F '\t' '
    $1=="high" && $3=="TEST01:g1" && $5==161 { positive=1 }
    $1=="high" && $3=="TEST01:g2" && $5==0 { negative=1 }
    END { exit !(positive && negative) }
' "$summary/genes.tsv"
awk -F '\t' '
    $1=="high" && $3=="TEST01:g1" && $8=="present" && $9>1 { present=1 }
    $1=="high" && $3=="TEST01:g2" && $8=="absent" && $9==0 { absent=1 }
    $1=="low" && $8=="uncertain" && $9=="NA" && $10=="low_coverage" { uncertain++ }
    END { exit !(present && absent && uncertain==2) }
' "$summary/pav_cnv.tsv"
awk -F '\t' 'NR>1 && $2<1701 && $3>=1701 { covered=1 } END { exit !covered }' "$scratch/combined/snv/high/callable.bed"
./syncom-metaflow merge --outdir "$scratch/merged" "$scratch/low" "$scratch/high"
for matrix in gene_fpkm function_fpkm gene_pav gene_cnv strains; do
    cmp "$summary/${matrix}_matrix.tsv" "$scratch/merged/summary/${matrix}_matrix.tsv"
done
if ./syncom-metaflow merge --outdir "$scratch/duplicated" "$scratch/high" "$scratch/high"; then
    printf 'Duplicate sample IDs were accepted.\n' >&2
    exit 1
fi
printf 'Genetics and merge test passed.\n'
