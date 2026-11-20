COMMANDS = r'''
M04() {
    local library_size help
    local -a pairing=(-p)
    help=$("$featureCounts" -h 2>&1 || true)
    if [[ $help == *--countReadPairs* ]]; then pairing+=(--countReadPairs); fi

    # Remove ambiguous, non-primary and flagged records before fragment counting.
    "$samtools" view -c -f 64 -F 2304 "$bam" > "$tmp/clean_fragments.txt"
    read -r library_size < "$tmp/clean_fragments.txt"
    "$samtools" view -b -f 3 -F 3844 -q 20 -o "$tmp/countable.bam" "$bam"
    "$featureCounts" -T "$((threads > 64 ? 64 : threads))" "${pairing[@]}" -B -C \
        -F SAF -a "$genes_saf" -o "$tmp/featurecounts.tsv" "$tmp/countable.bam"
    awk -F '\t' -v s="$sample" -v db="$db_version" -v n="$library_size" '
        BEGIN { OFS=FS; print "sample_id","db_version","gene_id","length_bp","fragment_count","fpkm" }
        /^#/ || $1=="Geneid" { next }
        NF!=7 || $6<=0 || $7 !~ /^[0-9]+$/ || seen[$1]++ { exit 1 }
        { print s,db,$1,$6,$7,(n ? sprintf("%.9g",1e9*$7/($6*n)) : "NA"); rows++ }
        END { if (!rows) exit 1 }
    ' "$tmp/featurecounts.tsv" > "$tmp/genes.tsv"
    awk -F '\t' -v s="$sample" -v db="$db_version" '
        BEGIN { OFS=FS; print "sample_id","db_version","function_key","source","function_id","fragment_count","fpkm_sum","reference_genes" }
        FNR==NR {
            if(FNR>1) { counts[$3]=$5; fpkm[$3]=$6 }
            next
        }
        FNR==1 { next }
        {
            if (!($1 in counts)) exit 1
            key=$2 ":" $3
            if ($1!=last_gene) { delete seen; last_gene=$1 }
            if (seen[key]++) next
            source[key]=$2; term[key]=$3; total[key]+=counts[$1]; genes[key]++
            if (fpkm[$1]=="NA") missing[key]=1
            else normalized[key]+=fpkm[$1]
        }
        END {
            for(key in genes) print s,db,key,source[key],term[key],total[key],(missing[key] ? "NA" : sprintf("%.9g",normalized[key])),genes[key]
        }
    ' "$tmp/genes.tsv" "$annotations" |
        { IFS= read -r header; printf '%s\n' "$header"; sort -t $'\t' -k3,3; } > "$tmp/functions.tsv"
}
'''

STEP = r'''
functions="$out/functions/$sample"
signature=$(printf '%s\n' "$bam_signature" "$gene_signature" "$featurecounts_version" "$(identity "$featureCounts")" | digest)
task M04 "$functions" "$signature" "$functions/genes.tsv" "$functions/functions.tsv" "$functions/featurecounts.tsv" "$functions/featurecounts.tsv.summary"
'''
