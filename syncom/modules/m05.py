COMMANDS = r'''
M05() {
    eligible_genomes
    awk -F '\t' 'BEGIN {OFS=FS} FILENAME==ARGV[1] {if ($6=="eligible") ok[$3]=1; next}
        FNR>1 && ($2 in ok) {print $1,0,$4}' "$tmp/eligibility.tsv" "$contig_map" > "$tmp/eligible.bed"
    local header
    header=$'sample_id\tdb_version\tgenome_id\tcontig_id\tposition_1based\tref\talt\tdepth\talt_depth\talt_fraction\talt_forward\talt_reverse'
    printf '%s\n' "$header" > "$tmp/snv.tsv"
    printf '#contig\tstart0\tend\n' > "$tmp/callable.bed"
    if [[ -s $tmp/eligible.bed ]]; then
        "$bcftools" mpileup --ignore-RG -I -q 20 -Q 20 -d 10000 --ns UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
            -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f "$reference" -R "$tmp/eligible.bed" -Ou "$bam" |
        "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%ADF\t%ADR]\n' |
        awk -F '\t' -v sample="$sample" -v db="$db_version" -v callable="$tmp/callable.bed" '
            BEGIN { OFS=FS }
            FILENAME==ARGV[1] { if(FNR>1) genome[$1]=$2; next }
            $3 ~ /^[ACGT]$/ {
                n=split($4, alts, ",")
                split($5, count, ","); split($6, fw, ","); split($7, rv, ",")
                depth=count[1]+0
                for(i=1;i<=n;i++) if(alts[i] ~ /^[ACGT]$/) depth+=count[i+1]
                if(depth<5) next
                if($1!=chrom || $2!=last+1) {
                    if(chrom!="") print chrom,start,last >> callable
                    chrom=$1; start=$2-1
                }
                last=$2
                for(i=1;i<=n;i++) {
                    alt=alts[i]; j=i+1
                    if(alt ~ /^[ACGT]$/ && alt!=$3 && count[j]>=3 && fw[j]>=1 && rv[j]>=1 && count[j]/depth>=0.05)
                        print sample,db,genome[$1],$1,$2,$3,alt,depth,count[j],sprintf("%.8g",count[j]/depth),fw[j],rv[j]
                }
            }
            END { if(chrom!="") print chrom,start,last >> callable }
        ' "$contig_map" - >> "$tmp/snv.tsv"
    fi
}
'''

STEP = r'''
snv="$out/snv/$sample"
signature=$(printf '%s\n' "$gate_signature" "$bcftools_version" "$(identity "$bcftools")" | digest)
task M05 "$snv" "$signature" "$snv/snv.tsv" "$snv/callable.bed" "$snv/eligibility.tsv"
'''
