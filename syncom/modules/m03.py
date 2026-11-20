COMMANDS = r'''
M03() {
    "$samtools" coverage -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
        "$bam" > "$tmp/contig_coverage.tsv"
    awk -F '\t' -v sample="$sample" -v db="$db_version" '
        BEGIN { OFS=FS; print "sample_id","db_version","genome_id","length_bp","reads_mapq20","covered_bp","breadth","mean_depth","target_depth_fraction","status" }
        FNR==NR {
            if (FNR == 1) next
            genome[$1]=$2; length_bp[$2]+=$4; contig_bp[$1]=$4
            next
        }
        FNR==1 { if ($1!="#rname" || $5!="covbases" || $7!="meandepth") exit 1; next }
        {
            if (!($1 in genome) || $3!=contig_bp[$1] || observed[$1]++) exit 1
            g=genome[$1]; reads[g]+=$4; covered[g]+=$5; bases[g]+=$7*contig_bp[$1]
        }
        END {
            for (contig in genome) if (!(contig in observed)) exit 1
            for (g in length_bp) { depth[g]=bases[g]/length_bp[g]; total+=depth[g] }
            for (g in length_bp) {
                state=(reads[g] ? "reference_evidence" : "no_mapq20_evidence")
                printf "%s\t%s\t%s\t%d\t%d\t%d\t%.8g\t%.8g\t%s\t%s\n", sample,db,g,length_bp[g],reads[g],covered[g],covered[g]/length_bp[g],depth[g],(total ? sprintf("%.8g",depth[g]/total) : "NA"),state
            }
        }
    ' "$contig_map" "$tmp/contig_coverage.tsv" |
        { IFS= read -r header; printf '%s\n' "$header"; sort -t $'\t' -k3,3; } > "$tmp/strains.tsv"
    awk -F '\t' -v sample="$sample" -v db="$db_version" '
        BEGIN { OFS=FS; print "sample_id","db_version","strain_group_id","reference_members","reads_mapq20","target_depth_fraction","status" }
        FILENAME==ARGV[1] {
            if (FNR>1) { group[$1]=$2; members[$2]=$3 }
            next
        }
        FNR>1 {
            g=group[$3]; reads[g]+=$5
            if($9!="NA") { total[g]+=$9; measured=1 }
        }
        END {
            for (g in members) {
                state=reads[g] ? "group_reference_evidence" : "no_mapq20_evidence"
                printf "%s\t%s\t%s\t%d\t%d\t%s\t%s\n",sample,db,g,members[g],reads[g],(measured ? sprintf("%.8g",total[g]) : "NA"),state
            }
        }
    ' "$groups" "$tmp/strains.tsv" |
        { IFS= read -r header; printf '%s\n' "$header"; sort -t $'\t' -k3,3; } > "$tmp/strain_groups.tsv"
}
'''

STEP = r'''
strains="$out/strains/$sample"
signature=$(printf '%s\n' "$bam_signature" "$metadata_signature" | digest)
task M03 "$strains" "$signature" "$strains/contig_coverage.tsv" "$strains/strains.tsv" "$strains/strain_groups.tsv"
'''
