COMMANDS = r'''
M06() {
    eligible_genomes
    "$samtools" bedcov -Q 20 -G SUPPLEMENTARY -j -d 1 "$genes_bed" "$bam" > "$tmp/bedcov.tsv"
    awk -F '\t' -v sample="$sample" -v db="$db_version" '
        BEGIN { OFS=FS; print "sample_id","db_version","gene_id","genome_id","breadth","mean_depth","genome_mean_depth","pav","relative_copy_number","status" }
        FNR==NR { if(FNR>1) {state[$3]=$6; baseline[$3]=$5} next }
        NF!=7 || $3<=$2 || !($5 in state) { exit 1 }
        {
            length_bp=$3-$2; breadth=$7/length_bp; depth=$6/length_bp
            g=$5; quality=state[g]; pav="uncertain"; cn="NA"
            if(quality=="eligible" && breadth>=0.8 && depth>=1) {
                pav="present"; cn=sprintf("%.8g",depth/baseline[g])
            } else if(quality=="eligible" && breadth==0 && baseline[g]>=10) {
                pav="absent"; cn=0
            } else if(quality=="eligible") {
                quality="insufficient_gene_evidence"
            }
            printf "%s\t%s\t%s\t%s\t%.8g\t%.8g\t%s\t%s\t%s\t%s\n",sample,db,$4,g,breadth,depth,baseline[g],pav,cn,quality
            rows++
        }
        END {if(!rows) exit 1}
    ' "$tmp/eligibility.tsv" "$tmp/bedcov.tsv" > "$tmp/pav_cnv.tsv"
}
'''

STEP = r'''
pav="$out/pav_cnv/$sample"
signature=$(printf '%s\n' "$gate_signature" "$bed_signature" | digest)
task M06 "$pav" "$signature" "$pav/pav_cnv.tsv" "$pav/eligibility.tsv"
'''
