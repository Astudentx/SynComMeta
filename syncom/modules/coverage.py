COMMANDS = r'''
prepare_reference_groups() {
    local metadata="$plan/strains.tsv"
    [[ -f $metadata ]] || metadata=/dev/null
    awk -F '\t' '
        BEGIN { OFS=FS; print "genome_id","strain_group_id","reference_members" }
        FILENAME==ARGV[1] {
            if (FNR==1) {
                for (i=1;i<=NF;i++) {
                    if ($i=="ID") idcol=i
                    if ($i=="ANI >= 99.99 Group") anicol=i
                }
            } else if (idcol && anicol && $anicol ~ /^[1-9][0-9]*$/) {
                group[$idcol]=$anicol; size[$anicol]++
            }
            next
        }
        FNR>1 && !seen[$2]++ {
            g=$2
            if (g in group && size[group[g]]>1) print g,"ANI99.99:" group[g],size[group[g]]
            else print g,"genome:" g,1
        }
    ' "$metadata" "$contig_map" > "$groups"
}

eligible_genomes() {
    awk -F '\t' -v min_depth="$min_genome_depth" -v min_breadth="$min_genome_breadth" '
        BEGIN { OFS=FS }
        FILENAME==ARGV[1] { if(FNR>1) members[$1]=$3; next }
        FNR==1 { print "sample_id","db_version","genome_id","breadth","mean_depth","status"; next }
        {
            g=$3
            if ($7<min_breadth || $8<min_depth) state="low_coverage"
            else if (members[g]>1) state="near_identical_reference_group"
            else state="eligible"
            print $1,$2,g,$7,$8,state
        }
    ' "$groups" "$strains/strains.tsv" > "$tmp/eligibility.tsv"
}
'''

STEP = r'''
gate_signature=$(printf '%s\n' "$bam_signature" "$metadata_signature" "$min_genome_depth" "$min_genome_breadth" "$(identity "$strains/strains.tsv")" | digest)
'''
