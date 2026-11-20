COMMANDS = r'''
M02() {
    local left=reads_R1.fq right=reads_R2.fq
    [[ ${r1[0]} != *.gz ]] || left+=.gz
    [[ ${r2[0]} != *.gz ]] || right+=.gz
    ln -s -- "${r1[0]}" "$tmp/$left"
    ln -s -- "${r2[0]}" "$tmp/$right"
    ln -s -- "$motus_db" "$tmp/db_mOTU"
    (
        cd -- "$tmp"
        "$motus" profile -db db_mOTU -f "$left" -r "$right" -t "$threads" \
            -n "$sample" -o motus.tsv
    )
    awk -F '\t' -v sample="$sample" -v db="$db_version" '
        BEGIN { OFS=FS; print "sample_id","db_version","motu_id","label","relative_abundance" }
        /^#/ || NF==0 { next }
        {
            if (NF!=2 || $2 !~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/) exit 1
            id=$1
            if (id ~ /\[[^]]+\]$/) { sub(/^.*\[/,"",id); sub(/\]$/,"",id) }
            if (id in seen) {
                # The local extended database emits the same unassigned value twice.
                if (id=="unassigned" && $2==seen[id]) next
                print "Unexpected duplicate mOTU: " id > "/dev/stderr"
                exit 1
            }
            seen[id]=$2
            print sample,db,id,$1,$2
            rows++
        }
        END { if (!rows) exit 1 }
    ' "$tmp/motus.tsv" > "$tmp/taxonomy.tsv"
}
'''

STEP = r'''
tax="$out/taxonomy/$sample"
signature=$(printf '%s\n' "$method_sha" "$threads" "$db_version" "$reads_signature" "$motus_signature" "$motus_version" "$bwa_version" "$(identity "$motus" "$bwa")" | digest)
task M02 "$tax" "$signature" "$tax/motus.tsv" "$tax/taxonomy.tsv"
'''
