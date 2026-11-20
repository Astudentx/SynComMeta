COMMANDS = r'''
# Concatenate decompressed lanes in manifest order, preserving mate pairing.
merge_lanes() {
    local dest=$1 lane
    shift
    for lane in "$@"; do
        case $lane in
            *.gz) gzip -dc -- "$lane" ;;
            *) cat -- "$lane" ;;
        esac
    done > "$dest"
}

M01() {
    local left=${r1[0]} right=${r2[0]}
    if ((${#r1[@]} > 1)); then
        left="$tmp/merged_R1.fq"
        right="$tmp/merged_R2.fq"
        merge_lanes "$left" "${r1[@]}"
        merge_lanes "$right" "${r2[@]}"
    fi
    "$fastp" --in1 "$left" --in2 "$right" --out1 "$tmp/R1.fq.gz" --out2 "$tmp/R2.fq.gz" \
        --thread "$threads" --dont_eval_duplication --json "$tmp/fastp.json" --html "$tmp/fastp.html"
}
'''

STEP = r'''
signature=$(printf '%s\n' "$method_sha" "$threads" "$fastp_version" "$(identity "$fastp")" "$reads_signature" | digest)
task M01 "$qc" "$signature" "$qc/R1.fq.gz" "$qc/R2.fq.gz" "$qc/fastp.json" "$qc/fastp.html"
r1=("$qc/R1.fq.gz")
r2=("$qc/R2.fq.gz")
'''
