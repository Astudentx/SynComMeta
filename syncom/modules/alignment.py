COMMANDS = r'''
align_reads() {
    # BWA and the sorter overlap. At least one core is needed for each process.
    local align_threads=$((threads > 1 ? threads - 1 : 1))
    "$bwa" mem -t "$align_threads" "$reference" "${r1[0]}" "${r2[0]}" |
        "$samtools" sort -m 512M -o "$tmp/aligned.bam" -
    "$samtools" quickcheck "$tmp/aligned.bam"
    "$samtools" index "$tmp/aligned.bam"
}

'''

STEP = r'''
align="$out/work/$sample/alignment"
bam="$align/aligned.bam"
signature=$(printf '%s\n' "$method_sha" "$threads" "$db_version" "$reads_signature" "$ref_signature" "$bwa_version" "$samtools_version" "$(identity "$bwa" "$samtools")" | digest)
task align_reads "$align" "$signature" "$bam" "$bam.bai"
bam_signature=$(printf '%s\n' "$signature" "$(identity "$bam" "$bam.bai")" | digest)
'''
