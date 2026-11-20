{
    if ((M01)); then printf 'qc	qc.tsv	qc.tsv\n'; fi
    if ((M02)); then printf 'taxonomy	taxonomy.tsv	taxonomy.tsv	3	5	taxonomy_matrix.tsv\n'; fi
    if ((M03)); then
        printf 'strains	strains.tsv	strains.tsv	3	9	strains_matrix.tsv\n'
        printf 'strains	strain_groups.tsv	strain_groups.tsv	3	6	strain_groups_matrix.tsv\n'
    fi
    if ((M04)); then
        printf 'functions	genes.tsv	genes.tsv	3	6	gene_fpkm_matrix.tsv\n'
        printf 'functions	functions.tsv	functions.tsv	3	7	function_fpkm_matrix.tsv\n'
    fi
    if ((M05)); then
        printf 'snv	snv.tsv	snv.tsv\n'
        printf 'snv	eligibility.tsv	snv_eligibility.tsv\n'
    fi
    if ((M06)); then
        printf 'pav_cnv	pav_cnv.tsv	pav_cnv.tsv	3	8	gene_pav_matrix.tsv	9	gene_cnv_matrix.tsv\n'
        printf 'pav_cnv	eligibility.tsv	pav_cnv_eligibility.tsv\n'
    fi
} > "$out/provenance/summary_catalog.tsv"
while IFS=$'\t' read -r dir src file key value dest second_value second_dest; do
    if [[ $dir == qc ]]; then summarize_qc; else collect_tables "$dir" "$src" "$file"; fi
    catalog_matrices "$out/summary" "$file" "$key" "$value" "$dest" "$second_value" "$second_dest"
done < "$out/provenance/summary_catalog.tsv"
printf 'complete\t%s\n' "$plan" > "$out/provenance/status.tsv"
printf 'Completed %s sample(s).\n' "${#ids[@]}"
