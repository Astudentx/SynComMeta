summarize_qc() {
    local id
    tmp=$(mktemp -d "$out/summary/.summary.XXXXXXXX")
    {
        printf 'sample_id\tinput_pairs\tclean_pairs\n'
        for id in "${ids[@]}"; do
            awk -v id="$id" '
                /"before_filtering":/ { section="before" }
                /"after_filtering":/ { section="after" }
                /"total_reads":/ && section {
                    sub(/^.*:[[:space:]]*/, ""); sub(/,.*$/, "")
                    if ($0 !~ /^[0-9]+$/ || $0%2) exit 1
                    pairs[section]=$0/2
                    section=""
                }
                END {
                    if (!("before" in pairs) || !("after" in pairs)) exit 1
                    printf "%s\t%.0f\t%.0f\n",id,pairs["before"],pairs["after"]
                }
            ' "$out/qc/$id/fastp.json"
        done
    } > "$tmp/qc.tsv"
    mv -f -- "$tmp/qc.tsv" "$out/summary/qc.tsv"
    rmdir -- "$tmp"
    tmp=
}

collect_tables() {
    local module=$1 file=$2 name=${3:-$2} id
    local -a tables=()
    for id in "${ids[@]}"; do tables+=("$out/$module/$id/$file"); done
    tmp=$(mktemp -d "$out/summary/.summary.XXXXXXXX")
    awk 'FNR>1 || NR==1' "${tables[@]}" > "$tmp/$name"
    mv -f -- "$tmp/$name" "$out/summary/$name"
    rmdir -- "$tmp"
    tmp=
}

catalog_matrices() {
    local directory=$1 file=$2 key=$3 value=$4 dest=$5 second_value=$6 second_dest=$7
    if [[ -n $dest ]]; then matrix "$directory/$file" "$key" "$value" "$directory/$dest"; fi
    if [[ -n $second_dest ]]; then matrix "$directory/$file" "$key" "$second_value" "$directory/$second_dest"; fi
}

matrix() {
    local input=$1 entity_col=$2 value_col=$3 dest=$4
    tmp=$(mktemp -d "$out/summary/.summary.XXXXXXXX")
    tail -n +2 -- "$input" | sort -t $'\t' -k"${entity_col},${entity_col}" -k1,1 |
    awk -F '\t' -v entity="$entity_col" -v value="$value_col" -v samples="$(printf '%s ' "${ids[@]}")" '
        function emit(   i) {
            printf "%s",key
            for(i=1;i<=n;i++) printf "\t%s",(ids[i] in abundance ? abundance[ids[i]] : "NA")
            print ""
        }
        BEGIN { n=split(samples, ids, " "); printf "%s","entity_id"; for(i=1;i<=n;i++) printf "\t%s",ids[i]; print "" }
        {
            if (key!=$entity) {
                if (key!="") emit()
                delete abundance
                key=$entity
            }
            if ($1 in abundance) exit 1
            abundance[$1]=$value
        }
        END { if (key!="") emit() }
    ' > "$tmp/matrix.tsv"
    mv -f -- "$tmp/matrix.tsv" "$dest"
    rmdir -- "$tmp"
    tmp=
}
