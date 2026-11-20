#!/usr/bin/env python3
"""Import CDS and existing annotations against the exact SynCom reference."""
import argparse
import csv
import re
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path


def fasta(path):
    with path.open("rb") as stream:
        name, chunks = None, []
        for line in stream:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks).upper()
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, b"".join(chunks).upper()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    component = args.source / "04.Component" / "Gene"
    genes_by_contig = defaultdict(list)
    reports = {}
    with args.contigs.open() as file:
        contigs = list(csv.DictReader(file, delimiter="\t"))
    by_genome = defaultdict(dict)
    for row in contigs:
        by_genome[row["genome_id"]][row["original_contig_id"]] = row["contig_id"]

    complement = bytes.maketrans(b"ACGTRYMKBDHVSWN", b"TGCAYRKMVHDBSWN")
    with ExitStack() as stack:
        def output(name, header):
            handle = stack.enter_context((args.output / name).open("w", newline=""))
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            if header:
                writer.writerow(header.split("\t"))
            return writer

        saf = output("genes.saf", "GeneID\tChr\tStart\tEnd\tStrand")
        bed = output("genes.bed", None)
        genes = output("genes.tsv", "gene_id\tgenome_id\tcontig_id\tstart\tend\tstrand\tlength_bp\tsource_gene_id")
        annotations = output("annotations.tsv", "gene_id\tsource\tfunction_id")
        report = output("import_report.tsv", "genome_id\tgff_cds\tverified\trejected\tannotation_links\tstatus")
        rejected = output("rejected_cds.tsv", "genome_id\tsource_gene_id\treason")
        sequences = {}
        for genome, names in by_genome.items():
            gff, ffn = (component / f"{genome}.{ext}" for ext in ("gff", "ffn"))
            if not gff.is_file() or not ffn.is_file():
                reports[genome] = [0, 0, 0, 0, "missing_gff_or_ffn"]
                continue
            reports[genome] = [0, 0, 0, 0, "ok"]
            cds = defaultdict(list)
            for line in gff.read_text().splitlines():
                fields = line.split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                raw = attrs.get("Parent") or attrs.get("ID")
                if raw:
                    cds[raw].append((fields[0], int(fields[3]), int(fields[4]), fields[6]))
                else:
                    reports[genome][0] += 1
                    reports[genome][2] += 1
                    rejected.writerow([genome, ".", "missing_gene_id"])
            for raw, segments in cds.items():
                reports[genome][0] += 1
                locations = {(piece[0], piece[3]) for piece in segments}
                if len(locations) != 1 or len(segments) != 1:
                    reports[genome][2] += 1
                    rejected.writerow([genome, raw, "multipart_cds"])
                    continue
                original, start, end, strand = segments[0]
                contig = names.get(original)
                if contig is None or strand not in {"+", "-"}:
                    reports[genome][2] += 1
                    rejected.writerow([genome, raw, "unmatched_contig_or_strand"])
                    continue
                genes_by_contig[contig].append((start, end, strand, genome, raw))
            sequences[genome] = ffn

        verified = defaultdict(set)
        loaded = {}
        remaining = {genome: set(names.values()) for genome, names in by_genome.items()}
        contig_genome = {row["contig_id"]: row["genome_id"] for row in contigs}
        expected_length = {row["contig_id"]: int(row["length_bp"]) for row in contigs}
        seen = set()
        for contig, sequence in fasta(args.reference):
            if contig in seen or len(sequence) != expected_length[contig]:
                raise ValueError(f"Reference disagrees with contig map: {contig}")
            seen.add(contig)
            genome = contig_genome[contig]
            if genome in sequences and genome not in loaded:
                loaded[genome] = dict(fasta(sequences[genome]))
            for start, end, strand, genome, raw in sorted(genes_by_contig.pop(contig, [])):
                fragment = sequence[start - 1:end]
                if strand == "-":
                    fragment = fragment.translate(complement)[::-1]
                if not (1 <= start <= end <= len(sequence)) or loaded[genome].get(raw) != fragment:
                    reports[genome][2] += 1
                    rejected.writerow([genome, raw, "cds_sequence_or_coordinate_mismatch"])
                    continue
                gene = f"{genome}:{raw}"
                saf.writerow([gene, contig, start, end, strand])
                bed.writerow([contig, start - 1, end, gene, genome])
                genes.writerow([gene, genome, contig, start, end, strand, len(fragment), raw])
                verified[genome].add(raw)
                reports[genome][1] += 1
            remaining[genome].remove(contig)
            if not remaining[genome]:
                loaded.pop(genome, None)
        if seen != set(contig_genome):
            raise ValueError("Reference is missing contigs")

        patterns = {
            "GOs": ("GO", r"GO:\d{7}"),
            "COG/KOG": ("COG", r"(?:COG|KOG)\d{4}"),
            "KO/Gene_ID": ("KO", r"\bK\d{5}\b"),
        }
        for genome, ids in verified.items():
            folder = args.source / "05.Annotation" / genome
            table = folder / "annotation.table.xls"
            links = set()
            if table.is_file():
                with table.open(errors="replace") as stream:
                    for row in csv.DictReader(stream, delimiter="\t"):
                        raw = row.get("#Qeury_name")
                        if raw not in ids:
                            continue
                        for column, (source, pattern) in patterns.items():
                            for term in re.findall(pattern, row.get(column) or ""):
                                links.add((raw, source, term))
            else:
                reports[genome][4] = "missing_annotation_table"
            for pattern, source, column in [
                ("CAZy/*.anno.xls", "CAZy", "Subject"),
                ("pfam/*.pfam.xls", "Pfam", "hmm-acc"),
            ]:
                for path in folder.glob(pattern):
                    with path.open(errors="replace") as stream:
                        reader = csv.reader(stream, delimiter="\t")
                        header = next(reader, [])
                        if column not in header:
                            continue
                        index = header.index(column)
                        for row in reader:
                            if len(row) > index and row[0] in ids and row[index]:
                                links.add((row[0], source, row[index]))
            for raw, source, term in sorted(links):
                annotations.writerow([f"{genome}:{raw}", source, term])
            reports[genome][3] = len(links)
        for genome, numbers in sorted(reports.items()):
            report.writerow([genome, *numbers])
    if not verified:
        raise ValueError("No CDS passed sequence verification")


if __name__ == "__main__":
    main()
