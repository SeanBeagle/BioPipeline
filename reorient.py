from Bio import SeqIO
from Bio import Align
import pandas as pd


def get_kmers(string, k):
    kmers = set()
    for i in range(0, len(string) - k, k):
        kmers.add(string[i:i + k])
    return kmers



def pairwise_reorient(reference, fasta_in, fasta_out, csv):
    print("Reorienting using global pairwise alignment")
    aligner = Align.PairwiseAligner()
    reference = SeqIO.read(reference, "fasta")
    records_out = []
    metadata = []
    n = 0

    # Get metadata for pairwise alignments
    for record in SeqIO.parse(fasta_in, "fasta"):
        data = {
            "reference_id": reference.id,
            "reference_length": len(reference),
            "record_id": record.id,
            "record_length": len(record),
            "score_default": aligner.score(reference.seq, record.seq),                          # AS-IS score
            "score_reorient": aligner.score(reference.seq, record.seq.reverse_complement()),    # REV-COMP score
            "size_ratio": min(len(reference), len(record)) / max(len(reference), len(record)),
            "reoriented": "false"
        }

    # Filter, Reorient, and Write records to: fasta_out
        if data.get("size_ratio") >= 0.85:
            if data.get("score_reorient") > data.get("score_default"):
                record.seq = record.seq.reverse_complement()
                data['reoriented'] = "true"
                n += 1
            records_out.append(record)
            metadata.append(data)
    print(f"{n} records reoriented")
    print(f"{len(records_out)} records written to {fasta_out}")
    SeqIO.write(records_out, fasta_out, "fasta")

    # Write metadata to: csv
    print(f"{len(metadata)} records written to {csv} ")
    pd.DataFrame.from_records(metadata).to_csv(csv, index=False)


def kmer_reorient(K, reference, fasta_in, fasta_out, csv):
    print(f"Reorienting using {K}mers")
    reference = SeqIO.read(reference, "fasta")
    reference_kmers = get_kmers(str(reference.seq), K)
    records_out = []
    metadata = []
    n = 0

    # Get metadata for kmer matching
    for record in SeqIO.parse(fasta_in, "fasta"):
        fwd = get_kmers(str(record.seq), K)
        rev = get_kmers(str(record.seq.reverse_complement()), K)

        fwd_intersection = len(reference_kmers & fwd)
        rev_intersection = len(reference_kmers & rev)
        size_ratio = min(len(reference.seq), len(record.seq)) / max(len(reference.seq), len(record.seq))

        data = {
            "reference_id": reference.id,
            "record_id": record.id,
            "size_ratio": size_ratio,
            "K": K,
            "fwd_intersection": fwd_intersection,
            "rev_intersection": rev_intersection,
            "reoriented": "false"
        }

        # Filter, Reorient, and Write records to: fasta_out
        if data.get("size_ratio") >= 0.85:
            if data.get("rev_intersection") > data.get("fwd_intersection"):
                record.seq = record.seq.reverse_complement()
                data['reoriented'] = 'true'
                n += 1
            records_out.append(record)
            metadata.append(data)
    print(f"{n} records reoriented")
    print(f"{len(records_out)} records written to {fasta_out}")
    SeqIO.write(records_out, fasta_out, "fasta")

    print(f"{len(metadata)} records written to {csv} ")
    pd.DataFrame.from_records(metadata).to_csv(csv, index=False)


if __name__ == "__main__":
    REFERENCE = "/Users/beagles/dev/NJH/FlaskProjects/pipeline/data/ATCC19977.MAB_rpoB_partial.fasta"
    FASTA_IN = "/Users/beagles/dev/NJH/FlaskProjects/ntmhi/ntmhi/data/quintara_rpob.fna"
    FASTA_OUT = "out/quintara_filtered.fna"
    CSV = "out/rpob_pairwise.csv"

    pairwise_reorient(REFERENCE, FASTA_IN, FASTA_OUT, CSV)
    # kmer_reorient(5, REFERENCE, FASTA_IN, FASTA_OUT, CSV)

