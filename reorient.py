"""
AUTHOR: Sean Beagle
URL: https://SeanBeagle.github.io
"""
import argparse
from Bio import SeqIO, Align


_aligner = Align.PairwiseAligner()


def reorient_file(fasta_in, fasta_out):
    """Write records from `fasta_in` to `fasta_out`  """
    reference = next(SeqIO.parse(fasta_in, 'fasta'))
    records = SeqIO.parse(fasta_in, 'fasta')
    oriented_records = (reorient_record(record, reference) for record in records)
    SeqIO.write(oriented_records, fasta_out, 'fasta')


def reorient_record(record, reference):
    """Return the orientation of `record` with best global pairwise score to `reference`."""
    fwd_score = _aligner.score(reference.seq, record.seq)
    rev_score = _aligner.score(reference.seq, record.seq.reverse_complement())
    if rev_score > fwd_score:
        print(f'Reorienting {record.id}')
        record.seq = record.seq.reverse_complement()
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('-i', '--input', required=True, metavar='FASTA',
                        help='multi-fasta with reference as first record')
    parser.add_argument('-o', '--output', required=True, metavar='FASTA',
                        help='filename to write reoriented records')
    args = parser.parse_args()
    reorient_file(fasta_in=args.input, fasta_out=args.output)

