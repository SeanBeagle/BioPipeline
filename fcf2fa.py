from filehandler import Dir
import pandas as pd


"""
beagles[compbio] alma: grep , E-ALMA1.MAV_vs_NJH87.cf | more

CP018363.1	83972	C	T,G	0	34	34	SNP
CP018363.1	148672	A	G,C	0	14	14	SNP
CP018363.1	1381815	G	A,C	0	8	8	SNP
CP018363.1	2334023	G	A,T	0	30	30	SNP
CP018363.1	2948951	T	C,G	0	6	6	SNP
CP018363.1	3621547	T	C,G	0	35	35	SNP
CP018363.1	4243092	A	G,C	0	33	33	SNP
CP018363.1	4562156	T	C,G	0	15	15	SNP
"""


VCF_DIR = Dir("/Strong/proj/.data/alma")
FASTA_DIR = VCF_DIR.make_subdir("fasta")
MATRIX = FASTA_DIR.join("ALMA_matrix_N.fna")
CSV_OUT = FASTA_DIR.join("ALMA_stats.csv")
CSV_MUTATIONS = FASTA_DIR.join("ALMA_mutations.csv")
POSITIONS = 5626623

files = VCF_DIR.files(endswith='.cf')

records = []
multi_alleles = {}

print(f"Building Matrix: {MATRIX} from {len(files)} files.")
for i, file in enumerate(files):
    print(f"\t{i:02d} | {file.filename}")
    isolate = file.filename.split("_vs_")[0]
    seq = ['-'] * POSITIONS
    record = {'isolate': isolate, 'totalMulti': 0, 'N': 0, 'alt': 0, 'ref': 0}

    options=set()
    with open(file.path, 'r') as file_in, open(MATRIX, 'a+') as file_out:
        for line in file_in:
            # EXTRACT METADATA FROM LINE
            chrom, pos, ref, alt, ref_depth, alt_depth, best_depth, ref_or_snp = line.strip().split('\t')

            # TRANSFORM VALUES
            is_snp = (ref_or_snp == 'SNP')
            pos = int(pos) - 1                  # VCF index starts at 1
            ref_depth = int(ref_depth)
            alt_depth = int(alt_depth)

            # COUNT ALL MULTI-ALLELES
            record['totalMulti'] += (len(alt) > 1)

            # ASSIGN BASE TO POSITION
            if is_snp and len(alt) == 1:
                seq[pos] = alt
                record['alt'] += 1
            elif is_snp and len(alt) > 1:
                seq[pos] = 'N'
                record['N'] += 1
                multi_alleles[pos] = multi_alleles.get(pos, 0) + 1
            else:
                seq[pos] = ref
                record['ref'] += 1

        print(f"\t\t{record}")
        # WRITE SEQUENCE TO MATRIX FILE (2-line FastA format)
        print(f">{isolate}\n{''.join(seq)}", file=file_out)
    records.append(record)


pd.DataFrame.from_records(records).to_csv(CSV_OUT, index=False)
pd.DataFrame([multi_alleles.keys(), multi_alleles.values()], columns=['pos', 'freq']).to_csv(CSV_MUTATIONS, index=False)

with open(CSV_MUTATIONS, "a") as file:
    print("position,isolates", file=file)
    for (k, v) in multi_alleles.items():
        print(f"{k},{v}", file=file)
