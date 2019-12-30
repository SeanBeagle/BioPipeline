import pandas as pd
from Bio import SeqIO


class VCF:
    def __init__(self, filename):
        # ATTRIBUTES
        self.filename = filename
        self.version = ''
        self.contigs = {}
        self.columns = []
        self.records = []
        self.multi = set()

        # CONSTRUCT
        self.read_header()
        self.parse_positions()
        self.dataframe = pd.DataFrame.from_records(self.records)
        self.sequence = "-" * sum(self.contigs.values())

    def read_header(self):
        contigs = {}
        with open(self.filename) as fh:
            for line in fh:
                if line.startswith("##fileformat"):
                    self.version = line.strip().split("=")[1]
                elif line.startswith("##contig"):
                    id = line.split("ID=")[1].split(",")[0]
                    length = int(line.split("length=")[1].split(">")[0])
                    contigs[id] = length
                elif line.startswith("#CHROM"):
                    self.columns = line[1:].strip().split('\t')
                elif not line.startswith("#"):
                    self.contigs = contigs
                    return

    def filter_position(self, position, records):
        if len(records) > 1:
            self.multi.add(position)


    def parse_positions(self):
        position = -1
        records = []
        with open(self.filename) as fh:
            for line in fh:
                if not line.startswith("#"):
                    record = (self.clean(line))
                    if record['POS'] == position or position == -1:
                        records.append(record)
                        position = record['POS']
                    else:
                        self.filter_position(position, records)
                        position = record['POS']
                        records = [record]

    def clean(self, line):
        line = line.strip().split('\t')
        data = {
            self.columns[i]: line[i]
            for i in range(len(line))
            if len(self.columns) == len(line)
        }
        # data['INFO'] = data.get('INFO').split(";")
        # data['INFO']['DP'] = int(data.get('INFO').get('DP'))
        # data['INFO']['MMQ0F'] = float(data.get('INFO').get('MQ0F'))
        # data['INFO']['FQ'] = float(data.get('INFO').get('FQ'))
        # data['INFO']['AC1'] = int(data.get('INFO').get('AC1'))
        # data['INFO']['MQ0F'] = float(data.get('INFO').get('MQ0F'))
        # data['INFO']['MQ'] = int(data.get('INFO').get('MQ'))
        # data['INFO']['DP4'] = [int(i) for i in data.get('INFO').get('DP4').split(",")]
        data['DP4'] = data.get('INFO').split('DP4=')[1].split(';')[0].split(',')
        data['QUAL'] = float(data.get('QUAL'))
        data['POS'] = int(data.get('POS'))
        return data


if __name__ == '__main__':
    vcf = VCF('/Users/beagles/dev/NJH/FlaskProjects/pipeline/data/NJH87-ALMA-1-10_vs_NJH87.cf')
    print(vcf.sequence[:100])
    print(len(vcf.sequence))
    print(len(vcf.multi))




