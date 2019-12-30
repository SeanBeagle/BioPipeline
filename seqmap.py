from mmap import mmap
import json
import hashlib
from functools import wraps
from time import time

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import pandas as pd


FASTA_OLD = "/Users/beagles/dev/NJH/experimental/gmatrix/MAV_H87_matrix184.cfonly.fasta"
CSV = "data/matrixIndex.csv"
FASTA = "data/MAV_H87_matrix184.cfonly.CLEAN.fasta"
JSON = "data/matrixConfig.json"


def clean_fasta(fasta_in=FASTA_OLD, fasta_out=FASTA):
    """Remove newlines from sequences.
    This makes memory mapping easier since \n will no longer occupy positions within the sequence.
    TODO: Convert nucleotide pairs into single byte objects to reduce space by 0.5
    """
    with open(fasta_in, 'r') as file_in, open(fasta_out, 'a') as file_out:
        file_out.write(file_in.readline())
        for line in file_in:
            if line.startswith(">"):
                file_out.write("\n" + line)
            else:
                file_out.write(line.strip())


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        start = time()
        result = f(*args, **kw)
        end = time()
        elapsed = end-start
        if elapsed < 60:
            print(f'---> {f.__name__}() took {round(elapsed, 5)}s')
        elif 60 <= elapsed < 3600:
            print(f'---> {f.__name__}() took {round(elapsed/60, 5)}m')
        elif 3600 <= elapsed:
            print(f'---> {f.__name__}() took {round(elapsed/3600, 5)}h')
        return result
    return wrap


def checksum(file, algorithm, block=65536):
    try:
        if algorithm == 'md5':
            hasher = hashlib.md5()
        elif algorithm == 'sha1':
            hasher = hashlib.sha1()

        with open(file, 'rb') as fh:
            buffer = fh.read(block)
            while len(buffer) > 0:
                hasher.update(buffer)
                buffer = fh.read(block)
        return hasher.hexdigest()
    except Exception as e:
        print(e)


class FileMap:
    """Read mmap object using `with` statement."""
    def __init__(self, file, seek=0, access='ACCESS_DEFAULT'):
        file_access = {'ACCESS_READ': 'r+b', 'ACCESS_WRITE': 'w+b', 'ACCESS_COPY': 'r+b', 'ACCESS_DEFAULT': 'r+b'}
        mmap_access = ('ACCESS_READ', 'ACCESS_WRITE', 'ACCESS_COPY', 'ACCESS_DEFAULT')

        self.mmap_access = access if access in mmap_access else 'ACCESS_DEFAULT'
        self.file_access = file_access.get(access, 'r+b')
        self.filename = file
        self.seek = seek

    def __enter__(self):
        self.open_file = open(self.filename, self.file_access)
        self.mmap = mmap(self.open_file.fileno(), 0)
        self.mmap.seek(self.seek)
        return self.mmap

    def __exit__(self, exc_type, exc_value, tb):
        self.mmap.close()
        self.open_file.close()


class SeqMap:
    """Memory Mapped FastA file"""
    def __init__(self):
        self.fasta = ''
        self.md5 = ''
        self.sha1 = ''
        self._records = []
        self.header_map = {}

    def __repr__(self):
        return f"{__class__.__name__}('{self.fasta}', n={len(self._records)})"

    def __getitem__(self, item):
        if isinstance(item, int):
            return self._records[item]
        elif isinstance(item, str):
            if isinstance(self.header_map.get(item), int):
                return self._records[self.header_map.get(item)]

    @classmethod
    @timing
    def from_fasta(cls, fasta):
        print(f"...importing: {fasta}")
        obj = cls()
        obj.fasta = fasta
        obj.md5 = checksum(fasta, 'md5')
        obj.sha1 = checksum(fasta, 'sha1')

        with FileMap(fasta) as filemap:
            records = []
            record = None
            for i in range(filemap.size()):
                char = chr(filemap.read_byte())
                if char == '>':
                    if isinstance(record, MapRecord):
                        record.seq_size = i - record.seq_index - 1
                        records.append(record)
                        record = MapRecord(fasta, header_index=i)
                    else:
                        record = MapRecord(fasta, header_index=i)
                elif char == '\n' and not record.header_size:
                    record.header_size = i - record.header_index
                    record.seq_index = i + 1

            record.seq_size = filemap.size() - record.seq_index - 1
            records.append(record)
        obj._records = records
        obj._map_headers()
        obj._additional_tasks()
        return obj

    @classmethod
    def load(cls, file):
        data = json.load(open(file))
        if checksum(data['fasta'], 'md5') == data['md5'] and checksum(data['fasta'], 'sha1') == data['sha1']:
            obj = cls()
            obj.fasta = data.get('fasta')
            obj.sha1 = data.get('sha1')
            obj.md5 = data.get('md5')
            obj.header_map = data.get('header_map')
            for record in data['records']:
                obj.records.append(MapRecord(obj.fasta, **record))
            return obj
        else:
            raise ImportError(f"Checksum does not match!\nTry: {cls.__name__}.read_fasta() instead")

    def save(self, file):
        data = {
            'fasta': self.fasta,
            'md5': self.md5,
            'sha1': self.sha1,
            'header_map': self.header_map,
            'records': [record.to_dict() for record in self._records]
        }
        json.dump(data, open(file, 'w'))

    def _map_headers(self):
        self.header_map = {record.header: i for (i, record) in enumerate(self._records)}

    def _additional_tasks(self):
        return None

    def subset(self, records=None, indexes=None, headers=None, startswith='', endswith='', contains=''):
        obj = self.__class__()
        obj.fasta = self.fasta
        obj.md5 = self.md5
        obj.sha1 = self.sha1

        if records:
            obj._records = records
        elif indexes:
            obj._records = [self[i] for i in indexes if self[i]]
        elif headers:
            obj._records = [self[header] for header in headers if self[header]]
        else:
            obj._records = [
                record for record in self._records
                if record.header.startswith(startswith)
                and record.header.endswith(endswith)
                and contains in record.header
            ]

        obj._map_headers()
        obj._additional_tasks()
        return obj

    @property
    def records(self):
        return (record for record in self._records)

    def record(self, record):
        return self[record]

    def filter_records(self, startswith='', endswith='', contains=''):
        records = [
            record for record in self._records
            if record.header.startswith(startswith)
            and record.header.endswith(endswith)
            and contains in record.header
        ]
        return records


class SeqMatrix(SeqMap):
    """SeqMap of aligned FastA"""
    def __init__(self):
        super().__init__()

    def __repr__(self):
        return f"{__class__.__name__}('{self.fasta}', {self.shape})"

    @classmethod
    def load(cls, file=JSON):
        data = json.load(open(file))
        if checksum(data['fasta'], 'md5') == data['md5'] and checksum(data['fasta'], 'sha1') == data['sha1']:
            obj = cls()
            obj.fasta = data.get('fasta')
            obj.sha1 = data.get('sha1')
            obj.md5 = data.get('md5')
            obj.header_map = data.get('header_map')
            for record in data['records']:
                obj._records.append(MapRecord(matrix.fasta, **record))
            return obj
        else:
            raise ImportError("Checksum does not match!\nTry: SeqMatrix.read_fasta() instead")

    def save(self, file=JSON):
        data = {
            'fasta': self.fasta,
            'md5': self.md5,
            'sha1': self.sha1,
            'header_map': self.header_map,
            'records': [record.to_dict() for record in self.records]
        }
        json.dump(data, open(file, 'w'))

    @property
    def shape(self):
        return (len(self._records), self._records[0].seq_size) if len(self._records) > 0 else (0, 0)

    @property
    def positions(self):
        return (self.position(i) for i in range(self.shape[1]))

    def position(self, position):
        return Position(self, position)

    def is_matrix(self):
        return set(len(record) for record in self._records) == 1


class MapRecord:
    def __init__(self, fasta, **kwargs):
        self.fasta = fasta
        self.header_index = kwargs.get('header_index')
        self.header_size = kwargs.get('header_size')
        self.seq_index = kwargs.get('seq_index')
        self.seq_size = kwargs.get('seq_size')
        self._counts = kwargs.get('counts', {})
        self.unique_gaps = kwargs.get('unique_gaps')
        self.contribution_score = kwargs.get('contribution_score')

    def __getitem__(self, item):
        with FileMap(self.fasta) as filemap:
            obj = None
            if isinstance(item, int) and item < self.seq_size:
                filemap.seek(self.seq_index + item)
                obj = chr(filemap.read_byte())
            elif isinstance(item, slice):
                filemap.seek(self.seq_index + item.start)
                sequence = filemap.read(min(self.seq_size, item.stop) - item.start).decode('utf-8')
                seq = Seq(sequence, generic_dna)
                obj = SeqRecord(seq, id=self.header, name=self.header, description=f"{item.start + 1}:{item.stop + 1}")
            elif isinstance(item, str) and len(item) == 1:
                obj = self.counts.get(item, 0)
            return obj

    def __repr__(self):
        return self.header

    def __len__(self):
        return self.seq_size or 0

    @property
    def header(self):
        if isinstance(self.header_index, int) and isinstance(self.header_size, int):
            with FileMap(self.fasta, self.header_index) as filemap:
                return filemap.read(self.header_size).decode('utf-8')[1:]
        else:
            return ""

    @property
    def seq(self):
        with FileMap(self.fasta, self.seq_index) as filemap:
            sequence = filemap.read(self.seq_size).decode('utf-8')
            seq = Seq(sequence, generic_dna)
            return SeqRecord(seq, id=self.header, name=self.header)

    @property
    def counts(self):
        if not self._counts:
            with FileMap(self.fasta, self.seq_index) as fmap:
                for char in fmap.read(self.seq_size).decode('utf-8'):
                    self._counts[char] = self._counts.get(char, 0) + 1
        return self._counts

    def position(self, item):
        return self.__getitem__(item)

    @property
    def positions(self):
        return (self.position(i) for i in range(len(self)))

    def count(self, char):
        return self.counts.get(char) or 0

    def to_dict(self):
        return {'header_index': self.header_index, 'header_size': self.header_size,
                'seq_index': self.seq_index, 'seq_size': self.seq_size, 'counts': self._counts}


class Position:
    def __init__(self, matrix, index):
        self.matrix = matrix
        self.index = index
        self.counts = {}
        self.gap_records = []
        self.c_value = 0
        self.is_snp = None
        self.is_core_snp = None
        self._position_array = (record[index] for record in matrix.records)
        self._build_position()

    def __repr__(self):
        return f"{self.__class__.__name__}(index={self.index}, counts={self.counts})"

    def __str__(self):
        return f"{self.__class__.__name__}(index={self.index}, counts={self.counts})"

    def __contains__(self, item):
        return item in self.counts

    def _build_position(self):
        for char in self._position_array:
            self.counts[char] = self.counts.get(char, 0) + 1
            self.gap_records.append(char == '-')

        self.is_core_snp = (self.counts.get('-', 0) == 0 and len(self.counts) > 1)
        self.is_snp = self.is_core_snp or len(self.counts) > 2
        self.c_value = self.counts.get('-', 0) / len(self.gap_records)


class GFFMap:
    headers = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def __init__(self, file):
        self.gff = file
        # self.header = self._import_header()
        self.features = self._import_features()

    def _import_header(self):
        header = {}
        with open(self.gff) as gff:
            for line in gff:
                if line.startswith("#"):
                    if "gff-spec-version" in line:
                        print(line.split('version '))
                        header['version'] = line.strip().split('version ')[1]
                    elif "!processor" in line:
                        header['processor'] = line.strip().split('!processor ')[1]
                    elif "!genome-build " in line:
                        header['genome_build'] = line.split("!genome_build ")[1]
                    elif "!genome-build-accession " in line:
                        header['genome_build_accession'] = line.split('accession ')[1]
                    elif "!annotation-date" in line:
                        header['annotation_date'] = line.split('annotation-date ')[1]
                    elif "!annotation-source" in line:
                        header['annotation_source'] = line.split('annotation_source ')[1]
                    elif "sequence-region" in line:
                        header['sequence_region'] = line.split('sequence-region ')[1]
                    elif "species" in line:
                        header['species'] = line.split('species ')[1]
                else:
                    break
        return header

    def _import_features(self):
        records = []
        with open(self.gff) as gff:
            for line in gff:
                if not line.startswith('#') and line.count('\t') == 8:
                    seqname, source, feature, start, end, score, strand, frame, attribute = line.strip().split('\t')
                    if feature != 'region':
                        record = {
                            'seqname': seqname,
                            'source': source,
                            'feature': feature,
                            'start': int(start),
                            'end': int(end),
                            'score': score,
                            'strand': strand,
                            'frame': frame,
                        }
                        record.update({k: v for (k, v) in [pair.split('=') for pair in attribute.split(';')]})

                        records.append(record)
        return pd.DataFrame.from_records(records)

    def locate_snp(self, snp):
        snp += 1
        return self.features[(self.features['start'] < snp) & (snp < self.features['end'])]


if __name__ == '__main__':
    matrix = SeqMatrix().from_fasta(FASTA)




