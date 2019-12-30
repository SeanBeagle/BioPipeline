import os
import shutil
import sys
import re
import subprocess
import gzip
from mmap import mmap
import csv
import hashlib
from statistics import mean, stdev
from types import SimpleNamespace as namespace
import datetime
from statistics import StatisticsError, median
import multiprocessing as mp

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from seqmap import timing


class Dir:
    """ Base class for system directories """

    def __init__(self, path=os.getcwd()):
        self._path = None

        if isinstance(path, str):
            pass
        elif isinstance(path, self.__class__):
            path = path.path
        else:
            raise TypeError(f"path must be string, NOT {type(path)}")

        self.path = path.rstrip('/')

    def __repr__(self):
        return self.path

    def __str__(self):
        return self.path

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not os.path.isabs(value):
            value = os.path.join(os.getcwd(), value)
        if os.path.isdir(value):
            self._path = value
        else:
            raise NotADirectoryError(value)

    @property
    def dirname(self):
        """Return string of relative directory name"""
        return self.path.split('/')[-1]

    @property
    def parents(self):
        """Return list of parent directories to root"""
        components = list(os.path.dirname(self.path).strip('/').split('/'))
        components[0] = "/" + components[0]
        return [Dir("/".join(components[0:len(components) - i])) for i, parent in (enumerate(components))]

    @property
    def parent(self):
        return self.parents[0]

    @property
    def siblings(self):
        return self.parent.subdirs

    @property
    def subdirs(self):
        return [Dir(self.join(subdir)) for subdir in os.listdir(self.path) if os.path.isdir(self.join(subdir))]

    @property
    def recursive(self):
        """"Return list of all recursive subdirectories"""
        def get_recursive(directory):
            subdirs.extend(directory.subdirs)
            try:
                for subdir in directory.subdirs:
                    get_recursive(subdir)
            except PermissionError:
                return

        subdirs = [self]
        get_recursive(self)
        return subdirs

    @property
    def all_files(self):
        """Return all files in directory as File() objects"""
        try:
            return [File(self.join(file)) for file in os.listdir(self.path) if os.path.isfile(self.join(file))]
        except PermissionError:
            return

    def files(self, startswith='', endswith='', contains='', recursive=False, dataframe=False):
        """Return filtered list of File objects

        parameters
        ----------
        startswith : str

        :param startswith: str
        :param endswith:
        :param contains:
        :param recursive:
        :param dataframe:
        :return:
        """

        def meets_conditions(f):
            return f.filename.startswith(startswith) and f.filename.endswith(endswith) and contains in f.filename

        if recursive:
            files = [file for subdir in self.recursive for file in subdir.all_files if meets_conditions(file)]
        else:
            files = [file for file in self.all_files if meets_conditions(file)]

        if dataframe:
            return objects_as_dataframe(files, columns=['filename', 'path', 'filesize', 'inode'])
        else:
            return files

    def join(self, *args):
        return os.path.join(self.path, *args)

    def move(self, path):
        if isinstance(path, Dir):
            path = path.path
        elif isinstance(path, str):
            path = path
        else:
            raise TypeError(f"path must be str or Dir. NOT {type(path)}")
        try:
            _newpath = os.path.join(path, self.dirname)
            shutil.move(self.path, _newpath)
            self.path = _newpath
        except Exception as e:
            print(e)

    def rmtree(self):
        shutil.rmtree(self.path)

    def make_subdir(self, *args):
        """ Makes recursive subdirectories from 'os.path.join' like arguments """
        subdir = self.join(*args)
        return self.make(subdir)

    @classmethod
    def make(cls, path):
        try:
            os.makedirs(path)
            return cls(path)
        except FileExistsError:
            if os.path.isfile(path):
                raise FileExistsError(path)
            else:
                return cls(path)

    @classmethod
    def verify_or_make(cls, path):
        """ Verifies that ``path`` is an existing directory or tries to make it """
        if not os.path.isabs(path):
            path = os.path.join(os.getcwd(), path)

        if os.path.isdir(path):
            try:
                test = os.path.join(path, 'testing_permissions')
                subprocess.call(['touch', f'{test}'])
                subprocess.call(['rm', f'{test}'])
                return cls(path)
            except Exception as e:
                print(e)
                sys.exit(f"Unable to use directory: {path}, exception: {e}")
        else:
            try:
                subprocess.call(['mkdir', '-p', f"{path}"])
                return cls(path)
            except Exception as e:
                print(e)
                sys.exit(f"Unable to make directory: {path}")


def objects_as_dataframe(files, columns):
    return pd.DataFrame.from_records([{attribute: getattr(f, attribute) for attribute in columns} for f in files])


class File:
    """ Base class for all file-types"""
    extensions = ''

    def __init__(self, path, file_type=None):
        self._path = None
        self.path = path
        self.file_type = file_type

    def __repr__(self):
        return self.path

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.md5 == other.md5) & (self.sha1 == other.sha1)
        else:
            return False

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not os.path.isabs(value):
            value = os.path.join(os.getcwd(), value)
        if os.path.isfile(value):
            self._path = value
        else:
            raise FileNotFoundError(value)

    @property
    def dir(self):
        return Dir(os.path.dirname(self.path))

    @property
    def dirname(self):
        return self.dir.dirname

    @property
    def md5(self):
        hasher = hashlib.md5()
        blocksize = 65536
        with open(self.path, 'rb') as fh:
            buf = fh.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = fh.read(blocksize)
        return hasher.hexdigest()

    @property
    def sha1(self):
        hasher = hashlib.sha1()
        blocksize = 65536
        with open(self.path, 'rb') as fh:
            buf = fh.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = fh.read(blocksize)
        return hasher.hexdigest()

    @property
    def relpath(self):
        return self.dir.path

    @property
    def filename(self):
        return os.path.basename(self.path)

    @property
    def file_prefix(self):
        return self.filename.split(".")[0]

    @property
    def extension(self):
        return self.filename.split(".")[-1]

    @property
    def device(self):
        return os.stat(self.path).st_dev

    @property
    def inode(self):
        return os.stat(self.path).st_ino

    @property
    def hardlinks(self):
        return os.stat(self.path).st_nlink

    @property
    def is_symlink(self):
        return os.path.islink(self.path)

    @property
    def filesize(self):
        return os.stat(self.path).st_size

    @property
    def filesize_hr(self, suffix='B'):
        num = self.filesize
        for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
            if abs(num) < 1024.0:
                return "%3.1f%s%s" % (num, unit, suffix)
            num /= 1024.0
        return "%.1f%s%s" % (num, 'Yi', suffix)

    def rename(self, filename):
        """ Rename file as``filename`` and return File(object) """
        new_path = os.path.join(self.relpath, filename)
        if os.path.isfile(new_path):
            raise FileExistsError(new_path)
        else:
            os.rename(self.path, new_path)
            self.__init__(new_path)
        return self

    def copy(self, filename, directory=None):
        """ Copy file as ``filename`` and return File(object) """
        if directory:
            directory = Dir(directory)

        new_file = os.path.join(directory or self.relpath, filename)
        if not os.path.isfile(new_file):
            shutil.copyfile(self.path, new_file)
        else:
            raise FileExistsError(new_file)
        return self.__class__(new_file)

    def move(self, path):
        """ Move file to directory: ``path`` """
        if isinstance(path, Dir):
            path = path.path
        elif not os.path.isabs(path):
            path = os.path.join(self.relpath, path)
        new_path = os.path.join(path, self.filename)
        Dir.verify_or_make(path)
        subprocess.call(["mv", self.path, new_path])
        self.path = new_path
        return self

    def hash(self, hasher, blocksize=65536):
        with open(self.path, 'rb') as fh:
            buf = fh.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = fh.read(blocksize)
        return hasher.hexdigest()

    def gunzip(self):
        subprocess.run(['gunzip', self.path])
        self.path = ".".join(self.path.split(".")[:-1])

    def gzip(self):
        subprocess.run(['gzip', self.path])
        self.path = self.path + ".gz"

    def compare_filesize(self, other):
        """ Compare file sizes to another file or File(object): ``other`` """
        try:
            if isinstance(other, self.__class__):
                return round(self.filesize / other.filesize, 2)
            else:
                return round(self.filesize / File(other).filesize, 2)
        except Exception as e:
            sys.exit(e)

    @classmethod
    def get_all(cls, path=os.getcwd(), recursive=False):
        """Return a list of files in ``path`` as ``File`` objects."""
        if isinstance(path, Dir):
            directory = path
        elif isinstance(path, str):
            directory = Dir(path)
        else:
            raise TypeError(f"path must be of type Dir or str. NOT {type(path)}")
        return [cls(file.path) for file in directory.files(endswith=cls.extensions, recursive=recursive)]


class Fastq(File):
    """Create FastQ file objects."""
    extensions = ('fastq', 'fastq.gz', 'fq', 'fq.gz')

    def __init__(self, path):
        super().__init__(path=path, file_type="fastq")
        self._stats = {'set': False, 'read_count': None, 'read_lengths': None, 'avg_read_len': None, 'std_read_len': None}
        self._reads = None

        # SEE IF THE FILE'S NAMING CONVENTION MATCHES THAT OF RAW ILLUMINA READS
        try:
            r_illumina = r'\b.*_S\d+_L.{3}_R[12]_\d{3}\.fastq\.gz'
            match = re.search(r_illumina, self.filename)
            if match[0] == self.filename:
                self.sample_name = self.file_prefix.split("_")[0]
                self.sample_number = self.file_prefix.split("_")[1]
                self.lane = self.file_prefix.split("_")[2]
                self.read = self.file_prefix.split("_")[3]
                self.last_segment = self.file_prefix.split("_")[4]
        except TypeError:
            pass

    @property
    def reads(self):
        if not self._reads:
            n = 0
            for _ in gzip.open(self.path, "rt"):
                n += 1
            self._reads = n / 4
        return self._reads

    @property
    def read_lengths(self, stat='read_lengths'):
        self.get_stats()
        return self._stats.get(stat)

    @property
    def avg_read_len(self, stat='avg_read_len'):
        self.get_stats()
        return self._stats.get(stat)

    @property
    def std_read_len(self, stat='std_read_len'):
        self.get_stats()
        return self._stats.get(stat)

    def get_stats(self):
        if not self._stats['set']:
            with gzip.open(self.path, "rt") as handle:
                n = [len(record) for record in SeqIO.parse(handle, self.file_type)]
                self._stats['reads'] = len(n)
                self._stats['avg_read_len'] = mean(n)
                self._stats['std_read_len'] = stdev(n, self._stats['avg_read_len'])
                self._stats['read_lengths'] = n
                self._stats['set'] = True

    @classmethod
    def get_pairs(cls, directory=os.getcwd(), r1="_R1_", r2="_R2_",):
        """Return list of Fastq objects of paired-end reads in directory: `path`"""
        if isinstance(directory, Dir):
            pass
        elif isinstance(directory, str):
            directory = Dir(directory)
        else:
            raise TypeError(f"path must be of type str or Dir. NOT: {type(path)}")

        read1 = [Fastq(file.path) for file in directory.files(endswith=cls.extensions, contains=r1)]
        return [
            namespace(pair1=pair1, pair2=Fastq(pair1.path.replace(r1, r2)))
            for pair1 in read1 if os.path.exists(pair1.path.replace(r1, r2))]

    @classmethod
    def get_trimmed_pairs(cls, path=os.getcwd()):
        return cls.get_pairs(path=path, r1="pair1", r2="pair2")

    def compare_read_count(self, other):
        """ Compare read_counts to another file or Fastq(object)"""
        try:
            if isinstance(other, self.__class__):
                return round(self.reads / other.reads, 2)
            else:
                return round(self.reads / Fastq(other).reads, 2)
        except ZeroDivisionError:
            print(f"ZeroDivisionError: {self.file_prefix}:{self.reads} {self.file_prefix}: {other.reads}")
        except Exception as e:
            sys.exit(e)

    def compare_avg_len(self, other):
        """ Compare average read length to another fastq file or existing Fastq(object) """
        try:
            if isinstance(other, self.__class__):
                return round(self.avg_read_len / other.avg_read_len, 2)
            else:
                return round(self.avg_read_len / Fastq(other).avg_read_len, 2)
        except ZeroDivisionError:
            print(f"ZeroDivisionError: {self.file_prefix}:{self.avg_read_len} {self.file_prefix}: {other.avg_read_len}")

    def compare_std_len(self, other):
        """ Compare standard deviation of read lengths to another fastq file or existing Fastq(object) """
        try:
            if isinstance(other, self.__class__):
                return round(self.std_read_len / other.std_read_len, 2)
            else:
                return round(self.std_read_len / Fastq(other).std_read_len, 2)
        except ZeroDivisionError:
            print(f"ZeroDivisionError: {self.file_prefix}:{self.std_read_len} {self.file_prefix}: {other.std_read_len}")

    @classmethod
    def combine_lanes(cls, path=os.getcwd(), minreads=1500000):
        for pair in cls.get_pairs(path):
            if pair.pair1.lane == "L001" and os.path.isfile(pair.pair1.path.replace("L001", "L002")):

                # FILES AND DIRECTORIES
                L001_R1, L001_R2 = pair.pair1, pair.pair2
                lane1 = L001_R1.reads + L001_R2.reads

                L002_R1, L002_R2 = Fastq(L001_R1.path.replace("L001", "L002")), Fastq(L001_R2.path.replace("L001", "L002"))
                lane2 = L002_R1.reads + L002_R2.reads
                archive = Dir.verify_or_make(os.path.join(path, "000_Archive"))
                quarantine = Dir.verify_or_make(os.path.join(path, "000_Quarantine"))

                if lane1 >= minreads:
                    # If L001 meets minimum read count
                    with open('README', 'a+') as readme:
                        print("L001 ONLY:", L001_R1.sample_name, file=readme)
                    L001_R1.copy(L001_R1.filename.replace("L001", "LCP1"))
                    L001_R2.copy(L001_R2.filename.replace("L001", "LCP1"))
                    L001_R1.move(archive.path)
                    L001_R2.move(archive.path)
                    L002_R1.move(archive.path)
                    L002_R2.move(archive.path)

                elif lane1 < minreads and lane1 + lane2 >= minreads:
                    # L001 doesn't meet minimum read count, but L001 + L002 do... merge L001+L002
                    print(f"combining: {L001_R1.sample_name}")
                    LCAT_R1 = L001_R1.copy(L001_R1.filename.replace("L001", "L1X2"))
                    LCAT_R2 = L001_R2.copy(L001_R2.filename.replace("L001", "L1X2"))
                    with open('README', 'a+') as readme:
                        print("L001+L002:", L001_R1.sample_name, file=readme)
                    subprocess.Popen(
                        f"bsub 'cat {L002_R1.path} >> {LCAT_R1.path}; mv {L001_R1.path} {L002_R1.path} {archive.path}'",
                        shell=True)
                    subprocess.Popen(
                        f"bsub 'cat {L002_R2.path} >> {LCAT_R2.path}; mv {L001_R2.path} {L002_R2.path} {archive.path}'",
                        shell=True)



                elif lane1 < minreads and lane1 + lane2 < minreads:
                    # If lane1 doesn't meet minimum and lane1 + lane2 don't either, quarantine all reads.
                    L001_R1.move(quarantine.path)
                    L001_R2.move(quarantine.path)
                    L002_R1.move(quarantine.path)
                    L002_R2.move(quarantine.path)

    @classmethod
    def did_not_assemble(cls, path=os.getcwd()):
        for pair in Fastq.get_pairs():
            if not os.path.isfile("/Strong/proj/.data/Morty/data/02_assembled/2019-03-14_hiseq" + pair.pair1.filename.split('_')[1] + ".000.000.fasta"):
                print(pair.pair1.filename.split('_')[1])


class Fasta(File):
    """Fasta File Type"""
    extensions = ('fasta', 'fa', 'fna')

    def __init__(self, path):
        super().__init__(path=path, file_type="fasta")

        self.stats = {
            'contigs': None,
            'bp': None,
            'avg_contig_len': None,
            'std_contig_len': None,
            'nt_freq': None
        }

    def _get_stat(self, stat):
        if not self.stats.get(stat):
            self.get_stats()
        return self.stats.get(stat)

    @property
    def contigs(self, stat='contigs'):
        if self.stats[stat] is None:
            self.get_stats()
        return self.stats[stat]

    @property
    def bp(self, stat='bp'):
        if self.stats[stat] is None:
            self.get_stats()
        return self.stats[stat]

    @property
    def avg_contig_len(self, stat='avg_contig_len'):
        if self.stats[stat] is None:
            self.get_stats()
        return self.stats[stat]

    @property
    def std_contig_len(self, stat='std_contig_len'):
        if not self.stats.get(stat):
            self.get_stats()
        return self.stats[stat]

    @property
    def nt_freq(self, stat='nt_freq'):
        if not self.stats.get(stat):
            self.get_stats()
        return self.stats[stat]

    def get_stats(self):
            n = []
            nt_freq = {'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N': 0, '-': 0}
            for record in SeqIO.parse(self.path, self.file_type):
                n.append(len(record))
                for nt in record.seq:
                    nt_freq[nt.upper()] += 1

            self.stats['bp'] = sum(n)
            self.stats['contigs'] = len(n)
            self.stats['avg_contig_len'] = int(mean(n))
            self.stats['std_contig_len'] = int(stdev(n, self.stats['avg_contig_len']))
            self.stats['nt_freq'] = nt_freq

    def chunk(self, out, chunksize=1024):
        """ Chunk Fasta file into records of size ``base_pairs`` and save as file ``out`` """
        with open(out, 'a+') as fh:
            for r in SeqIO.parse(self.path, self.file_type):
                start = 0
                chunk = 1
                while start < len(r):
                    if start + chunksize < len(r):
                        end = start + chunksize
                    elif start + chunksize > len(r):
                        end = len(r)

                    record = SeqRecord(
                        r.seq[start:end],
                        id=f"{r.id}.{chunk}",
                        description=f"({start}-{end}) of {len(r)}bp")
                    SeqIO.write(record, fh, self.file_type)
                    start = end + 1
                    chunk += 1

    @staticmethod
    def get_reports():
        for fa in Fasta.get_all():
            with open('fasta_report.log', 'a+') as fh:
                print(f"{fa.filename.split('.')[0]},{fa.contigs},{fa.bp}", file=fh)

    @staticmethod
    def rename_reads():
        for fasta in Fasta.get_all():
            project = fasta.dir.dirname
            trimmed_root = Dir("/Strong/proj/.data/Morty/data/01_trimmed/")
            pairs = Fastq.get_pairs(os.path.join(trimmed_root.path, project))
            sample, taxon, level, ext = fasta.filename.split(".")
            if taxon != '000':
                for pair in pairs:
                    if pair.pair1.sample_name == sample:
                        trim1 = pair.pair1
                        trim2 = pair.pair2
                        trim1.rename(f"{trim1.sample_name}.{taxon}.trimmed-pair1.fastq.gz")
                        trim2.rename(f"{trim2.sample_name}.{taxon}.trimmed-pair2.fastq.gz")


class ABI(File):
    extensions = ("ab1", )
    DIR = "/Users/beagles/dev/NJH/FlaskProjects/ntmhi/ntmhi/data/2019-11-12_577059_99831_dataab1/"
    FILE = "/Users/beagles/dev/NJH/FlaskProjects/ntmhi/ntmhi/data/2019-11-12_577059_99831_dataab1/Mycobacterial_rpoB_A01_Myco_rpoB_R_2019-11-15_A01.ab1"

    def __init__(self, path=None):
        super().__init__(path=path or ABI.FILE, file_type="abi")
        self.record = SeqIO.read(self.path, "abi")
        self.seq = self.record.seq
        self.well = self.record.annotations.get('sample_well')
        self.q_scores = self.record.letter_annotations.get('phred_quality')
        self.trim_range = [0, len(self.record.seq)]
        hour1, min1, sec1 = self.record.annotations.get('abif_raw').get('RUNT1').split(":")
        hour2, min2, sec2 = self.record.annotations.get('abif_raw').get('RUNT2').split(":")
        t1 = datetime.timedelta(hours=int(hour1), minutes=int(min1), seconds=int(sec1))
        t2 = datetime.timedelta(hours=int(hour2), minutes=int(min2), seconds=int(sec2))
        self.duration = round((t2-t1).total_seconds() / 60**2, 2)

    def __len__(self):
        return len(self.seq)

    def quality_plot(self):
        q_scores = pd.Series(self.record.letter_annotations.get('phred_quality'))
        q_mean = q_scores.median()
        l_trim = self.quality_trim('left')
        r_trim = self.quality_trim('right')
        q_trim_mean = pd.Series([q for i, q in enumerate(q_scores) if l_trim < i < r_trim]).median()

        fig, ax = plt.subplots()
        ax.axhspan(0, 10, facecolor='r', alpha=0.2)
        ax.axhspan(10, 20, facecolor='y', alpha=0.2)
        ax.axhspan(20, 60, facecolor='g', alpha=0.2)
        ax.set_ylabel("Q-score")
        ax.set_xlabel("position")
        ax.set_title(self.well)
        q_scores.plot(ax=ax, color='black', linewidth=0.5)
        ax.axhline(y=q_mean, color='red', linestyle='-', linewidth=1.5)             # AVERAGE Q
        ax.axhline(y=q_trim_mean, color='black', linestyle=':', linewidth=1.5)     # TRIMMED AVERAGE Q
        # ax.axvline(x=self.quality_trim('left'), color='red')
        # ax.axvline(x=self.quality_trim('right'), color='red')
        ax.axvspan(0, self.quality_trim('left'), facecolor='red', alpha=0.5)                   # TRIM-ZONE: LEFT
        ax.axvspan(self.quality_trim('right'), len(self.q_scores), facecolor='r', alpha=0.5)   # TRIM-ZONE: RIGHT
        return fig

    @classmethod
    def parse_directory(cls, dir=os.getcwd(), out="quality_plots.pdf"):
        dir = dir or ABI.DIR
        records = ABI.get_all(dir)
        pdf = matplotlib.backends.backend_pdf.PdfPages(out)
        for record in records:
            if sum(record.q_scores) > 0:
                fig = record.quality_plot()
                pdf.savefig(fig)
                plt.close(fig)
        pdf.close()

    def quality_trim(self, side, threshold=10):
        window = 10
        pos = 0 if side == 'left' else len(self)
        q_values = []

        if side == 'left':
            for q in self.q_scores:
                q_values.append(q)
                if median(q_values[-window:]) < threshold:
                    pos = len(q_values) - 1
                else:
                    break

        elif side == 'right':
            for q in reversed(self.q_scores):
                q_values.append(q)
                if median(q_values[-window:]) < threshold:
                    pos = len(self) - len(q_values) + 1
                else:
                    break

        return pos

    def quality_trim_right(self, threshold=10):
        pos = 0
        for i in range(len(self.q_scores) - 5, 0, -5):
            mu = mean(self.q_scores)


class VCF(File):
    extensions = ("vcf", "cf",)

    def __init__(self, path):
        super().__init__(path=path, file_type="vcf")
        self.header = f">{self.filename.split('_vs_')[0]}"

    def to_fasta(self, matrix, lock=None, threshold=0.75):

        print(f"    ... adding {self.filename} to {matrix}")
        sequence = []
        with open(self.path) as fh:
            for line in fh:
                # CREATE SEQUENCE AS A LIST OF "-" THE SIZE OF THE REFERENCE GENOME
                if line.startswith("##contig"):
                    genome_length = int(line.split("length=")[1].split(">")[0])
                    sequence = ['-'] * genome_length
                # EXTRACT POSITION DATA FROM VCF LINES
                elif line[0] is not '#' and 'INDEL' not in line and len(line.split('\t')) == 10:
                    chrom, pos, id, ref, alt, qual, filter, info, format, etc = line.split('\t')
                    dp4 = line.split('DP4=')[1].split(';')[0].split(',')
                    reference = {'base': ref, 'depth': int(dp4[0]) + int(dp4[1])}
                    alternate = {'base': alt, 'depth': int(dp4[2]) + int(dp4[3])}
                    position = int(pos) - 1
                    # FILTER POSITIONS AND ASSIGN BASES TO SEQUENCE LIST
                    if reference['depth'] >= 4 and alternate['depth'] <= reference['depth'] * threshold:
                        base = reference['base']
                        # sequence[position] = reference['base'] if len(reference['base']) == 1 else 'N'
                    elif alternate['depth'] >= 4 and reference['depth'] <= alternate['depth'] * threshold:
                        base = alternate['base']
                        # if len(alternate['base']) == 1 and alternate['base'] is not '.':
                        #     sequence[position] = alternate['base']
                        # elif len(alternate['base']) > 1:
                        #     sequence[position] = 'N'
                    else:
                        base = '-'

                    if len(base) == 1 and base in "ATCG-":
                        sequence[position] = base
                    elif len(base) > 1:
                        sequence[position] = 'N'


        # RETURN TWO-LINE FASTA AS STRING
        if lock:
            lock.acquire() if lock else None
            with open(matrix, 'a+') as file:
                print(self.header + '\n' + ''.join(sequence), file=file)
            lock.release()
        else:
            with open(matrix, 'a+') as file:
                print(self.header + '\n' + ''.join(sequence), file=file)

    @classmethod
    @timing
    def dir_to_matrix(cls, directory, matrix):
        CPU_COUNT = mp.cpu_count()
        lock = mp.Lock()
        files = VCF.get_all(directory)
        all_workers = [mp.Process(target=vcf.to_fasta, args=(matrix, lock)) for vcf in files]
        print(f"==> Created worker pool of {len(all_workers)} workers")
        for i in range(0, len(all_workers), CPU_COUNT):
            workers = [mp.Process(target=vcf.to_fasta, args=(matrix, lock)) for vcf in files[i: i + CPU_COUNT]]
            print(f"==> Started {len(workers)} workers ({i}-{i + CPU_COUNT - 1} of {len(all_workers) - 1})")
            [worker.start() for worker in workers]
            [worker.join() for worker in workers]
            print("==> JOBS COMPLETE!")


        # print(f"==> Created {len(workers)} workers ")
        # [worker.start() for worker in workers]
        # print(f"==> Started workers")
        # [worker.join() for worker in workers]
        # print(f"==> Joined workers")
        # print("JOBS COMPLETE!")


        # with open(matrix, 'a+') as fh:
        #     print(self.header, file=matrix)
        #     print(self.sequence, file=matrix)


# DEERS 800-538-9552for f