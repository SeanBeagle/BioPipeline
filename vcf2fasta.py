"""vcf2fasta.py
AUTHOR : Sean Beagle
ORG    : http://StrongLab.org
NOTES  : Requires python 3.7 or higher

USAGE:
    # To add all files in current working directory
        python vcf2fasta.py MycobacteriumAviumMatrix.fa
    # To add all files in specified directory
        python vcf2fasta.py MycobacteriumAviumMatrix.fa /path/to/vcf_files
    # To create fasta from single vcf file
        python vcf2fasta.py MAV2201_aligned.fa MAV23201_vs_MavReference.vcf
"""

import multiprocessing as mp
import logging
import argparse
import os
from time import time
from functools import wraps


# DEFAULT PARAMETERS:
THRESHOLD = 0.75                 # ratio of call depth to total depth
MIN_DEPTH = 4                    # minimum depth of call
HEADER_DELIMITER = '_vs_'        # if "_vs_" and file name is "Isolate1_vs_Reference.vcf" then header is ">Isolate1"
USE_AMBIGUOUS_BASE = True        # if True and call is "A,T,C" then base is "N" else "-"
VCF_EXTENSION = (".vcf", ".cf")  # will find files with matching extensions (string or tuple of strings)


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


def logger(file, level="DEBUG", stdout=True):
    level = logging.getLevelName(level)
    # CONFIGURE LOGGING
    file = file
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    formatter = logging.Formatter(f'%(asctime)s:%(levelname)s:%(message)s', '%Y-%m-%d %H:%M')
    # SAVE LOGS TO: ``LOG_FILE``
    file_handler = logging.FileHandler(file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # SEND LOGS TO: ``sys.stdout``
    if stdout:
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setFormatter(formatter)
        logger.addHandler(consoleHandler)
    return logger


def get_vcf_files(path):
    """Return list of absolute filepaths for files ending with .cf and .vcf
    """
    try:
        return [
            os.path.join(path, file) for file in os.listdir(path)
            if os.path.isfile(os.path.join(path, file)) and file.endswith(VCF_EXTENSION)
        ]
    except FileNotFoundError:
        exit(f"{path} is not a valid directory")


def vcf_to_fasta(vcf, fasta, lock):
    """Return string of two line FastA
    """
    try:
        with open(vcf) as fin:
            header = f">{os.path.basename(vcf).split(HEADER_DELIMITER)[0]}"
            sequence = []
            for line in fin:
                # CREATE SEQUENCE AS A LIST OF "-" THE SIZE OF THE REFERENCE GENOME
                if line.startswith("##contig"):
                    genome_length = int(line.split("length=")[1].split(">")[0])
                    sequence = ['-'] * genome_length
                # EXTRACT POSITION DATA FROM VCF LINES
                elif line[0] != '#' and 'INDEL' not in line and len(line.split('\t')) == 10:
                    chrom, pos, id, ref, alt, qual, filter, info, format, etc = line.split('\t')
                    dp4 = line.split('DP4=')[1].split(';')[0].split(',')
                    ref_depth = int(dp4[0]) + int(dp4[1])
                    alt_depth = int(dp4[2]) + int(dp4[3])
                    total_depth = ref_depth + alt_depth
                    position = int(pos) - 1
                    # FILTER POSITIONS
                    if ref_depth >= MIN_DEPTH and ref_depth >= total_depth * THRESHOLD:
                        base = ref
                    elif alt_depth >= MIN_DEPTH and alt_depth >= total_depth * THRESHOLD:
                        base = alt
                    else:
                        base = '-'
                    # ASSIGN BASES TO SEQUENCE
                    if len(base) == 1 and base in "ATCG-":
                        sequence[position] = base
                    elif len(base) > 1 and USE_AMBIGUOUS_BASE:
                        sequence[position] = 'N'
        # WRITE TO FASTA
        lock.acquire()
        with open(fasta, 'a+') as file:
            print(header + '\n' + ''.join(sequence), file=file)
        lock.release()
        print(f"\t[+] {header[1:]} ({len(sequence)}bp)")

    except Exception as e:
        print(e)


def verify_matrix(path):
    headers = []
    lengths = []

    return len(set(headers)) == len(headers) & len(set(lengths)) == 1


@timing
def main(input, output, log, threads):
    log = log or output + '.log'
    manager = mp.Manager()
    lock = manager.Lock()
    if os.path.isdir(input):
        files = get_vcf_files(input)
        if len(files) > 0:
            pool = mp.Pool(threads)
            print(f"Creating matrix {output} from {len(files)} files in {input}")
            for file in files:
                pool.apply_async(vcf_to_fasta, (file, output, lock))
            pool.close()
            pool.join()
        else:
            print(f"No VCF files found in {input}")
    elif os.path.isfile(input):
        print(f"Creating fasta {output} from {input}")
        vcf_to_fasta(input, output, lock)
    else:
        raise FileNotFoundError(input)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert VCF files to FastA")
    parser.add_argument('--output', help='output file for fasta sequence(s)',
                        required=True, type=os.path.abspath
                        )
    parser.add_argument('--input', help='vcf file or directory of vcf files',
                        default=os.getcwd(), type=os.path.abspath
                        )
    parser.add_argument('--log', help='log file',
                        default=None
                        )
    parser.add_argument('--threshold', help='minimum ratio of call depth/total depth',
                        default=0.75, type=float
                        )
    parser.add_argument('--mindepth',
                        default=4, type=int, help='minimum depth of call'
                        )
    parser.add_argument('-n', action='store_true',
                        help='Use "N" instead of "-" for ambiguous base calls'
                        )
    parser.add_argument('--threads', help=f'number of parallel threads (max={mp.cpu_count()})',
                        default=mp.cpu_count(), type=(lambda x: min(int(x), mp.cpu_count()) or mp.cpu_count()),
                        )
    args = parser.parse_args()

    print(args.output)
    print(args.input)

    # main(input=args.input, output=args.output, log=args.log, threads=args.threads)

