#!/Strong/proj/.data/Morty/.config/venv/bin/python

# PYTHON STANDARD LIBRARY
import os
import sys
import subprocess
import platform
import logging
from types import SimpleNamespace as namespace
# EXTERNAL LIBRARY
import pandas as pd
# PROJECT LIBRARY
from filehandler import Fastq, Fasta, File, Dir
from clustertools import Module, LSF


def configure():
    global BASEDIR, REFERENCES, REFERENCE_GENOMES, SPECIES_GROUPS
    BASEDIR = Dir('/Strong/proj/.data/ProjectNTM')
    REFERENCES = Dir(BASEDIR.join("lib", "reference_genomes"))
    REFERENCE_GENOMES = {
        # Clinical
        'MAB': 'MAB.ATCC19977.fasta',
        'MBOL': 'MAB.ATCC19977.fasta',
        'MAV': 'MAV.HOM.H87.fasta',
        'MMAS': 'MMAS.BRAPA42FWDG01.fasta',
        'MCHIM': 'MCHIM.CDC2015-22-71.fasta',
        'MINT': 'MCHIM.CDC2015-22-71.fasta',
        'MCHE': 'MCHE.ATCC19237.fasta',
        'MTB': 'MTB.H37RV.fasta',
        # Environmental
        'MAROS': 'MAROS.DSM45069.fasta',
        'MASIA': 'MASIA.DSM44297.fasta',
        'MBOUCH': 'MBOUCH.DSM45439.fasta',
        'MBOV': 'MBOV.AF2122.fasta',
        'MCANE': 'MCANE.CIPT140070017.fasta',
        'MCHUB': 'MCHUB.NBB4.fasta',
        'MCOLOM': 'MCOLOM.CECT3035.fasta',
        'MELE': 'MELE.DSM44368.fasta',
        'MFORT': 'MFORT.CT6.fasta',
        'MFRANK': 'MFRANK.DSM45524.fasta',
        'MGILV': 'MGILV.SPYR1.fasta',
        'MGORD': 'MGORD.DSM44160.fasta',
        'MHAEM': 'MHAEM.DSM44634.fasta',
        'MIMMU': 'MIMMU.CCUG47286T.fasta',
        'MINDP': 'MINDP.MTCC9506.fasta',
        'MIRAN': 'MIRAN.DSM45541.fasta',
        'MKAN': 'MKAN.ATCC12478.fasta',
        'MKUB': 'MKUB.CIP106428.fasta',
        'MLENT': 'MLENT.CSURP1491.fasta',
        'MLEPR': 'MLEPR.TN.fasta',
        'MLIFL': 'MLIFL.128FXT.fasta',
        'MMANT': 'MMANT.DSM45255.fasta',
        'MMARI': 'MMARI.M.fasta',
        'MMARS': 'MMARS.DSM45437.fasta',
        'MMUCO': 'MMUCO.CSURP2099.fasta',
        'MNEOA': 'MNEOA.VKMAC-1815D.fasta',
        'MPORC': 'MPORC.CSURP1564.fasta',
        'MRHOD': 'MRHOD.NBB3.fasta',
        'MSALM': 'MSALM.D16Q15.fasta',
        'MSENE': 'MSENE.NCTC4524.fasta',
        'MSIMI': 'MSIMI.ATCC25275.fasta',
        'MSMEG': 'MSMEG.MC2155.fasta',
        'MTERR': 'MTERR.NCTC10856.fasta',
        'MTIM': 'MTIM.CCUG56329.fasta',
        'MTRIP': 'MTRIP.DSM44626.fasta',
        'MULCE': 'MULCE.AGY99.fasta',
        'MVANB': 'MVANB.PYR-1.fasta',
        'MVUL': 'MVUL.DSM45247T.fasta',
        'MXENO': 'MXENO.RIVM700367.fasta',
        'MYONG': 'MYONG.05-1390.fasta',
        'NFARC': 'NFARC.NCTC3000.fasta'
    }
    SPECIES_GROUPS = {
        'MAC': ['MAV', 'MCHIM', 'MINT', 'MTIM', 'MBOUCH', 'MMARS'],
        'MAB': ['MAB', 'MBOL', 'MMAS']
    }


class References:
    def __init__(self):
        directory = Dir(BASEDIR.join("lib", "reference_genomes"))
        references = Fasta.get_all(directory)


class Isolate:
    """Create isolate instance to manage metadata and files.
    """
    def __init__(self, pair1, pair2, log_level='DEBUG', *args, **kwargs):

        self.files = namespace()
        self.files.raw_pair = Fastq(pair1), Fastq(pair2)
        self.run = self.files.raw_pair[0].dir.dirname
        self.name = self.files.raw_pair[0].sample_name
        self.taxon = 'UNKNOWN'

        self._logger = self.configure_logging(level=log_level)
        self.log("Creating Isolate Instance", lvl="INFO")

        self.build_filesystem()

    def __repr__(self):
        return self.name

    def log(self, msg, lvl='DEBUG', *args, **kwargs):
        lvl = logging.getLevelName(lvl)
        self._logger.log(lvl, msg, *args, **kwargs)

    def build_filesystem(self):
        # DIRECTORY STRUCTURE
        global SCRATCH, REFERENCES, DATA, TRIM_DIR, ASSEMBLY_DIR, ANNOTATION_DIR, MAP_DIR
        SCRATCH = BASEDIR.make_subdir('tmp', 'scratch', self.run, self.name)
        # SCRATCH = Dir.make('/scratch/' + os.environ['USER'])
        DATA = BASEDIR.make_subdir("data")
        TRIM_DIR = DATA.make_subdir('trimmed_reads')
        ASSEMBLY_DIR = DATA.make_subdir('assemblies')
        ANNOTATION_DIR = DATA.make_subdir('annotations')
        MAP_DIR = DATA.make_subdir('mapped_reads')

    def configure_logging(self, level='DEBUG', to_stdout=True, to_file=True, to_database=False):
        node = platform.node()

        # CONFIGURE LOGGING
        logger = logging.getLogger(__name__)
        logger.setLevel(level)
        formatter = logging.Formatter(f'%(asctime)s:%(levelname)s:{node}:{self}:%(message)s', '%Y-%m-%d.%H%M')

        if to_file:
            # SAVE LOGS TO TEXT FILE
            file = BASEDIR.make_subdir('logs', self.run).join(f'{self}.log')
            file_handler = logging.FileHandler(file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        if to_stdout:
            # RETURN LOGS TO STANDARD OUTPUT
            consoleHandler = logging.StreamHandler(sys.stdout)
            consoleHandler.setFormatter(formatter)
            logger.addHandler(consoleHandler)
        if to_database:
            # SAVE LOGS TO DATABASE
            # TODO: CREATE LOGGING TO DATABASE OPTION
            pass
        return logger


def update_taxa(species_threshold=0.97, genus_threshold=0.80):
    logger = generic_logger('update_taxa.csv')

    run = Dir().dirname
    trim_dir = Dir(BASEDIR.join('data', 'trimmed_reads', run))
    assembly_dir = Dir(BASEDIR.join('data', 'assemblies', run))
    ani_dir = Dir(assembly_dir.join('ANI'))

    trimmed_reads = trim_dir.files(endswith="fq", dataframe=True)
    assemblies = assembly_dir.files(endswith='fna', dataframe=True)
    ani = ani_dir.files(endswith='.csv')

    for file in ani:
        try:
            df = pd.read_csv(file.path).sort_values('ani')
            sample_name = df.iloc[0].sample
            trim1 = trimmed_reads[(trimmed_reads.filename.str.contains(sample_name)) & (trimmed_reads.filename.str.contains('_R1'))].iloc[0].path
            trim2 = trimmed_reads[(trimmed_reads.filename.str.contains(sample_name)) & (trimmed_reads.filename.str.contains('_R2'))].iloc[0].path
            assembly = assemblies[(assemblies.filename.str.contains(sample_name))].iloc[0].path

            trim1 = File(trim1)
            trim2 = File(trim2)
            assembly = File(assembly)

            # ASSIGN TAXON
            taxon = 'UNKNOWN'
            possible_species = df[(df.ani >= species_threshold)]
            if len(possible_species) > 0:
                taxon = possible_species.iloc[0].taxon
            elif len(df[(df.ani >= genus_threshold)]) > 0:
                taxon = 'NTM'

            trim1_filename = trim1.filename
            trim2_filename = trim2.filename
            assembly_filename = assembly.filename

            trim1.rename(f'{sample_name}_{taxon}.fq.gz')
            trim2.rename(f'{sample_name}_{taxon}.fq.gz')
            assembly.rename(f'{sample_name}_{taxon}_000.fna')

            logger.info(f'renamed {trim1_filename} to {trim1.filename}')
            logger.info(f'renamed {trim2_filename} to {trim1.filename}')
            logger.info(f'renamed {assembly_filename} to {assembly.filename}')

        except Exception as e:
            logger.warning(e)





def process_directory(directory=os.getcwd()):
    """Run pipeline on all fastq pairs in directory"""
    [LSF.bsub(f"{sys.argv[0]} pipeline {pair.pair1} {pair.pair2}") for pair in Fastq.get_pairs(Dir(directory))]


def combine_lanes(directory=os.getcwd()):
    """
    REPORT:
        L001_R1 counts
        L002_R1 counts


    CF RULES 1e9 in one lane

    """


    directory = Dir(directory)
    temp = BASEDIR.make_subdir('tmp', directory.dirname)
    logger = generic_logger(temp.join('combine_lanes.log'))
    pairs = Fastq.get_pairs(directory)

    logger.info(f"Combining lanes of {len(pairs)} pairs in {directory}")
    for pair in pairs:
        try:
            if pair.pair1.lane == "L001" and os.path.isfile(pair.pair1.path.replace("L001", "L002")):
                # READ 1
                L001_R1 = pair.pair1
                L002_R1 = Fastq(L001_R1.path.replace("L001", "L002"))
                LCAT_R1 = temp.join(L001_R1.filename.replace("L001", "LCAT"))
                # READ 2
                L001_R2 = pair.pair2
                L002_R2 = Fastq(L001_R2.path.replace("L001", "L002"))
                LCAT_R2 = temp.join(L001_R2.filename.replace("L001", "LCAT"))
                try:
                    # COMBINE LANE 1
                    logger.info(f"Concatenating {L001_R1.filename} and {L002_R1.filename} as {LCAT_R1}")
                    subprocess.Popen(f'cat {L001_R1} {L002_R1} > {LCAT_R1}', shell=True)
                    logger.info("...PASS")
                except Exception as e:
                    logger.warning(f"...FAIL: {e}")
                try:
                    # COMBINE LANE 2
                    logger.info(f"Concatenating {L001_R2.filename} and {L002_R2.filename} as {LCAT_R2}")
                    subprocess.Popen(f'cat {L001_R2} {L002_R2} > {LCAT_R2}', shell=True)
                    logger.info("...PASS")
                except Exception as e:
                    logger.info(f"...FAIL: {e}")
        except Exception as e:
            logger.warning(f'Combining Failed: {e}')


def pipeline(pair1, pair2):
    isolate = Isolate(pair1, pair2, log_level='DEBUG')
    isolate.log('Starting Pipeline')

    try:
        trim(isolate)
    except RuntimeError as e:
        isolate.log(f'Trimming Failed: {e}', lvl='WARNING')
        exit('Trimming Failed')

    try:
        assemble(isolate)
    except RuntimeError as e:
        isolate.log(f'Assembling Failed: {e}', lvl='WARNING')
        exit('Assembly Failed')

    try:
        identify(isolate)
    except Exception as e:
        isolate.log(f'ANI Failed {e}', lvl='WARNING')
        exit('ANI Failed')

    try:
        rename_files(isolate)
    except Exception as e:
        isolate.log(f'Renaming Failed {e}', lvl='WARNING')

    # try:
    #     LOG.debug(f'{isolate.name}: Reordering Contigs')
    #     reorder(isolate)
    # except RuntimeError as e:
    #     LOG.warning(f'{isolate.name}: Reordering Failed')


def trim(isolate, **kwargs):
    """Trim raw fastq pair using skewer."""
    isolate.log('Calling Trimmer')

    try:
        pair1, pair2 = isolate.files.raw_pair
    except AttributeError as e:
        raise RuntimeError(e)

    # LOAD MODULES
    Module.load('skewer/0.2.2')
    # SET DEFAULT PARAMETERS
    adapter = kwargs.setdefault('adapter', 'CTGTCTCTTATACACATCT')
    quality = kwargs.setdefault('quality', 20)
    minlen = kwargs.setdefault('minlen', 40)
    sample = pair1.filename.split('_')[0]
    # FILES AND DIRS
    scratch = SCRATCH.make_subdir(f"{pair1.inode}")
    trim_dir = TRIM_DIR.make_subdir(pair1.dir.dirname)
    log_dir = trim_dir.make_subdir('logs')
    trim1 = pair1.filename.replace('_001', '_002')
    trim2 = pair2.filename.replace('_001', '_002')

    # TRIM READS
    cmd = f"skewer -x {adapter} -q {quality} -l {minlen} -m pe -z -o {scratch}/ {pair1} {pair2}"
    isolate.log(f'Running: {cmd}', lvl='INFO')
    try:
        job = subprocess.run(cmd.split(), capture_output=True)
        if job.returncode == 0:
            trim1 = Fastq(scratch.join("trimmed-pair1.fastq.gz")).rename(trim1).move(trim_dir)
            trim2 = Fastq(scratch.join("trimmed-pair2.fastq.gz")).rename(trim2).move(trim_dir)
            File(scratch.join("trimmed.log")).rename(f"{isolate}_skewer.log").move(log_dir)
            isolate.files.trimmed_pair = trim1, trim2
            return trim1, trim2
        elif job.returncode == 1:
            isolate.log(job.stderr, lvl='WARNING')
            return RuntimeError(cmd, job.stderr)
    except Exception as e:
        isolate.log(e, lvl='WARNING')
        return RuntimeError(cmd, e)


def assemble(isolate, **kwargs):
    isolate.log('Assembling Isolate')

    try:
        pair1, pair2 = isolate.files.trimmed_pair
    except AttributeError as e:
        raise RuntimeError(e)

    # LOAD ENVIRONMENT MODULES
    Module.load('spades/3.11.1', 'unicycler/git', 'racon/git', 'bowtie2/2.3.2',
                'samtools/1.5', 'pilon/testing', 'bcftools/1.3.1', 'blast+/2.6.0')
    # SET DEFAULT PARAMETERS
    threads = kwargs.setdefault('threads', 20)
    minlen = kwargs.setdefault('minlen', 1000)
    mode = kwargs.setdefault('mode', 'normal')
    verbosity = kwargs.setdefault('verbosity', 3)

    # FILES/DIRS
    assembly_dir = ASSEMBLY_DIR.make_subdir(pair1.dir.dirname, pair1.sample_name)
    archive_dir = ASSEMBLY_DIR.make_subdir(pair1.dir.dirname, 'archive', "unicycler")

    # ASSEMBLE READS
    cmd = f"unicycler-runner.py -1 {pair1} -2 {pair2} -o {assembly_dir} --min_fasta_length {minlen} --keep 2 --mode {mode} -t {threads} --vcf --verbosity {verbosity}"
    isolate.log(f'Running: {cmd}', lvl='INFO')
    try:
        job = subprocess.run(cmd.split(), capture_output=True)
        if job.returncode == 0:
            default_assembly = Fasta(assembly_dir.join('assembly.fasta'))
            new_name = f'{isolate}_000_000.fna'
            assembly = default_assembly.copy(new_name).move(default_assembly.dir.parent.path)
            assembly_dir.move(archive_dir.path)
            isolate.files.assembly = assembly
            return assembly
        elif job.returncode == 1:
            isolate.log(job.stderr, lvl='WARNING')
            return RuntimeError(cmd, job.stderr)
    except Exception as e:
        isolate.log(e, lvl='WARNING')
        return RuntimeError(cmd, e)


def identify(isolate, delimiter="_", species_threshold=0.97, genus_threshold=0.80):
    isolate.log('Identifying Isolate')
    assembly = isolate.files.assembly

    """Blast fasta to reference genomes and store values in db"""
    ani_script = "/Strong/proj/.data/Morty/.config/software/ani-script/ANI.pl"
    blastall = "/software/cgeh/blast/2.2.22/bin/blastall"
    formatdb = "/software/cgeh/blast/2.2.22/bin/formatdb"

    # COMPARE TO REFERENCES
    try:
        matches = []
        references = Fasta.get_all(REFERENCES)
        isolate.log(f"{isolate}: IDENTIFYING TAXON USING {len(references)} REFERENCES", lvl='INFO')
        for reference in references:
            ref_id = reference.filename.split('.')[0]
            scratch = SCRATCH.make_subdir('ani', isolate.name, f"{isolate}_vs_{ref_id}")
            command = f"perl {ani_script} -bl {blastall} -fd {formatdb} -qr {assembly} -sb {reference} -od {scratch}"
            output, error = subprocess.Popen(command.split(), stdout=subprocess.PIPE).communicate()

            try:
                ani = float(output) / 100
            except (ValueError, TypeError):
                ani = 0
            finally:
                record = {
                    'ani': ani,
                    'sample': isolate.name,
                    'reference': ref_id,
                    'taxon': ref_id.split('_')[0].split('-')[0]}
                matches.append(record)

        # WRITE CSV
        ani_csv = Dir.make(assembly.dir.join("ANI")).join(f"{isolate}_ANI.csv")
        df = pd.DataFrame.from_records(matches)
        df = df.sort_values('ani', ascending=False)
        df.to_csv(ani_csv, index=False)
        isolate.files.ani = File(ani_csv)

        # ASSIGN TAXON
        taxon = 'UNKNOWN'
        possible_species = df[(df.ani >= species_threshold)]
        if len(possible_species) > 0:
            taxon = possible_species.iloc[0].taxon
        elif len(df[(df.ani >= genus_threshold)]) > 0:
            taxon = 'NTM'

        isolate.taxon = taxon
        isolate.log(f"taxon={isolate.taxon}", lvl='INFO')
        return taxon

    except Exception as e:
        isolate.log(f"Identification failed: {e}", lvl='WARNING')
        subprocess.call(f"rm error.log formatdb.log".split())


def rename_files(isolate):
    isolate.log('Renaming Files')
    try:
        isolate.files.trimmed_pair[0].rename(f'{isolate}_{isolate.taxon}_R1.fq.gz')
        isolate.files.trimmed_pair[1].rename(f'{isolate}_{isolate.taxon}_R2.fq.gz')
    except AttributeError as e:
        raise RuntimeError(e)

    try:
        isolate.files.assembly.rename(f'{isolate}_{isolate.taxon}_000.fna')
    except AttributeError as e:
        raise RuntimeError(e)


# def reorder(isolate):
#     # LOAD MODULES
#     Module.load('mauve/2015-02-13')
#     mauve = "/software/cgeh/mauve/2015-02-13/install/Mauve.jar"
#
#     # INPUT
#     assembly = isolate.files.assembly
#
#     # FILES/DIRS
#     archive_dir = Dir.make(assembly.dir.join(ARCHIVE, "unordered"))
#
#     # FILES/DIRS
#     reference = Fasta(REFERENCES.join(REFERENCE_GENOMES.get(isolate.taxon)))
#     out = SCRATCH.make_subdir("mauve", isolate.name)
#
#     # REORDER CONTIGS
#     try:
#         LOG.info(f"{isolate.name} : REORDERING USING REFERENCE : {reference.filename}")
#         cmd = f"java -Xmx500m -cp {mauve} org.gel.mauve.contigs.ContigOrderer -output {out} -ref {reference} -draft {assembly}"
#         process = subprocess.run(cmd.split(), capture_output=True)
#         if process.returncode == 0:
#             LOG.info(f"{isolate.name} : REORDERING COMPLETE")
#         if process.returncode == 1:
#             LOG.info(f"{isolate.name} : REORDERING FAILED: {process.stderr}")
#             exit(process.stderr)
#     except Exception as e:
#         LOG.warning(f"{isolate.name} REORDERING FAILED: {e}")
#         raise RuntimeError(e)
#
#
#     # CLEAN UP
#     LOG.info(f"{isolate.name} : CLEANING UP")
#     final_dir = Dir(out.join(f"alignment{len(out.children)}"))
#     final_name = fa.filename.replace(".000.", ".finalAssembly.")
#     final_assembly = Fasta(final_dir.join(assembly.filename)).rename(final_name).move(assembly.dir.path)
#     assembly.move(archive_dir)



# def quick_id(fasta, out="ani.csv"):
#     """Blast fasta to reference genomes and store values in db"""
#     ani_script = "/Strong/proj/.data/Morty/.config/software/ani-script/ANI.pl"
#     blastall = "/software/cgeh/blast/2.2.22/bin/blastall"
#     formatdb = "/software/cgeh/blast/2.2.22/bin/formatdb"
#
#     # INPUT
#     fasta = Fasta(fasta)
#     references = Fasta.get_all(REFERENCES.path)
#
#     # COMPARE TO REFERENCES
#     matches = 0
#     try:
#         for reference in references:
#             ref_id = '.'.join(reference.filename.split('.')[:-1])
#             taxon = ref_id.split(".")[0]
#             scratch = SCRATCH.make_subdir('ani', fasta.filename, f"{fasta.filename}_vs_{ref_id}")
#             command = [
#                 'perl', ani_script, '-bl', blastall, '-fd', formatdb, '-qr', fasta.path,
#                 '-sb', reference.path, '-od', scratch.path]
#             (output, error) = subprocess.Popen(command, stdout=subprocess.PIPE).communicate()
#             try:
#                 output = float(output) / 100
#                 ani = float(format(output, '.7f'))
#             except ValueError:
#                 ani = 0.0
#             matches += 1
#
#             print(f"{fasta.filename},{taxon},{ani},{ref_id}")
#             # write to csv
#             with open(out, 'at+') as fh:
#                 fh.write(f"\n{fasta.filename},{taxon},{ani},{ref_id}")
#             subprocess.call(f"rm error.log formatdb.log".split())
#         print(f"{fasta.filename} : IDENTIFICATION COMPLETE: {matches} MATCHES")
#     except Exception as e:
#         print(f"Whoops! Something went wrong...\n{e}")

# class Reorder:
#     """Reorder contigs of WGS Assembly using a reference template."""
#     def __init__(self, fasta):
#
#         # PARAMETERS
#         self.mauve = "/software/cgeh/mauve/2015-02-13/install/Mauve.jar"
#
#         # NEXT
#         pipe = Annotate
#
#         # CALL SELF
#         self.__call__(fasta, pipe)
#
#     def __call__(self, fasta, pipe):
#
#         # LOAD MODULES
#         load_module('mauve/2015-02-13')
#
#         # INPUT
#         assembly = Fasta(fasta)
#
#         # FILES/DIRS
#         archive_dir = Dir.make(assembly.dir.join(ARCHIVE, "unordered"))
#
#         # METADATA
#         sample, taxon, process_level, extension = assembly.filename.split(".")
#         project = assembly.dir.dirname
#
#         # FILES/DIRS
#         reference = Fasta(REFERENCES.join(REFERENCE_GENOMES.get(taxon)))
#         out = SCRATCH.make_subdir("mauve", sample)
#
#         # REORDER CONTIGS
#         try:
#             LOG.info(f"{sample} : REORDERING USING REFERENCE : {reference.filename}")
#             command = [
#                 "java", "-Xmx500m", "-cp", self.mauve, "org.gel.mauve.contigs.ContigOrderer",
#                 "-output", out.path, "-ref", reference.path, "-draft", assembly.path]
#             process = subprocess.run(command, capture_output=True)
#             if process.returncode == 0:
#                 LOG.info(f"{sample} : REORDERING COMPLETE")
#             if process.returncode == 1:
#                 LOG.info(f"{sample} : REORDERING FAILED: {process.stderr}")
#                 exit(process.stderr)
#         except Exception as e:
#             LOG.warning(f"{sample} REORDERING FAILED: {e}")
#
#         # CLEAN UP
#         LOG.info(f"{sample} : CLEANING UP")
#         final_dir = Dir(out.join(f"alignment{len(out.children)}"))
#         final_name = assembly.filename.replace(".000.", ".finalAssembly.")
#         final_assembly = Fasta(final_dir.join(assembly.filename)).rename(final_name).move(assembly.dir.path)
#         assembly.move(archive_dir)
#
#         # NEXT
#         pipe(final_assembly.path)

def build_filesystem():
    # DIRECTORY STRUCTURE
    global SCRATCH, BASEDIR, REFERENCES, DATA, TRIM_DIR, ASSEMBLY_DIR, ANNOTATION_DIR, MAP_DIR
    SCRATCH = Dir.make('/scratch/' + os.environ['USER'])
    BASEDIR = Dir('/Strong/proj/.data/Project_NTM')
    REFERENCES = Dir(BASEDIR.join("lib", "reference_genomes"))
    DATA = BASEDIR.make_subdir("data")
    # RAW_DIR = DATA.make_subdir("00_raw")
    TRIM_DIR = DATA.make_subdir('trimmed_reads')
    ASSEMBLY_DIR = DATA.make_subdir('assemblies')
    ANNOTATION_DIR = DATA.make_subdir('annotations')
    MAP_DIR = DATA.make_subdir('mapped_reads')


def generic_logger(file, level="DEBUG", stdout=True):
    level = logging.getLevelName(level)
    # CONFIGURE LOGGING
    file = file
    node = platform.node()
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    formatter = logging.Formatter(f'%(asctime)s:%(levelname)s:{node}:%(message)s', '%Y-%m-%d %H:%M')
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


if __name__ == '__main__':
    configure()

    if len(sys.argv) == 2 and sys.argv[1] == 'process_directory':
        process_directory()
    elif len(sys.argv) == 4 and sys.argv[1] == 'pipeline':
        pipeline(pair1=sys.argv[2], pair2=sys.argv[3])
    elif len(sys.argv) == 2 and sys.argv[1] == 'combine_lanes':
        combine_lanes()
    elif len(sys.argv) == 2 and sys.argv[1] == 'update_taxa':
        update_taxa()
    # elif len(sys.argv) == 3 and sys.argv[1] == 'identify':
    #     identify(assembly=Fasta(sys.argv[2]))


    # # DIRECTORY STRUCTURE
    # SCRATCH = Dir.make('/scratch/' + os.environ['USER'])
    # BASEDIR = Dir('/Strong/proj/.data/Project_NTM')
    # REFERENCES = Dir(BASEDIR.join("lib", "reference_genomes"))
    # DATA = BASEDIR.make_subdir("data")
    # # RAW_DIR = DATA.make_subdir("00_raw")
    # TRIM_DIR = DATA.make_subdir('trimmed_reads')
    # ASSEMBLY_DIR = DATA.make_subdir('assemblies')
    # ANNOTATION_DIR = DATA.make_subdir('annotations')
    # MAP_DIR = DATA.make_subdir('mapped_reads')
    #
    # # NAMING CONVENTIONS
    # REPORTS = "reports"
    # ARCHIVE = "archive"
    # FAILED = "failed"
    # LOGS = "logs"
    #
    # # REPORTS
    # TRIM_STATS = "trim_stats.csv"
    # ASSEMBLY_STATS = "assembly_stats.csv"
    #
    # # LOGGING
    # LOG_LEVEL = logging.DEBUG
    # LOG_FILE = BASEDIR.join('log', 'pipeline.log')

    # REFERENCES


    # # CONFIGURE LOGGING
    # LOG = logging.getLogger(__name__)
    # LOG.setLevel(LOG_LEVEL)
    # formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s', '%Y%m%d:%H%M')
    # # SAVE LOGS TO: ``LOG_FILE``
    # file_handler = logging.FileHandler(LOG_FILE)
    # file_handler.setFormatter(formatter)
    # LOG.addHandler(file_handler)
    # # SEND LOGS TO: ``sys.stdout``
    # consoleHandler = logging.StreamHandler(sys.stdout)
    # consoleHandler.setFormatter(formatter)
    # LOG.addHandler(consoleHandler)













