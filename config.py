import os
import logging

from filehandler import Dir


class Config:

    # DIRECTORY STRUCTURE
    BASE_DIR = Dir('/Strong/proj/.data/Project_NTM')
    REFERENCES = Dir(BASE_DIR.join("lib", "reference_genomes"))

    DATA = BASE_DIR.make_subdir("data")
    TRIM_DIR = DATA.make_subdir("trimmed_reads")
    ASSEMBLY_DIR = DATA.make_subdir("assemblies")
    ANNOTATION_DIR = DATA.make_subdir("annotations")
    MAP_DIR = DATA.make_subdir("mapped_reads")

    # NAMING CONVENTIONS
    REPORTS = "reports"
    ARCHIVE = "archive"
    FAILED = "failed"
    LOGS = "logs"

    # REPORTS
    TRIM_STATS = "trim_stats.csv"
    ASSEMBLY_STATS = "assembly_stats.csv"

    # LOGGING
    LOG_LEVEL = logging.DEBUG
    LOG_FILE = BASE_DIR.join('log', 'pipeline.log')

    # REFERENCES
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
        'NFARC': 'NFARC.NCTC3000.fasta'}

    SPECIES_GROUPS = {
        'MAC': ['MAV', 'MCHIM', 'MINT']
    }

    @classmethod
    def declare_globals(cls):
        global BASEDIR, REFERENCES, REFERENCE_GENOMES, SPECIES_GROUPS

