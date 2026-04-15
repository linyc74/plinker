import plinker
import argparse
from typing import Any


__VERSION__ = '1.0.0-beta'


PROG = "python plinker"
DESCRIPTION = f'GWAS pipeline based on PLINK 1.9\nAuthor: Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['--bfile'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'plink binary files prefix, expects .bed/.bim/.fam files',
        },
    },
    {
        'keys': ['--id-link-xslx'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to ID link Excel file with TPMI ID and UUID',
        },
    },
    {
        'keys': ['--phenotype-xslx'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to phenotype Excel file with UUID and phenotype columns',
        },
    },
    {
        'keys': ['--uuid-column'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'UUID column name in the ID link Excel file and phenotype Excel file',
        },
    },
    {
        'keys': ['--tpmi-id-column'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'TPMI ID column name in the ID link Excel file',
        },
    },
    {
        'keys': ['--phenotype-columns'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'comma-separated phenotype column names in the phenotype Excel file',
        },
    },
]
OPTIONAL = [
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'plinker_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['--minimum-minor-allele-frequency'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.01,
            'help': 'minimum minor allele frequency (MAF) (default: %(default)s)',
        },
    },
    {
        'keys': ['--maximum-per-variant-missing-genotype-rate'],
        'properties': {
            "type": float,
            'required': False,
            'default': 0.05,
            'help': 'maximum missing genotype rate per variant (default: %(default)s)',
        },
    },
    {
        'keys': ['--maximum-per-sample-missing-genotype-rate'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.01,
            'help': 'maximum missing genotype rate per sample (default: %(default)s)',
        },
    },
    {
        'keys': ['--hardy-weinberg-p-value-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 1e-6,
            'help': 'Hardy-Weinberg equilibrium p-value threshold (default: %(default)s)',
        },
    },
    {
        'keys': ['--association-p-value-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 1e-5,
            'help': 'association analysis p-value threshold (default: %(default)s)',
        },
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running PLINKer version {__VERSION__}\n', flush=True)
        plinker.main(
            bfile=args.bfile,
            id_link_xslx=args.id_link_xslx,
            phenotype_xslx=args.phenotype_xslx,
            uuid_column=args.uuid_column,
            tpmi_id_column=args.tpmi_id_column,
            phenotype_columns=args.phenotype_columns,
            minimum_minor_allele_frequency=args.minimum_minor_allele_frequency,
            maximum_per_variant_missing_genotype_rate=args.maximum_per_variant_missing_genotype_rate,
            maximum_per_sample_missing_genotype_rate=args.maximum_per_sample_missing_genotype_rate,
            hardy_weinberg_p_value_threshold=args.hardy_weinberg_p_value_threshold,
            association_p_value_threshold=args.association_p_value_threshold,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)

if __name__ == "__main__":
    EntryPoint().main()
