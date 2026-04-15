import os
import shutil
import pandas as pd
from copy import copy
from typing import List
from os.path import basename, dirname, join
from .utils import edit_fpath
from .template import Processor


PLINK_EXE = f'{dirname(dirname(__file__))}\\plink.exe'


class Plinker(Processor):

    bfile: str
    id_link_xslx: str
    phenotype_xslx: str
    uuid_column: str
    tpmi_id_column: str
    phenotype_columns: List[str]
    minimum_minor_allele_frequency: float
    maximum_per_variant_missing_genotype_rate: float
    maximum_per_sample_missing_genotype_rate: float
    hardy_weinberg_p_value_threshold: float

    def main(
            self,
            bfile: str,
            id_link_xslx: str,
            phenotype_xslx: str,
            uuid_column: str,
            tpmi_id_column: str,
            phenotype_columns: List[str],
            minimum_minor_allele_frequency: float,
            maximum_per_variant_missing_genotype_rate: float,
            maximum_per_sample_missing_genotype_rate: float,
            hardy_weinberg_p_value_threshold: float):

        self.bfile = bfile
        self.id_link_xslx = id_link_xslx
        self.phenotype_xslx = phenotype_xslx
        self.uuid_column = uuid_column
        self.tpmi_id_column = tpmi_id_column
        self.phenotype_columns = phenotype_columns
        self.minimum_minor_allele_frequency = minimum_minor_allele_frequency
        self.maximum_per_variant_missing_genotype_rate = maximum_per_variant_missing_genotype_rate
        self.maximum_per_sample_missing_genotype_rate = maximum_per_sample_missing_genotype_rate
        self.hardy_weinberg_p_value_threshold = hardy_weinberg_p_value_threshold

        for phenotype_column in self.phenotype_columns:
            settings = copy(self.settings)
            settings.outdir = join(self.outdir, phenotype_column)
            settings.workdir = join(self.workdir, phenotype_column)
            for d in [settings.workdir, settings.outdir]:
                os.makedirs(d, exist_ok=True)
            OnePhenotypePipeline(settings).main(
                bfile=self.bfile,
                id_link_xslx=self.id_link_xslx,
                phenotype_xslx=self.phenotype_xslx,
                uuid_column=self.uuid_column,
                tpmi_id_column=self.tpmi_id_column,
                phenotype_column=phenotype_column,
                minimum_minor_allele_frequency=self.minimum_minor_allele_frequency,
                maximum_per_variant_missing_genotype_rate=self.maximum_per_variant_missing_genotype_rate,
                maximum_per_sample_missing_genotype_rate=self.maximum_per_sample_missing_genotype_rate,
                hardy_weinberg_p_value_threshold=self.hardy_weinberg_p_value_threshold)


class OnePhenotypePipeline(Processor):
    
    bfile: str
    id_link_xslx: str
    phenotype_xslx: str
    uuid_column: str
    tpmi_id_column: str
    phenotype_column: str
    minimum_minor_allele_frequency: float
    maximum_per_variant_missing_genotype_rate: float
    maximum_per_sample_missing_genotype_rate: float
    hardy_weinberg_p_value_threshold: float

    keep_samples_phen: str

    def main(
            self,
            bfile: str,
            id_link_xslx: str,
            phenotype_xslx: str,
            uuid_column: str,
            tpmi_id_column: str,
            phenotype_column: str,
            minimum_minor_allele_frequency: float,
            maximum_per_variant_missing_genotype_rate: float,
            maximum_per_sample_missing_genotype_rate: float,
            hardy_weinberg_p_value_threshold: float):

        self.bfile = bfile
        self.id_link_xslx = id_link_xslx
        self.phenotype_xslx = phenotype_xslx
        self.uuid_column = uuid_column
        self.tpmi_id_column = tpmi_id_column
        self.phenotype_column = phenotype_column
        self.minimum_minor_allele_frequency = minimum_minor_allele_frequency
        self.maximum_per_variant_missing_genotype_rate = maximum_per_variant_missing_genotype_rate
        self.maximum_per_sample_missing_genotype_rate = maximum_per_sample_missing_genotype_rate
        self.hardy_weinberg_p_value_threshold = hardy_weinberg_p_value_threshold

        self.copy_and_clean_bfile()
        self.build_keep_samples_phen()
        self.plink_keep()
        self.plink_pheno()
        self.plink_qc()
        self.plink_assoc()
        self.sort_by_association_p_values()

    def copy_and_clean_bfile(self):
        for ext in ['.bed', '.bim', '.fam']:
            shutil.copy(f'{self.bfile}{ext}', f'{self.workdir}\\')
        self.bfile = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='',
            dstdir=self.workdir
        )
        df = read_fam(path=f'{self.bfile}.fam')
        def remove_suffix(x: str) -> str:
            return x.split('_')[0]
        df['IID'] = df['IID'].apply(remove_suffix)
        write_fam(df=df, path=f'{self.bfile}.fam')

    def build_keep_samples_phen(self):
        id_link_df = pd.read_excel(self.id_link_xslx)
        phenotype_df = pd.read_excel(self.phenotype_xslx, usecols=[self.uuid_column, self.phenotype_column])
        fam_df = read_fam(path=f'{self.bfile}.fam')

        phenotype_df = phenotype_df.merge(
            right=id_link_df,
            how='inner',
            left_on=self.uuid_column,
            right_on=self.uuid_column,
        )
        phenotype_df.drop(columns=[self.uuid_column], inplace=True)
        
        fam_df = fam_df.merge(
            right=phenotype_df,
            how='inner',
            left_on='IID',
            right_on=self.tpmi_id_column,
        )
        fam_df['PHENO'] = fam_df[self.phenotype_column]
        fam_df.drop(columns=[self.tpmi_id_column, self.phenotype_column], inplace=True)

        self.keep_samples_phen = f'{self.workdir}\\keep_samples.phen'
        fam_df[['FID', 'IID', 'PHENO']].to_csv(self.keep_samples_phen, index=False, header=False, sep=' ')

    def plink_keep(self):
        out = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='_keep',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.bfile}',
            f'--keep {self.keep_samples_phen}',
            '--make-bed',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.bfile = out

    def plink_pheno(self):
        out = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='_pheno',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.bfile}',
            f'--pheno {self.keep_samples_phen}',
            f'--mpheno 1',
            '--make-bed',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.bfile = out

    def plink_qc(self):
        out = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='_qc',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.bfile}',
            f'--maf {self.minimum_minor_allele_frequency}',
            f'--geno {self.maximum_per_variant_missing_genotype_rate}',
            f'--mind {self.maximum_per_sample_missing_genotype_rate}',
            f'--hwe {self.hardy_weinberg_p_value_threshold}',
            '--make-bed',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.bfile = out
    
    def plink_assoc(self):
        out = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='_assoc',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.bfile}',
            '--assoc',
            '--adjust',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        for src_ext, dst_ext in [
            ('.assoc', '.txt'),
            ('.assoc.adjusted', '-adjusted.txt'),
            ('.log', '.log'),
        ]:
            src = f'{out}{src_ext}'
            dst = f'{self.outdir}\\association{dst_ext}'
            shutil.copy(src, dst)

    def sort_by_association_p_values(self):
        df = pd.read_csv(join(self.outdir, 'association.txt'), sep=r'\s+')
        df.dropna(subset=['P'], inplace=True)
        df.sort_values(by='P', ascending=True, inplace=True)
        df.to_csv(join(self.outdir, 'association-sorted.csv'), index=False)


def read_fam(path: str) -> pd.DataFrame:
    # '\t' or ' ' could be used as separator, determine it by the first line
    with open(path) as fh:
        first_line = fh.readline()
        if '\t' in first_line:
            sep = '\t'
        else:
            sep = ' '
    df = pd.read_csv(path, sep=sep, header=None)
    df.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
    return df


def write_fam(df: pd.DataFrame, path: str):
    df.to_csv(path, index=False, header=False, sep=' ')


def read_bim(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']
    return df
