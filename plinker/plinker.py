import os
import shutil
import pandas as pd
from copy import copy
from typing import List
from os.path import basename, dirname, join
from .utils import edit_fpath
from .template import Processor


PLINK_EXE = join(dirname(dirname(__file__)), 'plink.exe')


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
    association_p_value_threshold: float

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
            hardy_weinberg_p_value_threshold: float,
            association_p_value_threshold: float):

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
        self.association_p_value_threshold = association_p_value_threshold
        
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
                hardy_weinberg_p_value_threshold=self.hardy_weinberg_p_value_threshold,
                association_p_value_threshold=self.association_p_value_threshold)


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

    phenotype_type: str
    keep_samples_phen_file: str
    association_file_prefix: str 

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
            hardy_weinberg_p_value_threshold: float,
            association_p_value_threshold: float):

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
        self.association_p_value_threshold = association_p_value_threshold
        
        self.copy_and_clean_bfile()
        self.determine_phenotype_type()
        if self.phenotype_type == 'continuous':
            self.logger.info(f'Continuous phenotype "{self.phenotype_column}" is not supported yet, skipping')
            return
        self.build_keep_samples_phen_file()
        self.plink_keep()
        self.plink_pheno()
        self.plink_qc()
        self.plink_assoc()
        self.sort_and_filter_results()

    def copy_and_clean_bfile(self):
        for ext in ['.bed', '.bim', '.fam']:
            shutil.copy(f'{self.bfile}{ext}', self.workdir)
        self.bfile = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='',
            dstdir=self.workdir
        )

        self.logger.info(f'Cleaning up suffix in the IID column of {self.bfile}.fam')  # lousy problem caused by the TPMI team
        df = read_fam(path=f'{self.bfile}.fam')
        def remove_suffix(x: str) -> str:
            return x.split('_')[0]
        df['IID'] = df['IID'].apply(remove_suffix)
        write_fam(df=df, path=f'{self.bfile}.fam')

    def determine_phenotype_type(self):
        phenotypes = pd.read_excel(self.phenotype_xslx)[self.phenotype_column]
        unqiue_phenotypes = set(phenotypes.dropna().astype(float))
        if unqiue_phenotypes == {0, 1} or unqiue_phenotypes == {1, 2}:
            self.phenotype_type = 'binary'
        else:
            self.phenotype_type = 'continuous'

    def build_keep_samples_phen_file(self):
        self.logger.info(f'Building keep_samples.phen file')
        id_link_df = pd.read_excel(self.id_link_xslx)
        phenotype_df = pd.read_excel(self.phenotype_xslx, usecols=[self.uuid_column, self.phenotype_column])
        fam_df = read_fam(path=f'{self.bfile}.fam')

        phenotype_df = phenotype_df.dropna(subset=[self.phenotype_column])

        if set(phenotype_df[self.phenotype_column]) == {0, 1}:
            phenotype_df[self.phenotype_column] = phenotype_df[self.phenotype_column] + 1  # plink requires the phenotype to be 1=control, 2=case

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

        fam_df['PHENO'] = fam_df[self.phenotype_column]  # update with the user-specified phenotype column

        fam_df.drop(columns=[self.tpmi_id_column, self.phenotype_column], inplace=True)

        self.keep_samples_phen_file = join(self.workdir, 'keep_samples.phen')
        fam_df[['FID', 'IID', 'PHENO']].to_csv(self.keep_samples_phen_file, index=False, header=False, sep=' ')

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
            f'--keep {self.keep_samples_phen_file}',
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
            f'--pheno {self.keep_samples_phen_file}',
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
        self.association_file_prefix = out

    def sort_and_filter_results(self):
        df = pd.read_csv(f'{self.association_file_prefix}.assoc', sep=r'\s+')
        df = df[df['P'] <= self.association_p_value_threshold]
        df.sort_values(by='P', ascending=True, inplace=True)
        df.to_csv(join(self.outdir, 'association.csv'), index=False)

        snps = df['SNP'].tolist()

        df = pd.read_csv(f'{self.association_file_prefix}.assoc.adjusted', sep=r'\s+')
        df = df[df['SNP'].isin(snps)]
        df.to_csv(join(self.outdir, 'association-adjusted.csv'), index=False)

        shutil.copy(f'{self.association_file_prefix}.log', join(self.outdir, 'association.log'))


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
    assert len(df.columns) == 6, f'The number of columns in the dataframe must be 6, but got {len(df.columns)}: {list(df.columns)}'
    df.to_csv(path, index=False, header=False, sep=' ')


def read_bim(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']
    return df
