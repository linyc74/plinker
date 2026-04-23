import os
import shutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from copy import copy
from typing import List
from os.path import dirname, join
from qmplot import manhattanplot, qqplot
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
    pi_hat: float

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
            association_p_value_threshold: float,
            pi_hat: float):

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
        self.pi_hat = pi_hat

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
                association_p_value_threshold=self.association_p_value_threshold,
                pi_hat=self.pi_hat)


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
    association_p_value_threshold: float
    pi_hat: float

    phenotype_type: str
    keep_samples_phen_df: pd.DataFrame
    keep_samples_phen_file: str
    analysis_ready_bfile: str
    ld_pruned_variants_txt: str

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
            association_p_value_threshold: float,
            pi_hat: float):

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
        self.pi_hat = pi_hat
        
        self.copy_and_clean_bfile()
        self.determine_phenotype_type()
        if self.phenotype_type == 'continuous':
            self.logger.info(f'Continuous phenotype "{self.phenotype_column}" is not supported yet, skipping')
            return

        self.build_keep_samples_phen_file()
        self.plink_keep()
        self.plink_pheno()
        self.plink_qc()

        self.plink_linkage_disequilibrium_pruning()  # reduce the number of variants, for kinship and PCA

        self.plink_kinship()
        
        self.plink_assoc()

        self.plink_pca()
        # self.plink_build_covariate_file()
        # self.plink_logistic()

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

        if not id_link_df[self.tpmi_id_column].is_unique:
            n = id_link_df[self.tpmi_id_column].duplicated().sum()
            s = 's' if n > 1 else ''
            self.logger.warning(f'The "{self.tpmi_id_column}" column in the ID link Excel file is not unique, dropping {n} duplicate{s}')
            id_link_df = id_link_df.drop_duplicates(subset=[self.tpmi_id_column])
        if not id_link_df[self.uuid_column].is_unique:
            n = id_link_df[self.uuid_column].duplicated().sum()
            s = 's' if n > 1 else ''
            self.logger.warning(f'The "{self.uuid_column}" column in the ID link Excel file is not unique, dropping {n} duplicate{s}')
            id_link_df = id_link_df.drop_duplicates(subset=[self.uuid_column])

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

        self.keep_samples_phen_df = fam_df
        self.keep_samples_phen_file = join(self.workdir, 'keep_samples.phen')
        self.keep_samples_phen_df[['FID', 'IID', 'PHENO']].to_csv(self.keep_samples_phen_file, index=False, header=False, sep=' ')

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
            '--allow-no-sex',
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
            '--allow-no-sex',
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
            '--allow-no-sex',
            '--make-bed',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.analysis_ready_bfile = out
    
    def plink_linkage_disequilibrium_pruning(self):
        out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_indep',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.analysis_ready_bfile}',
            '--indep-pairwise 200 5 0.2',  # window size 200 variants, 5 variants shift each time, r^2 > 0.2 relatedness threshold
            '--allow-no-sex',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.ld_pruned_variants_txt = f'{out}.prune.in'

    def plink_kinship(self):
        kinship_out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_kinship',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.analysis_ready_bfile}',
            f'--extract {self.ld_pruned_variants_txt}',
            '--allow-no-sex',
            '--genome',
            f'--min {self.pi_hat}',
            f'--out {kinship_out}',
            f'1> {kinship_out}.stdout',
            f'2> {kinship_out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        src = f'{kinship_out}.genome'
        dst = join(self.outdir, 'kinship.txt')
        shutil.copy(src, dst)

    def plink_assoc(self):
        out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_assoc',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.analysis_ready_bfile}',
            '--allow-no-sex',
            '--assoc',
            '--adjust',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        self.qmplot(
            association_file=f'{out}.assoc',
            output_prefix=join(self.outdir, 'chi_square')
        )
        self.sort_and_filter_association_results(
            association_prefix=out,
            output_prefix=join(self.outdir, 'chi_square')
        )

    def plink_pca(self):
        out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_pca',
            dstdir=self.workdir
        )
        lines = [
            f'{PLINK_EXE}',
            f'--bfile {self.analysis_ready_bfile}',
            f'--extract {self.ld_pruned_variants_txt}',
            '--pca header',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        df = pd.read_csv(f'{out}.eigenvec', sep=r'\s+')
        df = df.merge(
            right=self.keep_samples_phen_df,
            how='inner',
            left_on='IID',
            right_on='IID',
        ).rename(
            columns={'PHENO': self.phenotype_column}
        ).drop(
            columns=['FID_x', 'FID_y', 'PID', 'MID', 'SEX'],
        )
        output_prefix = join(self.outdir, 'pca')
        df.to_csv(f'{output_prefix}.csv', index=False)

        plt.figure(figsize=(6/2.54, 6/2.54), dpi=600)
        sns.scatterplot(data=df, x='PC1', y='PC2', hue=self.phenotype_column, marker='o', alpha=0.8, palette=['#478EC9', '#EA4242'])
        plt.savefig(f'{output_prefix}.png', dpi=600, bbox_inches='tight')
        plt.legend().remove()
        plt.savefig(f'{output_prefix}_clean.png', dpi=600, bbox_inches='tight')
        plt.close()
    
    def plink_build_covariate_file(self):
        pass
    
    def plink_logistic(self):
        pass

    def qmplot(self, association_file: str, output_prefix: str):
        df = pd.read_csv(association_file, sep=r'\s+')
        f, ax = plt.subplots(figsize=(22/2.54, 9/2.54), dpi=600)
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 7
        manhattanplot(
            data=df,
            chrom='CHR',
            pos='BP',
            pv='P',
            snp='SNP',
            color='#3B5488,#53BBD5',
            suggestiveline=1e-5,
            genomewideline=5e-8,
            sign_line_cols='#D62728,#2CA02C',  # colors used `suggestiveline` and `genomewideline`
            sign_marker_color='red'
        )
        plt.savefig(f'{output_prefix}_manhattan.png', dpi=600, bbox_inches='tight')
        plt.close()

        f, ax = plt.subplots(figsize=(6/2.54, 6/2.54), dpi=600)
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 8
        qqplot(
            data=df['P'],
            ax=ax,
            marker='o',
            xlabel=r"Expected $-log_{10}{(P)}$",
            ylabel=r"Observed $-log_{10}{(P)}$",
        )
        plt.savefig(f'{output_prefix}_qq.png', dpi=600, bbox_inches='tight')
        plt.close()

    def sort_and_filter_association_results(self, association_prefix: str, output_prefix: str):
        association_file = f'{association_prefix}.assoc'
        df = pd.read_csv(association_file, sep=r'\s+')
        df = df[df['P'] <= self.association_p_value_threshold]
        df.sort_values(by='P', ascending=True, inplace=True)
        df.to_csv(f'{output_prefix}.csv', index=False)

        snps = df['SNP'].tolist()

        adjusted_association_file = f'{association_prefix}.assoc.adjusted'
        df = pd.read_csv(adjusted_association_file, sep=r'\s+')
        df = df[df['SNP'].isin(snps)]
        df.to_csv(f'{output_prefix}_adjusted.csv', index=False)


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
