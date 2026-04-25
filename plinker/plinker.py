import os
import shutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from copy import copy
from typing import List
from os.path import abspath, join
from qmplot import manhattanplot, qqplot
from .utils import edit_fpath
from .template import Processor


class Plinker(Processor):

    plink_path: str
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
    covariate_columns: List[str]
    num_pc_covariates: int

    def main(
            self,
            plink_path: str,
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
            pi_hat: float,
            covariate_columns: List[str],
            num_pc_covariates: int):

        self.plink_path = abspath(plink_path)
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
        self.covariate_columns = covariate_columns
        self.num_pc_covariates = num_pc_covariates

        for phenotype_column in self.phenotype_columns:
            settings = copy(self.settings)
            settings.outdir = join(self.outdir, phenotype_column)
            settings.workdir = join(self.workdir, phenotype_column)
            for d in [settings.workdir, settings.outdir]:
                os.makedirs(d, exist_ok=True)
            OnePhenotypePipeline(settings).main(
                plink_path=self.plink_path,
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
                pi_hat=self.pi_hat,
                covariate_columns=self.covariate_columns,
                num_pc_covariates=self.num_pc_covariates)


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
    covariate_columns: List[str]
    num_pc_covariates: int

    sample_df: pd.DataFrame
    keep_samples_phen_file: str
    analysis_ready_bfile: str
    ld_pruned_variants_txt: str
    covariates_txt: str

    def main(
            self,
            plink_path: str,
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
            pi_hat: float,
            covariate_columns: List[str],
            num_pc_covariates: int):

        self.plink_path = plink_path
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
        self.covariate_columns = copy(covariate_columns)  # copy to avoid leaking the list to other iterations
        self.num_pc_covariates = num_pc_covariates
        
        self.copy_and_clean_bfile()
        phenotype_is = self.get_phenotype_type()
        if phenotype_is == 'continuous':
            self.logger.info(f'Continuous phenotype "{self.phenotype_column}" is not supported yet. Skipping...')
            return
        elif phenotype_is == 'categorical':
            self.logger.info(f'Phenotype column "{self.phenotype_column}" contains string values (categorical), which is not supported. Skipping...')
            return

        self.build_sample_df()
        self.plink_keep()
        self.plink_pheno()
        self.plink_qc()

        self.plink_linkage_disequilibrium_pruning()  # reduce the number of variants, for kinship and PCA

        self.plink_kinship()
        
        self.plink_assoc()

        self.plink_pca()
        self.build_covariate_file()
        self.plink_logistic()

        self.sample_df.to_csv(join(self.outdir, 'sample_table.csv'), index=False)

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

    def get_phenotype_type(self) -> str:
        phenotypes = pd.read_excel(self.phenotype_xslx)[self.phenotype_column]
        phenotypes = phenotypes.dropna().unique()
        for value in phenotypes:
            try:
                float(value)  # anything that cannot be converted to float is categorical
            except ValueError:
                return 'categorical'
        
        phenotypes = set(map(float, phenotypes))  # now all values can be converted to float
        if phenotypes == {0.0, 1.0} or phenotypes == {1.0, 2.0}:
            return 'binary'
        else:
            return 'continuous'

    def build_sample_df(self):
        """
        self.sample_df defines the sample space for analysis
        It is the intersection of three data inputs: id_link_xslx, phenotype_xslx, and bfile.fam
        """
        self.logger.info(f'Building the sample dataframe and keep_samples.phen file')

        usecols = [self.uuid_column, self.phenotype_column] + self.covariate_columns
        phenotype_df = pd.read_excel(self.phenotype_xslx, usecols=usecols)
        id_link_df = pd.read_excel(self.id_link_xslx)
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
        
        fam_df = fam_df.merge(
            right=phenotype_df,
            how='inner',
            left_on='IID',
            right_on=self.tpmi_id_column,
        )

        fam_df['PHENO'] = fam_df[self.phenotype_column]  # update with the user-specified phenotype column

        self.sample_df = fam_df
        self.keep_samples_phen_file = join(self.workdir, 'keep_samples.phen')
        self.sample_df[['FID', 'IID', 'PHENO']].to_csv(self.keep_samples_phen_file, index=False, header=False, sep=' ')        

    def plink_keep(self):
        out = edit_fpath(
            fpath=self.bfile,
            old_suffix='',
            new_suffix='_keep',
            dstdir=self.workdir
        )
        lines = [
            f'{self.plink_path}',
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
            f'{self.plink_path}',
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
            f'{self.plink_path}',
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
            f'{self.plink_path}',
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
            f'{self.plink_path}',
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

        df = pd.read_csv(f'{kinship_out}.genome', sep=r'\s+')
        df.to_csv(join(self.outdir, 'kinship.csv'), index=False)

    def plink_assoc(self):
        out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_chi_square',
            dstdir=self.workdir
        )
        lines = [
            f'{self.plink_path}',
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
            association_file=f'{out}.assoc',
            adjusted_association_file=f'{out}.assoc.adjusted',
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
            f'{self.plink_path}',
            f'--bfile {self.analysis_ready_bfile}',
            f'--extract {self.ld_pruned_variants_txt}',
            '--pca header',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        pca_df = pd.read_csv(f'{out}.eigenvec', sep=r'\s+')
        usecols = ['IID']
        for i in range(1, self.num_pc_covariates + 1):
            if f'PC{i}' in pca_df.columns:
                usecols.append(f'PC{i}')
                self.covariate_columns.append(f'PC{i}')  # add PCs to the covariate columns
        pca_df = pca_df[usecols]

        self.sample_df = self.sample_df.merge(
            right=pca_df,
            how='left',
            left_on='IID',
            right_on='IID',
        )

        output_prefix = join(self.outdir, 'pca')
        plt.figure(figsize=(6/2.54, 6/2.54), dpi=600)
        sns.scatterplot(data=self.sample_df, x='PC1', y='PC2', hue=self.phenotype_column, marker='o', alpha=0.8, palette=['#478EC9', '#EA4242'])
        plt.savefig(f'{output_prefix}.png', dpi=600, bbox_inches='tight')
        plt.legend().remove()
        plt.savefig(f'{output_prefix}_clean.png', dpi=600, bbox_inches='tight')
        plt.close()
    
    def build_covariate_file(self):
        # if no sex covariate is provided, add it from the fam file
        has_sex_covar = False
        for c in self.covariate_columns:
            if c.lower() == 'sex':
                has_sex_covar = True
                break
        if not has_sex_covar:
            fam_file_sex = self.sample_df['SEX_']  # encoded as {1,2} in the fam file
            self.sample_df['Sex'] = fam_file_sex - 1  # convert to {0,1}
            self.covariate_columns.append('Sex')
        
        self.covariates_txt = join(self.workdir, 'covariates.txt')
        self.sample_df.to_csv(self.covariates_txt, index=False, sep='\t')
    
    def plink_logistic(self):
        covar_names = ','.join(self.covariate_columns)
        out = edit_fpath(
            fpath=self.analysis_ready_bfile,
            old_suffix='',
            new_suffix='_logistic',
            dstdir=self.workdir
        )
        lines = [
            f'{self.plink_path}',
            f'--bfile {self.analysis_ready_bfile}',
            f'--covar {self.covariates_txt}',
            f'--covar-name {covar_names}',
            '--logistic hide-covar',
            f'--out {out}',
            f'1> {out}.stdout',
            f'2> {out}.stderr',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
        self.qmplot(
            association_file=f'{out}.assoc.logistic',
            output_prefix=join(self.outdir, 'logistic')
        )
        self.sort_and_filter_association_results(
            association_file=f'{out}.assoc.logistic',
            adjusted_association_file=None,
            output_prefix=join(self.outdir, 'logistic')
        )

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
            data=df['P'].dropna(),
            ax=ax,
            marker='o',
            xlabel=r"Expected $-log_{10}{(P)}$",
            ylabel=r"Observed $-log_{10}{(P)}$",
        )
        plt.savefig(f'{output_prefix}_qq.png', dpi=600, bbox_inches='tight')
        plt.close()

    def sort_and_filter_association_results(
            self,
            association_file: str,
            adjusted_association_file: Optional[str],
            output_prefix: str):
        
        df = pd.read_csv(association_file, sep=r'\s+')
        df = df[df['P'] <= self.association_p_value_threshold]
        df.sort_values(by='P', ascending=True, inplace=True)
        df.to_csv(f'{output_prefix}.csv', index=False)

        if adjusted_association_file is None:
            return
        
        snps = df['SNP'].tolist()
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
    df.columns = ['FID', 'IID', 'PID', 'MID', 'SEX_', 'PHENO']  # "SEX_" is to avoid potential collision with the SEX column from user-provided covariates
    return df


def write_fam(df: pd.DataFrame, path: str):
    df = df[['FID', 'IID', 'PID', 'MID', 'SEX_', 'PHENO']]
    df.to_csv(path, index=False, header=False, sep=' ')
