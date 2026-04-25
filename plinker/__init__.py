import os
import shutil
from .plinker import Plinker
from .template import Settings
from .utils import get_temp_path


def main(
        plink_path: str,
        bfile: str,
        id_link_xslx: str,
        phenotype_xslx: str,
        uuid_column: str,
        tpmi_id_column: str,
        phenotype_columns: str,
        minimum_minor_allele_frequency: float,
        maximum_per_variant_missing_genotype_rate: float,
        maximum_per_sample_missing_genotype_rate: float,
        hardy_weinberg_p_value_threshold: float,
        association_p_value_threshold: float,
        pi_hat: float,
        covariate_columns: str,
        num_pc_covariates: int,
        outdir: str,
        threads: int,
        debug: bool):

    prefix = os.path.basename(outdir)
    for c in [' ', ',', '(', ')']:
        prefix = prefix.replace(c, '_')
    workdir = get_temp_path(prefix=f'./{prefix}_')

    settings = Settings(
        workdir=workdir,
        outdir=outdir,
        threads=threads,
        debug=debug,
        mock=False)

    for d in [workdir, outdir]:
        os.makedirs(d, exist_ok=True)

    Plinker(settings).main(
        plink_path=plink_path,
        bfile=bfile,
        id_link_xslx=id_link_xslx,
        phenotype_xslx=phenotype_xslx,
        uuid_column=uuid_column,
        tpmi_id_column=tpmi_id_column,
        phenotype_columns=phenotype_columns.split(','),
        minimum_minor_allele_frequency=minimum_minor_allele_frequency,
        maximum_per_variant_missing_genotype_rate=maximum_per_variant_missing_genotype_rate,
        maximum_per_sample_missing_genotype_rate=maximum_per_sample_missing_genotype_rate,
        hardy_weinberg_p_value_threshold=hardy_weinberg_p_value_threshold,
        association_p_value_threshold=association_p_value_threshold,
        pi_hat=pi_hat,
        covariate_columns=covariate_columns.split(',') if covariate_columns.lower() != 'none' else [],
        num_pc_covariates=num_pc_covariates,
    )

    if not debug:
        shutil.rmtree(workdir)
