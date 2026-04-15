import pandas as pd
from os.path import exists
from .setup import TestCase
from plinker.plinker import Plinker, read_fam, read_bim


class TestPlinker(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_small_dataset(self):
        Plinker(settings=self.settings).main(
            bfile=f'{self.indir}/wgas2',
            id_link_xslx=f'{self.indir}/wgas2_ID_link.xlsx',
            phenotype_xslx=f'{self.indir}/wgas2_phenotype.xlsx',
            uuid_column='uuid',
            tpmi_id_column='TPMINUM',
            phenotype_columns=['PHENOTYPE_A', 'PHENOTYPE_B'],
            minimum_minor_allele_frequency=0.01,
            maximum_per_variant_missing_genotype_rate=0.05,
            maximum_per_sample_missing_genotype_rate=0.01,
            hardy_weinberg_p_value_threshold=1e-6,
            association_p_value_threshold=1e-3,
        )

    def test_large_dataset(self):
        Plinker(settings=self.settings).main(
            bfile=f'{self.indir}/1kg_phase1_all',
            id_link_xslx=f'{self.indir}/1kg_phase1_all_ID_link.xlsx',
            phenotype_xslx=f'{self.indir}/1kg_phase1_all_phenotype.xlsx',
            uuid_column='uuid',
            tpmi_id_column='TPMINUM',
            phenotype_columns=['PHENOTYPE_A'],
            minimum_minor_allele_frequency=0.01,
            maximum_per_variant_missing_genotype_rate=0.05,
            maximum_per_sample_missing_genotype_rate=0.05,
            hardy_weinberg_p_value_threshold=0.001,
            association_p_value_threshold=1e-3,
        )


class TestFunctions(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_read_fam(self):
        df = read_fam(path=f'{self.indir}/wgas2.fam')
        self.assertEqual(df.shape, (90, 6))
    
    def test_read_bim(self):
        df = read_bim(path=f'{self.indir}/wgas2.bim')
        self.assertEqual(df.shape, (228694, 6))

