import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb_gen
import wc_kb
import unittest
import tempfile
import shutil
import os

class KbAndModelIOTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_generate_read_write(self):

        """ Generate KBs """
        kb_reduced = wc_kb_gen.random.RandomKbGenerator(options={
                     'component': {
                         'GenomeGenerator': {
                             'genetic_code': 'reduced',
                             'mean_num_genes': 30,
                             'mean_gene_len': 70,
                             'num_tRNA': 4,
                             'num_rRNA': 0,
                             'num_ncRNA': 0,
                             'min_prots': 5,
                             'translation_table': 4,
                             'mean_copy_number': 100,
                             'mean_half_life': 100},
                         'PropertiesGenerator': {
                             'mean_cell_cycle_length': 100},
                         'ObservablesGenerator': {
                             'genetic_code': 'reduced'},
                         }}).run()
        kb = wc_kb_gen.random.RandomKbGenerator(options={
                     'component': {
                         'GenomeGenerator': {
                             'mean_num_genes': 30,
                             'mean_gene_len': 70,
                             'num_rRNA': 0,
                             'num_ncRNA': 0,
                             'min_prots': 5,
                             'translation_table': 4,
                             'mean_copy_number': 100,
                             'mean_half_life': 100},
                         'PropertiesGenerator': {
                             'mean_cell_cycle_length': 100},
                         }}).run()

        self.assertIsInstance(kb_reduced, wc_kb.core.KnowledgeBase)
        self.assertIsInstance(kb, wc_kb.core.KnowledgeBase)

        """ Write KBs to disk """
        wc_kb.io.Writer().run(knowledge_base=kb_reduced,
                              core_path=os.path.join(self.dir, 'kb_reduced.xlsx'),
                              seq_path=os.path.join(self.dir, 'kb_reduced_seq.fna'),
                              set_repo_metadata_from_path=False)
        wc_kb.io.Writer().run(knowledge_base=kb,
                           core_path=os.path.join(self.dir, 'kb.xlsx'),
                           seq_path=os.path.join(self.dir, 'kb_seq.fna'),
                           set_repo_metadata_from_path=False)

        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'kb_reduced.xlsx')))
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'kb_reduced_seq.fna')))
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'kb.xlsx')))
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'kb_seq.fna')))

        """ Read back KBs from disk """
        kb_reduced_read = wc_kb.io.Reader().run(core_path=os.path.join(self.dir, 'kb_reduced.xlsx'),
                                                seq_path =os.path.join(self.dir, 'kb_reduced_seq.fna'),
                                                strict=False)
        kb_read = wc_kb.io.Reader().run(core_path=os.path.join(self.dir,'kb.xlsx'),
                                        seq_path=os.path.join(self.dir,'kb_seq.fna'),
                                        strict=False)

        self.assertIsInstance(kb_reduced_read, wc_kb.core.KnowledgeBase)
        self.assertIsInstance(kb_read, wc_kb.core.KnowledgeBase)
