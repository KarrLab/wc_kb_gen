from wc_model_gen import prokaryote
import unittest
import wc_lang
import wc_kb_gen
import wc_kb
import unittest
import os
import shutil
import tempfile

class KbIOTestCase(unittest.TestCase):

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
                             'num_genes': 30,
                             'num_tRNA': 4,
                             'translation_table': 4,
                             'seq_path': os.path.join(self.dir, 'kb_reduced_seq.fna')},
                         'ObservablesGenerator': {
                             'genetic_code': 'reduced'},
                         }}).run()

        kb = wc_kb_gen.random.RandomKbGenerator(options={
                     'component': {
                         'GenomeGenerator': {
                             'num_genes': 100,
                             'mean_gene_len': 70,
                             'seq_path': os.path.join(self.dir, 'kb_seq.fna'),
                             }}}).run()

        self.assertIsInstance(kb_reduced, wc_kb.core.KnowledgeBase)
        self.assertIsInstance(kb, wc_kb.core.KnowledgeBase)

        """ Write KBs to disk """
        wc_kb.io.Writer().run(knowledge_base=kb_reduced,
                              core_path=os.path.join(self.dir, 'kb_reduced.xlsx'),
                              set_repo_metadata_from_path=False)

        wc_kb.io.Writer().run(knowledge_base=kb,
                           core_path=os.path.join(self.dir, 'kb.xlsx'),
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
