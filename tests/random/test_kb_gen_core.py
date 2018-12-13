""" Tests of knowledge base generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen import random
import obj_model
import os
import shutil
import tempfile
import unittest
import wc_kb
import wc_utils.util.string


class KbGeneratorTestCase(unittest.TestCase):
    def test(self):
        self.dir = tempfile.mkdtemp()
        gen = random.RandomKbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_cell_cycle_len': 100,
                },
                'GenomeGenerator': {
                    'num_chromosomes': 10,
                    'num_genes': 200,
                    'mean_copy_number': 200,
                    'seq_path': os.path.join(self.dir, 'kb_seq.fna'),
                },
            },
        })

        kb = gen.run()
        self.assertEqual(len(kb.cell.species_types.get(
            __type=wc_kb.core.DnaSpeciesType)), 10)

        errors = obj_model.Validator().run(kb, get_related=True)
        self.assertEqual(
            errors, None, msg=wc_utils.util.string.indent_forest(errors))

        shutil.rmtree(self.dir)

    def test_clean_and_validate_options(self):
        gen = random.RandomKbGenerator()

        gen.options = {'seed': None}
        gen.clean_and_validate_options()

        gen.options = {'seed': int(2)}
        gen.clean_and_validate_options()

        gen.options = {'seed': 'ABC'}
        with self.assertRaises(ValueError):
            gen.clean_and_validate_options()
