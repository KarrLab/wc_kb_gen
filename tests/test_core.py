""" Tests of the knowledge base generator

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-05-04
:Copyright: 2018, Karr Lab
:License: MIT
"""

import unittest
import wc_kb
import wc_kb_gen


class TestKbGenerator(unittest.TestCase):
    def test_KbGenerator(self):
        gen = wc_kb_gen.KbGenerator()

        self.assertEqual(gen.component_generators, [])
        self.assertEqual(gen.options, {'id': None, 'name': None, 'version': None})

    def test_ModelGenerator_run(self):
        gen = wc_kb_gen.KbGenerator(options={
            'id': 'test_kb',
        })
        kb = gen.run()

        self.assertEqual(kb.id, 'test_kb')
        self.assertEqual(kb.version, None)

        gen = wc_kb_gen.KbGenerator(options={
            'id': 'test_kb',
            'version': '0.0.1',
        })
        kb = gen.run()

        self.assertEqual(kb.id, 'test_kb')
        self.assertEqual(kb.version, '0.0.1')

    def test_ModelGenerator_run_with_components(self):
        class TestComponentGenerator1(wc_kb_gen.KbComponentGenerator):
            def get_data(self):
                self.data = self.options['compartment_id']

            def process_data(self):
                self.data *= 2

            def gen_components(self):
                self.knowledge_base.cell.compartments.create(id=self.data, name='Cytosol')

        class TestComponentGenerator2(wc_kb_gen.KbComponentGenerator):
            def get_data(self):
                self.data = self.options['compartment_id']

            def process_data(self):
                self.data *= 3

            def gen_components(self):
                self.knowledge_base.cell.compartments.create(id=self.data, name='Membrane')

        component_generators = [
            TestComponentGenerator1,
            TestComponentGenerator2,
        ]

        options = {
            'id': 'test_kb',
            'version': '0.0.1',
            'component': {
                'TestComponentGenerator1': {
                    'compartment_id': 'c'
                },
                'TestComponentGenerator2': {
                    'compartment_id': 'm'
                },
            },
        }
        kb = wc_kb_gen.KbGenerator(component_generators=component_generators, options=options).run()

        self.assertEqual(kb.cell.compartments.get_one(id='cc').name, 'Cytosol')
        self.assertEqual(kb.cell.compartments.get_one(id='mmm').name, 'Membrane')

    def test_KbComponentGenerator(self):
        kb = wc_kb.KnowledgeBase()
        kb.cell = wc_kb.Cell()

        class TestKbComponentGenerator(wc_kb_gen.KbComponentGenerator):
            def get_data(self):
                self.data = self.options['value']

            def process_data(self):
                self.data *= 2

            def gen_components(self):
                self.knowledge_base.cell.compartments.create(id=self.data)

        gen = TestKbComponentGenerator(kb, options={'value': 'c'})
        gen.run()

        self.assertIsInstance(kb.cell.compartments.get_one(id='cc'), wc_kb.Compartment)
