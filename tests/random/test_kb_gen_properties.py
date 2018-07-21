""" Tests of generation of properties

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen.random import properties
import unittest
import wc_kb


class PropertiesGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = properties.PropertiesGenerator(kb, options={
            'mean_doubling_time': 3600,
        })
        gen.run()

        self.assertEqual(cell.properties.get_one(
            id='mean_doubling_time').value, 3600)