""" Tests of observables  generation"

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import unittest
from wc_kb_gen.random import observables, genome, properties, compartments, complex


class ObservablesGeneratorTestCase(unittest.TestCase):
    def setUp(self):
        kb = wc_kb.core.KnowledgeBase()
        self.cell = kb.cell = wc_kb.core.Cell()
        gen = properties.PropertiesGenerator(kb)
        gen.run()
        gen = compartments.CompartmentsGenerator(kb)
        gen.run()
        gen = genome.GenomeGenerator(
            kb, options={'num_chromosomes': 1, 'mean_num_genes': 200})
        gen.run()
        gen = complex.ComplexGenerator(kb)
        gen.run()
        gen = observables.ObservablesGenerator(
            kb, options={'assigned_proteins': ['a', 'b', 'c'], 'assigned_trnas': ['x', 'y', 'z']})
        gen.run()

    '''def test_assignment(self):
        cell = self.cell

        self.assertIsInstance(cell.species_types.get_one(
            id='a'), wc_kb.prokaryote_schema.ProteinSpeciesType)
        self.assertIsInstance(cell.species_types.get_one(
            id='x'), wc_kb.prokaryote_schema.RnaSpeciesType)'''

    def test_observables(self):
        cell = self.cell

        obs1 = cell.observables.get_one(id='a_obs')
        self.assertIsInstance(obs1, wc_kb.core.Observable)
        species_coefficient = obs1.species[0]
        self.assertIsInstance(species_coefficient, wc_kb.core.SpeciesCoefficient)
        species = species_coefficient.species
        species_type = species.species_type
        self.assertIsInstance(species, wc_kb.core.Species)
        self.assertIsInstance(species_type, wc_kb.prokaryote_schema.ProteinSpeciesType)
        '''self.assertEqual(species.id(), 'a[c]')
        self.assertEqual(species_type.id, 'a')'''
