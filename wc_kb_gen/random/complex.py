"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
::Date: 2018-08-02
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_kb_gen

class ComplexGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates ComplexSpeciesType objects for ribosomes (and any other complexes that user specifies in assigned_complexes). 

    Options:
        * assigned_complexes (:obj:'list'): A list of the names of complexes to be created
    """


    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options
        assigned_complexes = options.get('assigned_complexes', ['complex_70S_IA', 'complex_70S_A'])

        options['assigned_complexes'] = assigned_complexes

    def gen_components(self):
        """ Creates complex objects based on assigned complexes list"""
        cell = self.knowledge_base.cell
        assigned_complexes = self.options['assigned_complexes']
        for comp in assigned_complexes:
            comp_species = cell.species_types.get_or_create(id = comp, __type=wc_kb.ComplexSpeciesType)
            comp_species.concentration = 1e-2
            #comp_species.formation_process = wc_kb.ComplexFormationType.process_RibosomeAssembly
