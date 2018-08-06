"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author Bilal Shaikh <bilal.shaikh@columbia.edu>
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
        assigned_complexes = options.get(
            'assigned_complexes', ['subunit_30S', 'subunit_50S','complex_70S'])

        options['assigned_complexes'] = assigned_complexes

    def gen_components(self):
        """ Creates complex objects based on assigned complexes list"""
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        assigned_complexes = self.options['assigned_complexes']
        for comp in assigned_complexes:
            comp_species = cell.species_types.get_or_create(
                id=comp, __type=wc_kb.ComplexSpeciesType)
            
            comp_species.concentration = 1e-2
            species = comp_species.species.get_or_create(compartment=cytosol)
            if comp.startswith('subunit'):
                comp_species.formation_process = 7  # process_RibosomeAssembly
            elif '70S' in comp:
                comp_species.formation_process = 9  # process_translation
                '''species_30S = cell.species_types.get_one(id = 'subunit_30S').species.get_one(compartment = cytosol)
                species_50S = cell.species_types.get_one(id = 'subunit_50S').species.get_one(compartment = cytosol)
                comp_species.subunits.append(species_30S.species_coefficients.get_or_create(coefficient = 1))
                comp_species.subunits.append(species_50S.species_coefficients.get_or_create(coefficient = 1))'''
