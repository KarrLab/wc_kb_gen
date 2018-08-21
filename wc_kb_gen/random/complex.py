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
            'assigned_complexes', ['complex_70S'])

        options['assigned_complexes'] = assigned_complexes

    def gen_components(self):
        """ Creates complex objects based on assigned complexes list"""
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        assigned_complexes = self.options['assigned_complexes']
        for complex_name in assigned_complexes:
            complex = cell.species_types.get_or_create(
                id=complex_name, __type=wc_kb.ComplexSpeciesType)
            complex.formation_process = 7
            complex.concentration = 1e-2
            prot = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)[0]
            prot_species = prot.species.get_or_create(compartment=cytosol)
            prot_coeff = prot_species.species_coefficients.get_or_create(
                coefficient=1)
            complex.subunits.append(prot_coeff)

            complex_species = complex.species.get_or_create(
                compartment=cytosol)

