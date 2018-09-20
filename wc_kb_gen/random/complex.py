"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author Bilal Shaikh <bilal.shaikh@columbia.edu>
::Date: 2018-08-02
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_kb_gen
import random
import numpy
from numpy import random
import math
import scipy.stats as stats
import scipy.constants

class ComplexGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates ComplexSpeciesType objects for ribosomes (and any other complexes that user specifies in assigned_complexes).

    Options:
        * assigned_complexes (:obj:`list`): A list of the names of complexes to be created
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        mean_complex_copy_number = options.get('mean_complex_copy_number', 100)
        assert(mean_complex_copy_number > 0)
        options['mean_complex_copy_number'] = mean_complex_copy_number

        assigned_complexes = options.get('assigned_complexes', ['ribosome'])
        options['assigned_complexes'] = assigned_complexes

    def gen_components(self):
        """ Creates complex objects based on assigned complexes list"""
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        assigned_complexes = self.options['assigned_complexes']
        mean_complex_copy_number = self.options['mean_complex_copy_number']
        mean_volume = cell.properties.get_one(id='initial_volume').value

        for complex_name in assigned_complexes:
            cmplex_st = cell.species_types.get_or_create(id=complex_name, __type=wc_kb.core.ComplexSpeciesType)
            cmplex_specie = cmplex_st.species.get_or_create(compartment=cytosol)
            cmplex_st.formation_process = 7

            prot = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)[0]
            prot_species = prot.species.get_or_create(compartment=cytosol)
            prot_coeff = prot_species.species_coefficients.get_or_create(coefficient=1)
            cmplex_st.subunits.append(prot_coeff)

            conc = round(abs(random.normal(loc=mean_complex_copy_number,scale=15))) / scipy.constants.Avogadro / mean_volume
            cell.concentrations.get_or_create(species=cmplex_specie, value=conc, units='M')
