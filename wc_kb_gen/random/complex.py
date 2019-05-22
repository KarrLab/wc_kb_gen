"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author Bilal Shaikh <bilal.shaikh@columbia.edu>
::Date: 2018-08-02
:Copyright: 2018, Karr Lab
:License: MIT
"""

from numpy import random
from scipy import stats
from wc_utils.util.units import unit_registry
import math
import numpy
import scipy.constants
import wc_kb
import wc_kb_gen


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

        if self.knowledge_base.cell.parameters.get_one(id='mean_volume') is not None:
            mean_volume = self.knowledge_base.cell.parameters.get_one(id='mean_volume').value
        else:
            mean_volume = 0.000000000000000067
            print('"mean_volume" parameter is missing, using Mycoplasma pneumoniae value (6.7E-17L).')

        for complex_name in assigned_complexes:
            cmplex_st = cell.species_types.get_or_create(id=complex_name, __type=wc_kb.core.ComplexSpeciesType)
            cmplex_specie = cmplex_st.species.get_or_create(compartment=cytosol)
            cmplex_st.formation_process = 7

            prot = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)[0]
            prot_coeff = prot.species_type_coefficients.get_or_create(coefficient=1)
            cmplex_st.subunits.append(prot_coeff)

            conc = round(abs(random.normal(loc=mean_complex_copy_number,scale=15))) / scipy.constants.Avogadro / mean_volume
            cell.concentrations.get_or_create(species=cmplex_specie, value=conc, units=unit_registry.parse_units('M'))
