""" Generator for metabolites of random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import csv
import os
import pkg_resources
import wc_kb
import wc_kb_gen


class MetabolitesGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for metabolites for random in silico organisms

    Options:

    * data_path: path to CSV file with metabolite ids, name, InChI-encoded structures, and intracellular concentrations
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options
        data_path = options.get('data_path', pkg_resources.resource_filename(
            'wc_kb_gen', os.path.join('data', 'metabolites.csv')))

        assert(os.path.isfile(data_path))
        options['data_path'] = data_path

    def get_data(self):
        data_path = self.options.get('data_path')
        with open(data_path, 'r', encoding='utf-8') as file:
            self.data = []
            for met in csv.DictReader(file):
                met['Intracellular concentration (M)'] = float(
                    met['Intracellular concentration (M)'])
                self.data.append(met)

    def gen_components(self):
        """ Construct knowledge base components """
        cell = self.knowledge_base.cell

        # generate metabolites
        for met in self.data:

            stp = wc_kb.core.SpeciesTypeProperty(id='PROP(structure:{})'.format(met['Id']))
            met_species_type = cell.species_types.get_or_create(
                                __type=wc_kb.core.MetaboliteSpeciesType,
                                id=met['Id'],
                                name=met['Name'],
                                properties = [stp])

            #Intracellular specie
            met_species = wc_kb.core.Species(
                            species_type=met_species_type,
                            compartment=cell.compartments.get_or_create(
                            __type=wc_kb.core.Compartment, id='c', name='cytosol'))

            cell.concentrations.get_or_create(
                id = 'CONC({})'.format(met_species.id),
                species=met_species,
                value=met['Intracellular concentration (M)'])

            #Extracellular specie
            met_species = wc_kb.core.Species(
                            species_type=met_species_type,
                            compartment=cell.compartments.get_or_create(
                            __type=wc_kb.core.Compartment, id='e'))

            cell.concentrations.get_or_create(
                id = 'CONC({})'.format(met_species.id),
                species=met_species,
                value=met['Extracellular concentration (M)'])
