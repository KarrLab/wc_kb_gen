"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
::Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import wc_kb_gen
import numpy


class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates observable objects for proteins and tRNAs that are assigned to specific functions. Adds these observables to the knowledge base.

    Options:
        * assigned_trnas (:obj:'list'): A list of the names of trnas to be created
        * assigned_proteins (:obj: 'list'): A list of the names of proteins to be created
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        options = self.options
        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]
        default_trnas = []
        for codon in codons:
            default_trnas.append('tRNA_'+codon)

        assigned_trnas = options.get('assigned_trnas', default_trnas)

        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.RnaSpeciesType)

        count = 0
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                count += 1

        assert (len(assigned_trnas) <= count)
        options['assigned_trnas'] = assigned_trnas

        assigned_proteins = options.get('assigned_proteins', ['IF', 'EF', 'RF',
                                                              'deg_ATPase', 'deg_protease', 'deg_rnase',
                                                              'rna_poly'])

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)

        assert(len(assigned_proteins) <= len(prots))
        options['assigned_proteins'] = assigned_proteins

    def gen_components(self):
        """ Takes random samples of the generated rnas and proteins and assigns them functions based on the included list of proteins and rnas"""
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')

        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.RnaSpeciesType)

        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                trnas.append(rna)

        sampled_trnas = numpy.random.choice(
            trnas, len(assigned_trnas), replace=False)

        assigned_trnas = iter(assigned_trnas)

        for rna in sampled_trnas:
            rna_name = next(assigned_trnas)
            observable = cell.observables.get_or_create(id=rna_name+'_obs')
            observable.name = rna_name
            observable.species.append(
                wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=rna, compartment=cytosol), coefficient=1))

        sampled_proteins = numpy.random.choice(
            prots, len(assigned_proteins), replace=False)

        assigned_proteins = iter(assigned_proteins)
        for protein in sampled_proteins:
            protein_name = next(assigned_proteins)
            observable = cell.observables.get_or_create(id=protein_name+'_obs')
            observable.name = protein_name
            observable.species.append(
                wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
