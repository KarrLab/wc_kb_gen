"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
::Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import wc_kb_gen
import numpy
from wc_onto import onto as wcOntology
from wc_utils.util.ontology import are_terms_equivalent

class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates observable objects for proteins and tRNAs and complexes that are assigned to specific functions. Adds these observables to the knowledge base.

    Options:
        * assigned_trnas (:obj:`list`): A list of the names of trnas to be created
        * assigned_proteins (:obj:`list`): A list of the names of proteins to be created
        * assigned_complexes (:obj:`list`): A list of the names of complexes to be created
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        options = self.options

        genetic_code = options.get('genetic_code', 'normal')
        assert(genetic_code in ['normal', 'reduced'])
        options['genetic_code'] = genetic_code

        genetic_code = options['genetic_code']
        if genetic_code=='normal':
            bases = "TCAG"
            codons = [a + b + c for a in bases for b in bases for c in bases]
            default_trnas = []
            for codon in codons:
                default_trnas.append('tRNA_'+codon)

        elif genetic_code=='reduced':
            default_trnas = ['tRNA_ATC','tRNA_CTG','tRNA_ATG','tRNA_ACG']

        assigned_trnas = options.get('assigned_trnas', default_trnas)
        options['assigned_trnas'] = assigned_trnas

        assigned_proteins = options.get('assigned_proteins', [
            'translation_init_factors',
            'translation_elongation_factors',
            'translation_release_factors',
            'degrade_rnase',
            'degrade_protease',
            'rna_polymerase',
            'aminoacyl_synthetase'])

        prots = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote.ProteinSpeciesType)
        assert(len(assigned_proteins) <= len(prots))
        options['assigned_proteins'] = assigned_proteins

        assigned_complexes = options.get('assigned_complexes', ['ribosome'])
        options['assigned_complexes'] = assigned_complexes

    def gen_components(self):
        """ Takes random samples of the generated rnas and proteins and assigns them functions based on the included list of proteins and rnas"""

        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        genetic_code = self.options['genetic_code']
        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']
        assigned_complexes = self.options['assigned_complexes']
        prots = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)

        trnas = []
        for rna in rnas:
            if are_terms_equivalent(rna.type, wcOntology['WC:tRNA']):
                trnas.append(rna)

        if genetic_code=='normal':
            codons = {
                'I': ['ATT', 'ATC', 'ATA'],
                'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
                'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                'F': ['TTT', 'TTC'],
                'M': ['ATG'],
                'C': ['TGT', 'TGC'],
                'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                'Y': ['TAT', 'TAC'],
                'W': ['TGG'],
                'Q': ['CAA', 'CAG'],
                'N': ['AAT', 'AAC'],
                'H': ['CAT', 'CAC'],
                'E': ['GAA', 'GAG'],
                'D': ['GAT', 'GAC'],
                'K': ['AAA', 'AAG'],
                'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']}

        elif genetic_code=='reduced':
            codons = {
                'I': ['ATC'],
                'L': ['CTG'],
                'M': ['ATG'],
                'T': ['ACG']}

        for aa in codons:
            rna = numpy.random.choice(trnas)
            trnas.remove(rna)
            species = rna.species.get_or_create(compartment=cytosol)
            expression = wc_kb.core.ObservableExpression(
                    expression=species.id(), species=[species])

            for i in range(len(codons[aa])):
                codon = codons[aa][i]
                rna_name = 'tRNA_'+codon
                observable = cell.observables.get_or_create(id=rna_name+'_obs')
                observable.name = rna_name
                observable.expression = expression

        sampled_proteins = numpy.random.choice(
            prots, len(assigned_proteins), replace=False)
        assigned_proteins = iter(assigned_proteins)
        for protein in sampled_proteins:
            protein_name = next(assigned_proteins)
            observable = cell.observables.get_or_create(id=protein_name+'_obs')
            observable.name = protein_name
            species = protein.species.get_or_create(compartment=cytosol)
            observable.expression = wc_kb.core.ObservableExpression(
                    expression=species.id(), species=[species])

        for comp in assigned_complexes:
            comp_species = cell.species_types.get_or_create(
                id=comp, __type=wc_kb.core.ComplexSpeciesType)
            observable = cell.observables.get_or_create(id=comp+'_obs')
            observable.name = comp
            species = comp_species.species.get_one(compartment=cytosol)
            observable.expression = wc_kb.core.ObservableExpression(
                    expression=species.id(), species=[species])
