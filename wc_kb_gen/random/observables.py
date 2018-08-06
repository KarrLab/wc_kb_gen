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


class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates observable objects for proteins and tRNAs and complexes that are assigned to specific functions. Adds these observables to the knowledge base.

    Options:
        * assigned_trnas (:obj:'list'): A list of the names of trnas to be created
        * assigned_proteins (:obj: 'list'): A list of the names of proteins to be created
        * assigned_complexes (:obj:'list'): A list of the names of complexes to be created
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
        options['assigned_trnas'] = assigned_trnas

        assigned_proteins = options.get('assigned_proteins', ['IF', 'EF', 'RF',
                                                              'deg_ATPase', 'deg_protease', 'deg_rnase',
                                                              'rna_poly', 'aminoacyl_synthetase'])

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)

        assert(len(assigned_proteins) <= len(prots))
        options['assigned_proteins'] = assigned_proteins

        assigned_complexes = options.get('assigned_complexes', ['complex_70S'])

        options['assigned_complexes'] = assigned_complexes
                

    def gen_components(self):
        """ Takes random samples of the generated rnas and proteins and assigns them functions based on the included list of proteins and rnas"""
 
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']
        assigned_complexes = self.options['assigned_complexes']
        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.RnaSpeciesType)
        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                trnas.append(rna)
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
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],}
         
        for aa in codons:
            rna = numpy.random.choice(trnas)
            trnas.remove(rna)
            species_coefficient = wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=rna, compartment=cytosol), coefficient=1)
            for i in range(len(codons[aa])):
                codon = codons[aa][i]
                rna_name = 'tRNA_'+codon
                observable = cell.observables.get_or_create(id=rna_name+'_obs')
                observable.name = rna_name
                observable.species.append(species_coefficient)
                    
        sampled_proteins = numpy.random.choice(
            prots, len(assigned_proteins), replace=False)
        assigned_proteins = iter(assigned_proteins)
        for protein in sampled_proteins:
            protein_name = next(assigned_proteins)
            if protein_name.startswith('IF'):
                observable = cell.observables.get_or_create(id='IF_obs')
                observable.name = 'IF'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
            if protein_name.startswith('EF'):
                observable = cell.observables.get_or_create(id='EF_obs')
                observable.name = 'EF'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
            if protein_name.startswith('RF'):
                observable = cell.observables.get_or_create(id='RF_obs')
                observable.name = 'RF'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
            if protein_name.startswith('rna_poly'):
                observable = cell.observables.get_or_create(id='rna_poly_obs')
                observable.name = 'rna_poly'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
            if protein_name.startswith('deg_protease'):
                observable = cell.observables.get_or_create(id='deg_protease_obs')
                observable.name = 'deg_protease'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
            if protein_name.startswith('deg_rnase'):
                observable = cell.observables.get_or_create(id='deg_rnase_obs')
                observable.name = 'deg_rnase'
                observable.species.append(
                    wc_kb.SpeciesCoefficient(species=wc_kb.Species(species_type=protein, compartment=cytosol), coefficient=1))
                
        for comp in assigned_complexes:
            comp_species = cell.species_types.get_or_create(id = comp, __type=wc_kb.ComplexSpeciesType)
            observable = cell.observables.get_or_create(id=comp+'_obs')
            observable.name = comp
            observable.species.append(comp_species.species.get_one(compartment=cytosol).species_coefficients.get_or_create(coefficient=1))



