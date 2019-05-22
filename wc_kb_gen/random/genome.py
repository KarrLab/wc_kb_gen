"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""

import Bio.SeqIO
import Bio.SeqRecord
import math
import numpy
import random
import scipy.constants
import wc_kb
import wc_kb_gen
from Bio.Data import CodonTable
from Bio.Seq import Seq, Alphabet
from numpy import random
from scipy import stats
from wc_utils.util.units import unit_registry
from wc_onto import onto as wcOntology
from wc_utils.util.ontology import are_terms_equivalent


class GenomeGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Generate synthetic chromosome with randomized genes/intergenic
    regions. Creates RNA and protein objects corresponding to the genes
    this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute).

    Options:

    * num_chromosomes (:obj:`int`): number of chromosomes
    * mean_gc_frac (:obj:`float`): fraction of nucleotides which are G or C
    * num_genes (:obj:`float`): mean number of genes
    * mean_gene_len (:obj:`float`): mean codon length of a gene
    * mean_coding_frac (:obj:`float`): mean coding fraction of the genome
    * translation_table (:obj:`int`): The NCBI standard genetic code used
    * num_ncRNA (:obj:`float`): The proportion of non coding RNAs
    * num_rRNA  (:obj:`float`): The proportion of ribosomal RNAs
    * tRNA_prop (:obj:`float`): The proportion of transfer RNAs
    * five_prime_len (:obj:`int`): Average 5' UTR length for transcription units
    * three_prime_len (:obj:`int`): Average 3' UTR length for transcription units
    * operon_prop (:obj:`float`): Proportion of genes that should be in an operon (polycistronic mRNA)
    * operon_gen_num (:obj:`int`): Average number of genes in an operon
    * mean_copy_number (:obj:`float`): mean copy number of each RNA
    * mean_half_life (:obj:`float`): mean half-life of RNAs
    * genetic_code (:obj:`str`): 'normal' / 'reduced', if reduced only 'I': ['ATC'], 'L': ['CTG'],
        'M': ['ATG'], 'T': ['ACG'] codons in genome
    * seq_path (:obj:`str`): path to save genome sequence
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        # Default options are loosely  based on Escherichia coli K-12
        # Nucleic Acids Research 41:D605-12 2013

        options = self.options

        genetic_code = options.get('genetic_code', 'normal')
        assert(genetic_code in ['normal', 'reduced'])
        options['genetic_code'] = genetic_code

        num_chromosomes = options.get('num_chromosomes', 1)
        assert(num_chromosomes >= 1 and int(num_chromosomes) == num_chromosomes)
        options['num_chromosomes'] = num_chromosomes

        chromosome_topology = options.get('chromosome_topology', 'circular')
        assert(chromosome_topology in ['circular', 'linear'])
        options['chromosome_topology'] = chromosome_topology

        num_genes = options.get('num_genes', 4500)
        assert(num_genes >= 1)
        options['num_genes'] = int(num_genes)

        num_ncRNA = options.get('num_ncRNA', 10)  # not sure
        assert(isinstance(num_ncRNA, int))
        options['num_ncRNA'] = num_ncRNA

        # http://book.bionumbers.org/how-many-ribosomal-rna-gene-copies-are-in-the-genome/
        num_rRNA = options.get('num_rRNA', 7)
        assert(isinstance(num_rRNA, int))
        options['num_rRNA'] = num_rRNA

        num_tRNA = options.get('num_tRNA', 20)
        assert(isinstance(num_tRNA, int))
        options['num_tRNA'] = num_tRNA

        min_prots = options.get('min_prots', 8)
        assert(isinstance(min_prots, int))
        options['min_prots'] = min_prots

        assert((num_ncRNA + num_rRNA + num_tRNA + min_prots) <= num_genes)

        mean_gc_frac = options.get('mean_gc_frac', 0.58)
        assert(mean_gc_frac >= 0 and mean_gc_frac <= 1)
        options['mean_gc_frac'] = mean_gc_frac

        # DOI: 10.1093/molbev/msk019
        mean_gene_len = options.get('mean_gene_len', 308)  # codon length (924 bp)
        assert(mean_gene_len >= 1)
        options['mean_gene_len'] = mean_gene_len

        # DOI: 10.1007/s10142-015-0433-4
        mean_coding_frac = options.get('mean_coding_frac', 0.88)
        assert(mean_coding_frac > 0 and mean_coding_frac < 1)
        options['mean_coding_frac'] = mean_coding_frac

        translation_table = int(options.get('translation_table', 1))
        assert(translation_table in range(1, 32))
        options['translation_table'] = translation_table

        five_prime_len = int(options.get('five_prime_len', 7))
        assert(five_prime_len >= 0)
        options['five_prime_len'] = five_prime_len

        three_prime_len = int(options.get('three_prime_len', 5))  # guess
        assert(three_prime_len >= 0)
        options['three_prime_len'] = three_prime_len

        operon_prop = (options.get('operon_prop', 0.2))  # guess
        assert(operon_prop >= 0 and operon_prop <= 1)
        options['operon_prop'] = operon_prop

        operon_gen_num = int(options.get('operon_gen_num', 3))
        assert(operon_gen_num >= 2)
        options['operon_gen_num'] = operon_gen_num

        mean_rna_half_life = options.get('mean_rna_half_life', 8 * 60)
        assert(mean_rna_half_life > 0)
        options['mean_rna_half_life'] = mean_rna_half_life

        # DOI: 10.1073/pnas.0308747101
        mean_protein_half_life = options.get('mean_protein_half_life', 750 * 60)
        assert(mean_protein_half_life > 0)
        options['mean_protein_half_life'] = mean_protein_half_life

        # DOI: 10.1038/ismej.2012.94
        mean_rna_copy_number = options.get('mean_rna_copy_number', 0.4)
        assert(mean_rna_copy_number > 0)
        options['mean_rna_copy_number'] = mean_rna_copy_number

        # DOI: 10.1038/ismej.2012.94
        mean_protein_copy_number = options.get('mean_protein_copy_number', 75)
        assert(mean_protein_copy_number > 0)
        options['mean_protein_copy_number'] = mean_protein_copy_number

        seq_path = options.get('seq_path', 'rand_seq.fna')
        options['seq_path'] = seq_path

    def gen_components(self):
        self.gen_genome()
        self.gen_tus()
        self.gen_rnas_proteins()
        self.gen_concentrations()
        self.reduce_model()

    def gen_genome(self):
        '''Construct knowledge base components and generate the DNA sequence'''

        # get options
        options = self.options
        genetic_code = options.get('genetic_code')
        num_chromosomes = options.get('num_chromosomes')
        mean_gene_len = options.get('mean_gene_len')
        translation_table = options.get('translation_table')
        num_genes = options.get('num_genes')
        mean_coding_frac = options.get('mean_coding_frac')
        mean_gc_frac = options.get('mean_gc_frac')
        chromosome_topology = options.get('chromosome_topology')
        num_ncRNA = options.get('num_ncRNA')
        num_rRNA = options.get('num_rRNA')
        num_tRNA = options.get('num_tRNA')
        min_prots = options.get('min_prots')
        seq_path = options.get('seq_path')

        cell = self.knowledge_base.cell
        self.knowledge_base.translation_table = translation_table
        codon_table = CodonTable.unambiguous_dna_by_id[translation_table]

        BASES = ['A', 'C', 'G', 'T']
        PROB_BASES = [(1 - mean_gc_frac) / 2, mean_gc_frac /2, mean_gc_frac/2, (1-mean_gc_frac)/2]

        if genetic_code == 'normal':
            START_CODONS = codon_table.start_codons
            STOP_CODONS = codon_table.stop_codons

        elif genetic_code == 'reduced':
            START_CODONS = ['TTA']
            STOP_CODONS = ['TAA']

        num_genes_all = num_genes
        assignList = num_tRNA*[wcOntology['WC:tRNA']] + \
                     num_rRNA*[wcOntology['WC:rRNA']] + \
                     num_ncRNA*[wcOntology['WC:ncRNA']] + \
                     (num_genes_all-(num_ncRNA + num_tRNA + num_rRNA))*[wcOntology['WC:mRNA']]

        random.shuffle(assignList)

        # Create a chromosome n times
        dna_seqs = []
        for i_chr in range(num_chromosomes):

            num_genes = math.ceil(num_genes_all / num_chromosomes)
            gene_lens = 3 * self.rand(mean_gene_len, count=num_genes, min=2)
            intergene_lens = 3 * self.rand(mean_gene_len / mean_coding_frac * (1 - mean_coding_frac), count=num_genes)
            seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)

            # Generate seq based on random codons (NOT start/stop codons)
            seq_str = []

            if genetic_code=='normal':
                for i in range(0, seq_len, 3):
                    codon_i = random.choice(STOP_CODONS)
                    codon_i = "".join(random.choice(BASES, p=PROB_BASES, size=(3,)))
                    seq_str.append(codon_i)

            elif genetic_code=='reduced':
                for i in range(0, seq_len, 3):
                    codon_i = STOP_CODONS[0]
                    codon_i = "".join(random.choice(['ATC', 'CTG', 'ATG', 'ACG']))
                    seq_str.append(codon_i)


            seq_str = "".join(seq_str)
            seq = Seq(seq_str, Alphabet.DNAAlphabet())

            chro = cell.species_types.get_or_create(id='chr_{}'.format(i_chr + 1), __type=wc_kb.core.DnaSpeciesType)
            chro.name = 'Chromosome {}'.format(i_chr + 1)
            chro.circular = chromosome_topology == 'circular'
            chro.double_stranded = True
            chro.sequence_path = seq_path

            gene_starts = numpy.int64(numpy.cumsum(numpy.concatenate(([0], gene_lens[0:-1])) +
                                                   numpy.concatenate((numpy.round(intergene_lens[0:1] / 2), intergene_lens[1:]))))

            # creates GeneLocus objects for the genes and labels their GeneType (which type of RNA they transcribe)
            for i_gene, gene_start in enumerate(gene_starts):
                gene = self.knowledge_base.cell.loci.get_or_create(
                    id='gene_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.GeneLocus)
                gene.start = gene_start + 1  # 1-indexed
                gene.polymer = chro
                gene.end = gene.start + gene_lens[i_gene] - 1  # 1-indexed
                gene.name = 'gene {} {}'.format(i_chr+1, i_gene+1)

                if len(assignList) > 0:
                    gene.type = assignList.pop()
                else:
                    gene.type = wcOntology['WC:mRNA']

                if gene.type == wcOntology['WC:mRNA']:  # if mRNA, then set up start/stop codons in the gene

                    start_codon = random.choice(START_CODONS)
                    stop_codon = random.choice(STOP_CODONS)
                    seq_str = str(seq)

                    seq_str = seq_str[:gene.start-1] + \
                              start_codon + \
                              seq_str[gene.start+2: gene.end-3] + \
                              stop_codon + seq_str[gene.end:]

                    for i in range(gene.start+2, gene.end-3, 3):
                        while seq_str[i:i+3] in START_CODONS or seq_str[i:i+3] in STOP_CODONS:

                            if genetic_code == 'normal':
                                codon_i = "".join(random.choice(BASES, p=PROB_BASES, size=(3,)))
                            elif genetic_code == 'reduced':
                                codon_i = "".join(random.choice(['ATC', 'CTG', 'ATG', 'ACG']))

                            seq_str = seq_str[:i]+codon_i+seq_str[i+3:]

                    seq = Seq(seq_str, Alphabet.DNAAlphabet())

            dna_seqs.append(Bio.SeqRecord.SeqRecord(seq, chro.id))

        with open(seq_path, 'w') as file:
            writer = Bio.SeqIO.FastaIO.FastaWriter(
                file, wrap=70, record2title=lambda record: record.id)
            writer.write_file(dna_seqs)

    def gen_tus(self):
        """ Creates transcription units with 5'/3' UTRs, polycistronic mRNAs, and other types of RNA (tRNA, rRNA, sRNA) """

        options = self.options
        five_prime_len = options.get('five_prime_len') # 7 bp default (E. coli, wikipedia)
        three_prime_len = options.get('three_prime_len')  # 5 bp default guess
        operon_prop = options.get('operon_prop')  # 0.2 default guess
        operon_gen_num = options.get('operon_gen_num') # 3 genes default (https://academic.oup.com/gbe/article/5/11/2242/653613)

        for i_chr, chromosome in enumerate(self.knowledge_base.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType)):
            seq = chromosome.get_seq()
            i_gene = 0
            transcription_loci = []

            # Todo make this into a proper for loop that deals with repeats/additional loci
            while i_gene < len(chromosome.loci):

                gene = chromosome.loci[i_gene]

                if gene.type == wcOntology['WC:mRNA']:
                    # polycistronic mRNA (multiple GeneLocus objects per TranscriptionUnitLocus)

                    five_prime = self.rand(five_prime_len)[0]
                    three_prime = self.rand(three_prime_len)[0]
                    operon_prob = random.random()

                    # make an operon (polycistronic mRNA, put multiple genes in one TransUnitLocus)
                    if operon_prob <= operon_prop:
                        operon_genes = self.rand(operon_gen_num, min=2)[0]

                        # add 3', 5' UTRs to the ends of the transcription unit (upstream of first gene, downstream of last gene)
                        tu = self.knowledge_base.cell.loci.get_or_create(
                            id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
                        tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)

                        five_prime_start = gene.start - five_prime
                        if five_prime_start < 0:
                            five_prime_start = 0
                        tu.genes.append(gene)
                        tu.start = five_prime_start

                        for k in range(operon_genes-1):
                            i_gene += 1

                            if i_gene >= len(chromosome.loci):
                                break

                            if (chromosome.loci[i_gene]).type == wcOntology['WC:mRNA']:
                                gene = chromosome.loci[i_gene]
                                tu.genes.append(gene)
                            else:
                                break

                        three_prime_end = gene.end + three_prime
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1
                        tu.end = three_prime_end
                        transcription_loci.append(tu)

                    else:  # make an individual transcription unit for the gene
                        five_prime_start = gene.start - five_prime
                        three_prime_end = gene.end + three_prime

                        if five_prime_start < 0:
                            five_prime_start = 0
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1

                        tu = self.knowledge_base.cell.loci.get_or_create(
                            id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
                        tu.start = five_prime_start
                        tu.end = three_prime_end
                        tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)
                        tu.genes.append(gene)
                        transcription_loci.append(tu)
                        i_gene += 1

                # make a transcription unit that transcribes other types of RNA (tRNA, rRNA, sRNA)
                else:
                    tu = self.knowledge_base.cell.loci.get_or_create(
                        id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
                    tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)
                    tu.start = gene.start
                    tu.end = gene.end
                    tu.genes.append(gene)
                    transcription_loci.append(tu)
                    i_gene += 1

            for locus in transcription_loci:
                locus.polymer = chromosome

    def gen_rnas_proteins(self):
        """ Creates RNA and protein objects corresponding to genes on chromosome. """

        cell = self.knowledge_base.cell
        options = self.options
        mean_rna_half_life = options.get('mean_rna_half_life')
        mean_protein_half_life = options.get('mean_protein_half_life')

        for chromosome in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType):
            for i in range(len(chromosome.loci)):
                locus = chromosome.loci[i]

                if type(locus) == wc_kb.prokaryote_schema.TranscriptionUnitLocus:
                    tu = locus

                    # creates RnaSpeciesType for RNA sequence corresponding to gene
                    rna = self.knowledge_base.cell.species_types.get_or_create(id='rna_{}'.format(tu.id), __type=wc_kb.prokaryote_schema.RnaSpeciesType)
                    rna.name = 'rna {}'.format(tu.id)

                    # GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
                    if are_terms_equivalent(tu.genes[0].type, wcOntology['WC:mRNA']):
                        rna.type = wcOntology['WC:mRNA']
                    elif are_terms_equivalent(tu.genes[0].type, wcOntology['WC:rRNA']):
                        rna.type = wcOntology['WC:rRNA']
                    elif are_terms_equivalent(tu.genes[0].type, wcOntology['WC:tRNA']):
                        rna.type = wcOntology['WC:tRNA']
                    elif are_terms_equivalent(tu.genes[0].type, wcOntology['WC:ncRNA']):
                        rna.type = wcOntology['WC:ncRNA']

                    rna.half_life = random.normal(mean_rna_half_life, numpy.sqrt(mean_rna_half_life))
                    rna.transcription_units.append(tu)

                    if are_terms_equivalent(rna.type, wcOntology['WC:mRNA']):
                        for gene in tu.genes:
                            # creates ProteinSpecipe object for corresponding protein sequence(s)
                            prot = self.knowledge_base.cell.species_types.get_or_create(
                                id='prot_{}'.format(gene.id),
                                __type=wc_kb.prokaryote_schema.ProteinSpeciesType)
                            prot.name = 'prot_{}'.format(gene.id)

                            prot.cell = cell
                            prot.cell.knowledge_base = self.knowledge_base
                            prot.gene = gene
                            prot.rna = rna
                            prot.half_life = random.normal(mean_protein_half_life, numpy.sqrt(mean_protein_half_life))

    def gen_concentrations(self):
        """ Creates the concentration objects of RNA and protein objects """

        options = self.options
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')
        mean_rna_copy_number     = options.get('mean_rna_copy_number')
        mean_protein_copy_number = options.get('mean_protein_copy_number')

        if self.knowledge_base.cell.parameters.get_one(id='mean_volume') is not None:
            mean_volume = self.knowledge_base.cell.parameters.get_one(id='mean_volume').value
        else:
            mean_volume = 0.000000000000000067
            print('"mean_volume" parameter is missing, using Mycoplasma pneumoniae value (6.7E-17L).')

        for rna in cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType):
            rna_specie = rna.species.get_or_create(compartment=cytosol)
            conc = round(abs(random.normal(loc=mean_rna_copy_number,scale=15))) / scipy.constants.Avogadro / mean_volume
            cell.concentrations.get_or_create(id='CONC({})'.format(rna_specie.id), species=rna_specie, value=conc, units=unit_registry.parse_units('M'))

        for prot in cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            prot_specie = prot.species.get_or_create(compartment=cytosol)
            conc = round(abs(random.normal(loc=mean_protein_copy_number,scale=15))) / scipy.constants.Avogadro / mean_volume
            cell.concentrations.get_or_create(id='CONC({})'.format(prot_specie.id), species=prot_specie, value=conc, units=unit_registry.parse_units('M'))

    def reduce_model(self):
        options = self.options
        genetic_code = options.get('genetic_code')
        seq_path = options.get('seq_path')

        if genetic_code == 'normal':
            pass

        elif genetic_code == 'reduced':
            cell = self.knowledge_base.cell
            kb = self.knowledge_base

            bases = "TCAG"
            codons = [a + b + c for a in bases for b in bases for c in bases]

            dna  = kb.cell.species_types.get(__type = wc_kb.core.DnaSpeciesType)[0]
            seq_str  = str(dna.get_seq())
            seq_list = list(seq_str)

            for prot in kb.cell.species_types.get(__type = wc_kb.prokaryote_schema.ProteinSpeciesType):
                for base_num in range(prot.gene.start+2,prot.gene.end-3,3):
                    new_codon = random.choice(['ATC', 'CTG', 'ATG', 'ACG'])
                    seq_list[base_num]=new_codon[0]
                    seq_list[base_num+1]=new_codon[1]
                    seq_list[base_num+2]=new_codon[2]

            seq_str_new = ''.join(seq_list)
            seq=Seq(seq_str_new)

            dna_seqs = [Bio.SeqRecord.SeqRecord(seq, dna.id)]

            with open(seq_path, 'w') as file:
                writer = Bio.SeqIO.FastaIO.FastaWriter(
                    file, wrap=70, record2title=lambda record: record.id)
                writer.write_file(dna_seqs)

    def rand(self, mean, count=1, min=0, max=numpy.inf):
        """ Generated 1 or more random normally distributed integer(s) with standard deviation equal
        to the square root of the mean value.

        Args:
            mean (:obj:`float`): mean value
            count (:obj:`int`): number of random numbers to generate

        Returns:
            :obj:`int` or :obj:`numpy.ndarray` of :obj:`int`: random normally distributed integer(s)
        """
        a = (min-mean)/numpy.sqrt(mean)
        b = (max - mean)/numpy.sqrt(mean)
        return numpy.int64(numpy.round(stats.truncnorm.rvs(a, b, loc=mean, scale=numpy.sqrt(mean), size=count)))
