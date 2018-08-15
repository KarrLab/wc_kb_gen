"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import scipy.stats as stats
import scipy.constants
import wc_kb
import wc_kb_gen
import numpy
from numpy import random
from Bio.Seq import Seq, Alphabet
from Bio.Data import CodonTable
from types import *


class GenomeGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates synthetic chromosome with randomized genes/intergenic regions. Creates RNA and protein objects corresponding to the genes on this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute)

    Options:

    * variation_function (:obj: 'function'): a function used to calculate the standard deviation of the normal function which produces variability in the parameters
    * num_chromosomes (:obj:`int`): number of chromosomes
    * mean_gc_frac (:obj:`float`): fraction of nucleotides which are G or C
    * mean_num_genes (:obj:`float`): mean number of genes
    * mean_gene_len (:obj:`float`): mean codon length of a gene
    * mean_coding_frac (:obj:`float`): mean coding fraction of the genome
    * translation_table (:obj:'int'): The NCBI standard genetic code used
    * num_ncRNA (:obj:'float'): The proportion of non coding RNAs
    * num_rRNA  (:obj:'float'): The proportion of ribosomal RNAs
    * tRNA_prop (:obj:'float'): The proportion of transfer RNAs
    * five_prime_len (:obj:'int'): Average 5' UTR length for transcription units
    * three_prime_len (:obj:'int'): Average 3' UTR length for transcription units
    * operon_prop (:obj:'float'): Proportion of genes that should be in an operon (polycistronic mRNA)
    * operon_gen_num (:obj:'int'): Average number of genes in an operon
    * mean_copy_number (:obj:`float`): mean copy number of each RNA
    * mean_half_life (:obj:`float`): mean half-life of RNAs
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        # Default options are loosely  based on Escherichia coli K-12
        # Nucleic Acids Research 41:D605-12 2013

        options = self.options

        variation_function = options.get(
            'variation_function', lambda x: numpy.sqrt(x))
        options['variation_function'] = variation_function

        num_chromosomes = options.get('num_chromosomes', 1)
        assert(num_chromosomes >= 1 and int(
            num_chromosomes) == num_chromosomes)
        options['num_chromosomes'] = num_chromosomes

        chromosome_topology = options.get('chromosome_topology', 'circular')
        assert(chromosome_topology in ['circular', 'linear'])
        options['chromosome_topology'] = chromosome_topology

        mean_gc_frac = options.get('mean_gc_frac', 0.58)
        assert(mean_gc_frac >= 0 and mean_gc_frac <= 1)
        options['mean_gc_frac'] = mean_gc_frac

        mean_num_genes = options.get('mean_num_genes', 4500)
        assert(mean_num_genes >= 1)
        options['mean_num_genes'] = mean_num_genes

        num_ncRNA = options.get('num_ncRNA', 10)  # not sure
        options['num_ncRNA'] = num_ncRNA

        # http://book.bionumbers.org/how-many-ribosomal-rna-gene-copies-are-in-the-genome/
        num_rRNA = options.get('num_rRNA', 7)
        options['num_rRNA'] = num_rRNA

        num_tRNA = options.get('num_tRNA', 20)
        assert(num_tRNA >= 20)
        options['num_tRNA'] = num_tRNA

        # Proteins that are to be assigned to a gene
        assigned_prots = options.get('assigned_prots', [])
        options['assigned_prots'] = assigned_prots

        assert (num_ncRNA + num_rRNA + num_tRNA +
                len(assigned_prots) <= mean_num_genes)

        # DOI: 10.1093/molbev/msk019
        mean_gene_len = options.get('mean_gene_len', 300)
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

        operon_prop = (options.get('operon_prop', 0.8))  # guess
        assert(operon_prop >= 0 and operon_prop <= 1)
        options['operon_prop'] = operon_prop

        mean_operon_len = int(options.get('mean_operon_len', 4))
        assert(mean_operon_len >= 2)
        options['mean_operon_len'] = mean_operon_len

        max_operon_len = int(options.get('max_operon_len', 7))
        assert(max_operon_len >= 2)
        options['max_operon_len'] = max_operon_len

        # DOI: 10.1038/ismej.2012.94
        mean_copy_number = options.get('mean_copy_number', 0.4)
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        # DOI: 10.1073/pnas.0308747101
        mean_half_life = options.get('mean_half_life', 2.1 * 60)
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

    def gen_components(self):
        self.gen_genome2()
        # self.gen_tus()
        # self.gen_rnas_proteins()

    def gen_genome2(self):
        options = self.options
        num_chromosomes = options.get('num_chromosomes')
        num_ncRNA = options.get('num_ncRNA')
        num_rRNA = options.get('num_rRNA')
        num_tRNA = options.get('num_tRNA')
        mean_num_genes = options.get('mean_num_genes')
        assigned_prots = options.get('assigned_prots')

        min_num_genes = num_ncRNA + num_rRNA + num_tRNA + len(assigned_prots)

        num_genes = self.rand(mean=mean_num_genes,
                              min=min_num_genes, round=True)

        tu_types, operon_lens = self.gen_tu_types(num_genes)
        print(tu_types)
        print(operon_lens)

        for i in range(num_chromosomes):
            pass

    def gen_tu_types(self, num_genes: int) -> list:

        options = self.options
        operon_prop = options.get("operon_prop")
        mean_operon_len = options.get('mean_operon_len')
        max_operon_len = options.get('max_operon_len')
        min_operon_len = options.get('min_operon_len', 2)
        num_ncRNA = options.get('num_ncRNA')
        num_rRNA = options.get('num_rRNA')
        num_tRNA = options.get('num_tRNA')

        num_monocistronic_mRNA = int((
            1-operon_prop)*(num_genes - num_tRNA - num_rRNA - num_ncRNA))

        # This is the number of genes that are contained within some operon
        num_polycistronic_genes = int((
            num_genes-num_monocistronic_mRNA-num_tRNA-num_ncRNA-num_rRNA))

        operon_lens = []

        while(numpy.sum(operon_lens) <= num_polycistronic_genes-max_operon_len):
            operon_len = self.rand(
                mean=mean_operon_len, min=min_operon_len, max=max_operon_len, round=True)
            operon_lens.append(operon_len)

        operon_len = num_polycistronic_genes - numpy.sum(operon_lens)
        operon_lens.append(operon_len)

        assert (numpy.sum(operon_lens)+num_monocistronic_mRNA + num_ncRNA +
                num_rRNA + num_tRNA) == num_genes, "Total number of genes is incorrect"

        num_operons = len(operon_lens)

        tu_types = []
        for i in range(num_tRNA):
            tu_types.append("tRNA")
        for i in range(num_ncRNA):
            tu_types.append("ncRNA")
        for i in range(num_rRNA):
            tu_types.append("rRNA")
        for i in range(num_monocistronic_mRNA):
            tu_types.append("mRNA")
        for i in range(num_operons):
            tu_types.append('operon')

        numpy.random.shuffle(tu_types)

        return tu_types, operon_lens

    def gen_chromosomes(self, id_num: int, num_genes: int):
        option = self.options
        chromosome_topology = options.get('chromosome_topology')
        cell = self.kb.cell

        chromosome = cell.species_types.get_or_create(
            id='chr_{}'.format(id_num + 1), __type=wc_kb.DnaSpeciesType)

        chromosome.name = 'Chromosome {}'.format(id_num + 1)
        chromosome.circular = chromosome_topology == 'circular'
        chromosome.double_stranded = True

    def gen_transcription_units(self):
        pass

    def gen_genes(self):
        pass

    def gen_sequence(self):
        pass

    def rand(self, mean, count=1, min=0, max=numpy.inf, varfunc=None, round=False):
        """ Generated 1 or more random normally distributed values with standard deviation calculated by the provided variation function. If round flag is set, the method will return only integers

        Args:
            mean (:obj:`float`): mean value
            count (:obj:`int`): number of random numbers to generate
            min (:obj:`float`): the minimum value of the numbers to include
            max (:obj:`float`): the maximum value of the numbers to include
            varfunc (:obj: 'function'): the function to use to calculate the standard deviation
            round (:obj: 'boolean'): Whether or not to round the numbers to integers

        Returns:
            :obj:`float` or :obj:`numpy.ndarray` of :obj:`float`: random normally distributed number(s)
        """
        if not varfunc:
            varfunc = self.options.get('variation_function')
        a = (min-mean)/varfunc(mean)
        b = (max - mean)/varfunc(mean)
        if round:
            if count == 1:
                return numpy.int64(numpy.round(stats.truncnorm.rvs(a, b, size=count, loc=mean, scale=varfunc(mean))))[0]
            else:
                return numpy.int64(numpy.round(stats.truncnorm.rvs(a, b, size=count, loc=mean, scale=varfunc(mean))))
        else:
            if count == 1:
                return (stats.truncnorm.rvs(a, b, size=count, loc=mean, scale=varfunc(mean))[0])
            else:
                return (stats.truncnorm.rvs(a, b, size=count, loc=mean, scale=varfunc(mean)))
