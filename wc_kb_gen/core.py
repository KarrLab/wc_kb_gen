""" Base classes for generating :obj:`wc_kb`-formatted knowledge bases for whole-cell models.

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-05-04
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
from Bio.Seq import Seq, Alphabet
import random

class KbGenerator(object):
    """ Generator for knowledge bases of experimental data for whole-cell models

    Options:

    * id
    * name
    * component

    Attributes:
        component_generators (:obj:`list` of :obj:`KbComponentGenerator`): component
            generators of the knowledge base
        options (:obj:`dict`, optional): dictionary of options whose keys are the names
            of component generator classes and whose values are dictionaries of options
            for the component generator classes
    """

    DEFAULT_COMPONENT_GENERATORS = ()

    def __init__(self, component_generators=None, options=None):
        """
        Args:
            component_generators (:obj:`list` of :obj:`KbComponentGenerator`, optional): component
                generators of the knowledge base
            options (:obj:`dict`, optional): dictionary of options whose keys are the names
                of component generator classes and whose values are dictionaries of options
                for the component generator classes
        """
        if not component_generators:
            component_generators = list(self.DEFAULT_COMPONENT_GENERATORS)
        self.component_generators = component_generators

        self.options = options or {}
        self.clean_and_validate_options()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        id = options.get('id', None)
        assert(isinstance(id, str) or id is None)
        options['id'] = id

        name = options.get('name', None)
        assert(isinstance(name, str) or name is None)
        options['name'] = name

        version = options.get('version', None)
        assert(isinstance(version, str) or version is None)
        options['version'] = version

    def run(self):
        """ Generate a knowledge base of experimental data for a whole-cell model

        Returns:
            :obj:`wc_kb.core.KnowledgeBase`: knowledge base
        """
        kb = wc_kb.core.KnowledgeBase()
        kb.id = self.options.get('id')
        kb.name = self.options.get('name')
        kb.version = self.options.get('version')
        kb.cell = wc_kb.core.Cell(id='cell')

        component_options = self.options.get('component', {})
        for component_generator in self.component_generators:
            options = component_options.get(component_generator.__name__, {})
            component_generator(kb, options=options).run()

        return kb

class KbComponentGenerator(object):
    """ Base class for knowledge base component generators

    Attributes:
        knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
        options (:obj:`dict`, optional): options
    """

    def __init__(self, knowledge_base, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
            options (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.options = options or {}
        self.clean_and_validate_options()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        pass  # pragma: no cover

    def get_data(self):
        """ Get data for knowledge base components """
        pass  # pragma: no cover

    def process_data(self):
        """ Process data for knowledge base components """
        pass  # pragma: no cover

    def gen_components(self):
        """ Construct knowledge base components """
        pass  # pragma: no cover

    def run(self):
        """ Generate knowledge base components """
        self.get_data()
        self.process_data()
        self.gen_components()
