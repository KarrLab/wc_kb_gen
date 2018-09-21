from .core import RandomKbGenerator
from .genome import GenomeGenerator
from .metabolites import MetabolitesGenerator
from .properties import PropertiesGenerator
from .observables import ObservablesGenerator
from .compartments import CompartmentsGenerator
from .complex import ComplexGenerator
import logging
from os import path, remove
import sys

if path.isfile("RandomKbGenerator.log"):
    remove("RandomKbGenerator.log")
root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.DEBUG)

root_logger_file_handler = logging.FileHandler("RandomKbGenerator.log")
root_logger_file_handler.setLevel(logging.DEBUG)

root_logger_output_handler = logging.StreamHandler(sys.stdout)
root_logger_output_handler.setLevel(logging.DEBUG)

formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
root_logger_file_handler.setFormatter(formatter)
root_logger_output_handler.setFormatter(formatter)

root_logger.addHandler(root_logger_file_handler)
root_logger.addHandler(root_logger_output_handler)
