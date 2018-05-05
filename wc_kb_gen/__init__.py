import pkg_resources

with open(pkg_resources.resource_filename('wc_kb_gen', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import KbGenerator, KbComponentGenerator
