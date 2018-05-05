""" Tests of the API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-05-04
:Copyright: 2018, Karr Lab
:License: MIT
"""

import unittest
import wc_kb_gen


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(wc_kb_gen.KbGenerator, type)
