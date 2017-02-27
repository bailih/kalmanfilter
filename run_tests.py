'''
    Small script to run all unit tests.

    We could use loadTestFromModule to make this script shorter, but then there is no visual
    distinction between tests from different classes.
'''

import unit_tests
import unittest
import re

print('Running tests!\n\n')
for name, cls in unit_tests.__dict__.items():
    if isinstance(cls, type):
        # Separate and put to upper-case
        title = re.sub('(.)([A-Z][a-z]+)', r'\1 \2', name)
        title = re.sub('([a-z0-9])([A-Z])', r'\1 \2', title).upper()
        # Pad title with '-', fitting to 70 characters total
        print('{:-^70}\n'.format(title))
        tests = unittest.TestLoader().loadTestsFromTestCase(cls)
        unittest.TextTestRunner(verbosity=2).run(tests)
        print()
