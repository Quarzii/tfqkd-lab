"""Long scenario/reference checks run only in the explicit full group."""
import os
import unittest
full_only = unittest.skipUnless(os.environ.get('TFQKD_FULL_TESTS') == '1',
                                'full group: set TFQKD_FULL_TESTS=1')
