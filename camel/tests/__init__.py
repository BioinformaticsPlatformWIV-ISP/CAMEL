import unittest
from typing import Callable

import os


# noinspection PyPep8Naming
def longRunningTest() -> Callable:
    is_skipped = (os.environ.get('CAMEL_SKIP_LONG_TESTS')) == '1'
    if is_skipped:
        return unittest.skip("Skipping long test")
    return lambda func: func
