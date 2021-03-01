import unittest
from typing import Callable, Optional

import os


# noinspection PyPep8Naming
def longRunningTest() -> Callable:
    is_skipped = (os.environ.get('CAMEL_SKIP_LONG_TESTS')) == '1'
    if is_skipped:
        return unittest.skip("Skipping long test")
    return lambda func: func


# noinspection PyPep8Naming
def resourceIntensiveTest(reason: Optional[str] = None) -> Callable:
    is_skipped = (os.environ.get('CAMEL_SKIP_RESOURCE_INTENSIVE_TESTS')) == '1'
    if is_skipped:
        return unittest.skip(f"Resource intensive test ({reason})")
    return lambda func: func
