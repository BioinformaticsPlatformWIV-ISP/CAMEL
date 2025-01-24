from pathlib import Path
from typing import Union


def absolute_path_by_pathlib(relative_path: Union[str, Path]) -> Path:
    """
    Takes a relative path and returns the absolute path.
    :param relative_path: relative path
    :return: absolute path
    """
    return Path(relative_path).absolute()
