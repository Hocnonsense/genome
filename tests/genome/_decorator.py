# -*- coding: utf-8 -*-
"""
 * @Date: 2023-10-22 20:59:23
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-04-27 18:28:22
 * @FilePath: /genome/tests/genome/_decorator.py
 * @Description:
"""
# """

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable

import pytest

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


MARK_LIMIT_RESOURCE = "-m not limit_resource" not in " ".join(sys.argv)
pytest_mark_resource = pytest.mark.skipif(
    MARK_LIMIT_RESOURCE, reason="no mark '-m \"not limit_resource\"'"
)


def temp_output(f: Callable[[Path], None]):
    def _f():
        with TemporaryDirectory(prefix=str(test_temp)) as _td:
            f(Path(_td))

    return _f
