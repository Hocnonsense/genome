# -*- coding: utf-8 -*-
"""
* @Date: 2023-10-22 20:59:23
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-14 10:43:11
* @FilePath: /genome/tests/__init__.py
* @Description:
"""
# """

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable


sys.path.insert(0, str(Path(__file__).parent.parent))

test_temp = Path(__file__).parent / "temp"
test_files = Path(__file__).parent / "file"


def temp_output(f: Callable):
    def _f(*args, **kwargs):
        with TemporaryDirectory(prefix=str(test_temp)) as _td:
            return f(*args, **kwargs, test_temp=Path(_td))

    _f._func = f  # type: ignore[attr-defined]
    return _f
