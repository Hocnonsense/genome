# -*- coding: utf-8 -*-
"""
 * @Date: 2023-10-22 20:59:23
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-07-02 22:59:42
 * @FilePath: /genome/tests/__init__.py
 * @Description:
"""
# """

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable


sys.path.insert(0, str(Path(__file__).parent.parent.parent))

test_temp = Path(__file__).parent / "temp"
test_files = Path(__file__).parent / "file"


def temp_output(f: Callable[[Path], None]):
    def _f():
        with TemporaryDirectory(prefix=str(test_temp)) as _td:
            f(Path(_td))

    return _f
