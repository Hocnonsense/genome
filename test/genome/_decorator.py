# -*- coding: utf-8 -*-
"""
 * @Date: 2023-10-22 20:59:23
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-22 21:03:41
 * @FilePath: /genome/test/genome/__init__.py
 * @Description:
"""
# """

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable


test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def temp_output(f: Callable[[Path], None]):
    def _f():
        with TemporaryDirectory(prefix=str(test_temp)) as _td:
            f(Path(_td))

    return _f
