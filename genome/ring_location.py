# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-15 18:39:44
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-17 17:28:33
 * @FilePath: /genome/genome/ring_location.py
 * @Description:
"""

from Bio import Seq


class RingSequenceData(Seq.SequenceDataAbstractBaseClass):
    def __init__(self, data):
        self._data = bytes(data, encoding="ASCII")
        # print(self._data)
        self._length = len(self._data)
        super().__init__()

    def __len__(self):
        return self._length

    def __getitem__(self, key):
        if isinstance(key, slice):
            # ignore key.step
            # [:], [x:], [:y],
            if key.start is None or key.stop is None:
                # print("none")
                return self._data[key]
            # [x:-y],
            if key.start >= 0 and key.stop <= 0:
                # print("+->-")
                return self._data[key]
            mod_start = key.start % self._length
            mod_end = key.stop % self._length
            # [x:y]
            # >seq length=l
            # -------
            #        >seq dup
            #        -------
            # --y  x---y' y' = y+l
            if key.start >= key.stop:
                # print("start>end")
                return self[mod_start : mod_end + self._length]
            # [x:y]
            # >seq
            # -------
            #  x---y
            if mod_start < mod_end and key.stop - key.start < self._length:
                # print("0<start<end<length")
                return self._data[mod_start:mod_end]
            # [x:y]
            # >seq
            # -------
            #        >seq dup
            #        -------
            #      x---y
            # [x:y]
            # >seq          >seq
            # -------       -------
            #        >seq          >seq dup
            #        -------       -------
            #      x-----------------y
            # print("normal ring")
            return (
                self._data[mod_start:]
                + self._data * ((key.stop - key.start) // self._length - 1)
                + self._data[:mod_end]
            )
        else:
            # print("others")
            return self._data[key]

    @property
    def _2_data(self):
        return self._data * 2

    def __contains__(self, item):
        assert len(item) <= self._length
        return self._2_data.__contains__(item)

    def find(self, sub):
        """Return the lowest index in data where subsection sub is found.

        Once search cross the end, will map continue from the first of the sequence.
        Do not search the third time.

        If want search with range, please refer like data[start,end].find(sub)

        Return -1 on failure.
        """
        return self._2_data.find(sub)

    def rfind(self, sub):
        """Return the highest index in data where subsection sub is found.

        Return -1 on failure.
        """
        return self._2_data.rfind(sub)

    def split(self, sep=None, maxsplit=-1):
        """Return a list of the sections in the data, using sep as the delimiter.

        sep
          The delimiter according which to split the data.
          None (the default value) means split on ASCII whitespace characters
          (space, tab, return, newline, formfeed, vertical tab).
        maxsplit
          Maximum number of splits to do.
          -1 (the default value) means no limit.
        """
        return bytes(self).split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        """Return a list of the sections in the data, using sep as the delimiter.

        sep
          The delimiter according which to split the data.
          None (the default value) means split on ASCII whitespace characters
          (space, tab, return, newline, formfeed, vertical tab).
        maxsplit
          Maximum number of splits to do.
          -1 (the default value) means no limit.

        Splitting is done starting at the end of the data and working to the front.
        """
        return bytes(self).rsplit(sep, maxsplit)


def test_rsd_get_item():
    s1 = Seq.Seq(RingSequenceData("AAAAAAATCG"))
    s2 = Seq.Seq("AAAAAAATCG")

    assert s1[::-1] == s2[::-1]
    assert s1[1:10] == s2[1:10]
    assert s1[2:-1] == s2[2:-1]
    assert not s1[2 : -1 - len(s1)]
    assert s1[6:3] == s2[6:] + s2[:3]
    assert s1[-1 : -1 + len(s1)] == s2[-1:] + s2[: -1 + len(s2)]
    assert s1[2 : -1 + len(s1) * 0] == s1[2 : -1 + len(s1) * 1]
    assert s1[2 : -1 + len(s1) * 2] == s2[2:-1] + s1[-1 : -1 + len(s1)] * 1
    assert s1[2 : -1 + len(s1) * 3] == s2[2:-1] + s1[-1 : -1 + len(s1)] * 2
    assert s1[2 - len(s1) : 9 - len(s1)] == s2[2:9]
    assert s1[6:9] == s2[6:9]
    assert s1[6 : 9 + len(s1) * 0] == s2[6:9]
    assert s1[6 : 9 + len(s1) * 1] == s2[6:9] + s1[9 : 9 + len(s1)] * 1
    assert s1[6 : 9 + len(s1) * 2] == s2[6:9] + s1[9 : 9 + len(s1)] * 2
    assert s1[6 : 9 + len(s1) * 3] == s2[6:9] + s1[9 : 9 + len(s1)] * 3
    assert s1[6 - len(s1) : 9 - len(s1)] == s2[6:9]
