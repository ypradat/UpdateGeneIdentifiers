# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 05 2021
@created: Jan 05 2021

Test for utils module.
"""

import numpy as np
import pandas as pd
from geneids import discard_rows_with_less_info

def test_discard_rows_with_less_info():
    df = pd.DataFrame({"A": [0,0,1,1,2], "B": ["A", "A", "B", np.nan, "C"], "C": [np.nan, -1, 0, 0, 1]})
    df = _discard_rows_with_less_info(df, fields_subset=["A"], fields_info=["B", "C"])
    assert (df.values == pd.DataFrame({"A": [0,1,2], "B": ["A","B","C"], "C": [-1, 0, 1]}).values).all()
