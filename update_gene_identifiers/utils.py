# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 06 2021
@created: Jan 06 2021

Useful functions.

Functions
    discard_rows_with_less_info
"""

def discard_rows_with_less_info(df, fields_subset, fields_info, verbose=True):
    df_new = df.copy()
    df_new["info"] = sum([(~df_new[field_info].isnull()).astype(int) for field_info in fields_info])
    df_new = df_new.sort_values("info")
    df_new = df_new.drop_duplicates(subset=fields_subset, keep="last")
    if verbose:
        print("discarded %d/%d rows with lesser info" % (df.shape[0]-df_new.shape[0], df.shape[0]))
    del df_new["info"]
    return df_new
