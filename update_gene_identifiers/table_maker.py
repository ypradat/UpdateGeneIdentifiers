# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 12 2021
@created: Jan 06 2021

Class for making correspondence tables between gene identifiers from different file formats.
"""

import os
import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

class TableMaker(object):
    def __init__(self, fields_ids: list, fields_dup: list=None, save_folder: str=None,
                    save_all_fields_ids: bool=True, save_per_field_dup: bool=False):
        """
        Make a table of all unique correspondences between identifiers. As an option, the function can record in an excel
        file the tables of duplicated values for each identifier.

        Parameters
        ----------
        df: Input dataframe
        fields_ids: List of gene identifiers
        fields_dup: List of gene identifiers for which you want to get duplicates info
        save_folder: Used only if one of save_all_fields_ids or save_per_field_dup is True
        save_all_fields_ids: If True, writes an Excel table with the correspondence table
        save_per_field_dup: If True, writes an Excel table with one sheet per field showing its duplicated values

        """
        self.fields_ids = fields_ids
        self.fields_dup = fields_dup
        self.save_folder = save_folder
        self.save_all_fields_ids = save_all_fields_ids
        self.save_per_field_dup = save_per_field_dup

    def make(self, df):
        fields_ids = self.fields_ids
        fields_dup = self.fields_dup
        save_folder = self.save_folder
        save_all_fields_ids = self.save_all_fields_ids
        save_per_field_dup = self.save_per_field_dup

        df_table = df[fields_ids].drop_duplicates()
        df_table = df_table.replace(to_replace=["-"], value=np.nan)

        if save_all_fields_ids:
            filepath = os.path.join(save_folder, "all_fields_correspondences.xlsx")
            with pd.ExcelWriter(filepath) as writer:
                pd.DataFrame({"Correspondence tables on the fields": fields_ids}).to_excel(
                    excel_writer = writer,
                    sheet_name   = 'description',
                    header       = True,
                    index        = False,
                )

                df_table.to_excel(
                    excel_writer = writer,
                    sheet_name   = "table",
                    header       = True,
                    index        = False
                )

        # duplicates info
        if fields_dup is not None:
            dt_dup = {}
            for field in fields_dup:
                df_table_field = df_table.dropna(subset=[field])
                df_table_field = df_table_field[df_table_field.duplicated(subset=[field], keep=False)]
                df_table_field = df_table_field.sort_values([field] + [x for x in fields_dup if x!=field])

                dt_dup[field] = df_table_field

            if save_per_field_dup:
                filepath = os.path.join(save_folder, "per_field_duplicates.xlsx")
                with pd.ExcelWriter(filepath) as writer:
                    pd.DataFrame({"Duplicates tables on the fields": fields_dup}).to_excel(
                        excel_writer = writer,
                        sheet_name   = 'description',
                        header       = True,
                        index        = False,
                    )

                    for field in fields_dup:
                        dt_dup[field].to_excel(
                            excel_writer = writer,
                            sheet_name   = "dup_%s" % field,
                            header       = True,
                            index        = False
                        )

        return df_table
