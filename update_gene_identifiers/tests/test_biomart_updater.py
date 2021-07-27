# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 05 2021
@created: Jan 05 2021

Test for updater module
"""

import numpy as np

from update_gene_identifiers import discard_rows_with_less_info
from update_gene_identifiers import TableMaker
from update_gene_identifiers import BiomartUpdater

def test_Updater(sample_maf):
    table_maker = TableMaker(fields_ids=["Gene", "Hugo_Symbol", "Chromosome", "HGNC_ID", "Entrez_Gene_Id"],
                             fields_dup=["Gene", "Hugo_Symbol", "HGNC_ID", "Entrez_Gene_Id"],
                             save_per_field_dup=True,
                             save_folder="examples/output")

    table = table_maker.make(sample_maf)

    table = discard_rows_with_less_info(table,
                                        fields_subset=["Gene", "Hugo_Symbol", "Chromosome"],
                                        fields_info=["HGNC_ID", "Entrez_Gene_Id"])


    table["Entrez_Gene_Id"] = table["Entrez_Gene_Id"].replace(0, np.nan)
    fields_id = ["Gene", "Hugo_Symbol", "Entrez_Gene_Id", "HGNC_ID"]
    table = table.loc[table[fields_id].isnull().mean(axis=1) < 1]
    table["HGNC_ID"] = table["HGNC_ID"].apply(lambda x: x if np.isnan(x) else "HGNC:%d" % x)
    table = table.reset_index(drop=True)

    updater = BiomartUpdater()
    table_updated = updater.update_field(df=table,
                                         df2bio={"Hugo_Symbol": "hgnc_symbol", "HGNC_ID": "hgnc_id", "Entrez_Gene_Id":
                                                 "entrezgene_id", "Chromosome": "chromosome_name"},
                                         field_to_update="hgnc_symbol",
                                         fields_on_biomart=["hgnc_id", "entrezgene_id"],
                                         rules_discard_updates={"chromosome_name": "chromosome_name"},
                                         update_mode="add",
                                         suffix="Updated_Biomart")


    updater = BiomartUpdater(verbose=False)
    table_updated = updater.update_field(df=table,
                                         df2bio={"Hugo_Symbol": "hgnc_symbol", "HGNC_ID": "hgnc_id", "Entrez_Gene_Id":
                                                 "entrezgene_id", "Chromosome": "chromosome_name"},
                                         field_to_update="hgnc_symbol",
                                         fields_on_biomart=["hgnc_id", "entrezgene_id"],
                                         rules_discard_updates={"chromosome_name": "chromosome_name"},
                                         update_mode="add",
                                         suffix="Updated_Biomart")




    # df_table_original["HGNC_ID"] = df_table_original["HGNC_ID"].apply(lambda x: x if np.isnan(x) else "HGNC:%d" % x)

    # cols_unique = ["Gene", "Hugo_Symbol", "Chromosome"]
    # df_table_original["Row"] = df_table_original[cols_unique].fillna("NaN").agg(lambda x:"/".join(x), axis=1)
    # df_table_original = df_table_original.sort_values("Row").reset_index(drop=True)
