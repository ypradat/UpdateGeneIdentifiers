# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 05 2021
@created: Jan 05 2021

Test for make_gene_table module.
"""

from geneids import discard_rows_with_less_info
from geneids import make_gene_table

def test_make_gene_table(sample_maf):
    table = make_gene_table(sample_maf,
                            fields_ids=["Gene", "Hugo_Symbol", "Chromosome", "HGNC_ID", "Entrez_Gene_Id"],
                            fields_dup=["Gene", "Hugo_Symbol", "HGNC_ID", "Entrez_Gene_Id"],
                            save_per_field_dup=True,
                            save_folder="examples/output")

    table = discard_rows_with_less_info(table, fields_subset=["Gene", "Hugo_Symbol", "Chromosome"],
                                        fields_info=["HGNC_ID", "Entrez_Gene_Id"])

    df_table_original["HGNC_ID"] = df_table_original["HGNC_ID"].apply(lambda x: x if np.isnan(x) else "HGNC:%d" % x)

    cols_unique = ["Gene", "Hugo_Symbol", "Chromosome"]
    df_table_original["Row"] = df_table_original[cols_unique].fillna("NaN").agg(lambda x:"/".join(x), axis=1)
    df_table_original = df_table_original.sort_values("Row").reset_index(drop=True)

    #### add prefix HGNC: for matching ids from other tables
