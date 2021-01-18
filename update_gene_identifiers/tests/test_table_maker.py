# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 12 2021
@created: Jan 05 2021

Test for tablemaker module.
"""

from update_gene_identifiers import TableMaker

def test_TableMaker(sample_maf):
    table_maker = TableMaker(fields_ids=["Gene", "Hugo_Symbol", "Chromosome", "HGNC_ID", "Entrez_Gene_Id"],
                             fields_dup=["Gene", "Hugo_Symbol", "HGNC_ID", "Entrez_Gene_Id"],
                             save_per_field_dup=True,
                             save_folder="examples/output")

    table = table_maker.make(sample_maf)
