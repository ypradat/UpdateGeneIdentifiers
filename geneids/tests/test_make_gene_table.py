# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 05 2021
@created: Jan 05 2021

Test for make_gene_table module.
"""

from geneids import make_gene_table

def test_make_gene_table(sample_maf):
    table = make_gene_table(sample_maf,
                            fields_ids=["Gene", "Hugo_Symbol", "Chromosome", "HGNC_ID", "Entrez_Gene_Id"],
                            fields_dup=["Gene", "Hugo_Symbol", "HGNC_ID", "Entrez_Gene_Id"],
                            save_per_field_dup=True,
                            save_folder="examples/output")
