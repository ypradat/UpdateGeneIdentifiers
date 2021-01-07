"""
The :pkg:`geneids` package provides functions for establishing correspondences between different gene identifiers and
for updating hgnc symbols.
"""

from .make_gene_table import make_gene_table
from .utils import discard_rows_with_less_info

__all__ = [
    'discard_rows_with_less_info',
    'make_gene_table',
]
