"""
The :pkg:`update_gene_identifiers` package provides functions for establishing correspondences between different gene identifiers and
for updating hgnc symbols.
"""

from .table_maker import TableMaker
from .utils import discard_rows_with_less_info
from .biomart_updater import BiomartUpdater

__all__ = [
    'discard_rows_with_less_info',
    'TableMaker',
]
