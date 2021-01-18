# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 18 2021
@created: Jan 18 2021
"""

#     def update_hgnc_using_latest():
#         if field_biomart == "hgnc_symbol":
#             self.verbose:
#                 print("updating hgnc symbols using the latest version from EBI's database...", flush=True)
# 
#             df_new = self._update_hgnc_to_latest(df_new, )




# def _update_table_on(df_table, df_hgnc, cols_table_on, cols_hgnc_on, updates_done: str=None) -> tuple:
#     mask_nu = df_table["HGNC_ID_updated"].isnull()
#     df_table_nu = df_table.loc[mask_nu].copy()
# 
#     if updates_done is None:
#         print("%d/%d rows not yet updated" % (sum(mask_nu), len(mask_nu)))
#     else:
#         print("%d/%d rows not updated with %s" % (sum(mask_nu), len(mask_nu), updates_done))
#     print("trying to fill HGNC_ID_updated on %s ..." % "-".join(cols_hgnc_on))
# 
#     df_table_nu[cols_table_on] = df_table_nu[cols_table_on].fillna(-1)
# 
#     #### merge
#     df_table_mg = df_table_nu.merge(
#         right    = df_hgnc,
#         left_on  = cols_table_on,
#         right_on = cols_hgnc_on,
#         how      = "left"
#     )
# 
#     #### check one-to-one correspondence in merge
#     #### if not, drop duplicates arbitrarily
#     mask_nnl = df_table_mg[cols_table_on].mean(axis=1) == 0
#     mask_dup = df_table_mg.duplicated(subset=["Row"], keep=False)
#     mask_chk = mask_nnl & mask_dup
# 
#     if sum(mask_chk) > 0:
#         print("WARNING: there exist more than one entry in the HGNC table for the following")
#         print(df_table_mg.loc[mask_chk, cols_hgnc_on].drop_duplicates())
#         print("choose the symbol from the first in each set of duplicates")
#         df_table_mg = df_table_mg.drop_duplicates(subset=["Row"], keep="first")
# 
#     assert df_table_mg["Row"].equals(df_table.loc[mask_nu, "Row"].reset_index(drop=True))
# 
#     #### update columns
#     df_table.loc[mask_nu,"HGNC_ID_updated"]       = df_table_mg["hgnc_id"].values
#     df_table.loc[mask_nu,"HGNC_Symbol_updated"]   = df_table_mg["symbol"].values
#     df_table.loc[mask_nu,"Gene_HGNC_updated"]     = df_table_mg["ensembl_gene_id"].values
#     df_table.loc[mask_nu,"location_HGNC_updated"] = df_table_mg["location"].values
# 
#     if updates_done is None:
#         updates_done = "%s" % "-".join(cols_hgnc_on)
#     else:
#         updates_done += ";%s" % "-".join(cols_hgnc_on)
# 
#     return df_table, updates_done
# 
# def _update_table_hgnc(df_table: DataFrame, dt_cols: str) -> DataFrame:
#     #### downloaded from ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt 
#     #### redownload for updated table
# 
#     current_wd = setwd_to_data()
#     filepath   = "./raw/div/gene_symbols/hgnc_complete_set.txt"
#     df_hgnc = pd.read_csv(
#         filepath_or_buffer = filepath,
#         sep                = "\t",
#         low_memory         = False
#     )
#     os.chdir(current_wd)
# 
#     extract_chr = lambda x: x if type(x)==float else re.split("p|q", x)[0]
#     df_hgnc.insert(0, "chromosome", df_hgnc["location"].apply(extract_chr))
# 
#     #### init columns with updates values
#     df_table.insert(df_table.shape[1], "HGNC_ID_updated", np.nan)
#     df_table.insert(df_table.shape[1], "HGNC_Symbol_updated", np.nan)
#     df_table.insert(df_table.shape[1], "Gene_HGNC_updated", np.nan)
#     df_table.insert(df_table.shape[1], "location_HGNC_updated", np.nan)
#     updates_done = None
# 
#     #### try updating with HGNC_ID and HGNC_ID_biomart
#     for col_hgnc_id in dt_cols["hgnc_id"]:
#         df_table, updates_done = _update_table_on(
#             df_table      = df_table,
#             df_hgnc       = df_hgnc,
#             cols_table_on = [col_hgnc_id],
#             cols_hgnc_on  = ["hgnc_id"],
#             updates_done  = updates_done
#         )
# 
#     #### try with ensembl_id
#     df_table, updates_done = _update_table_on(
#         df_table      = df_table,
#         df_hgnc       = df_hgnc,
#         cols_table_on = [dt_cols["ensembl_id"]],
#         cols_hgnc_on  = ["ensembl_gene_id"],
#         updates_done  = updates_done
#     )
# 
#     #### try with current symbol 
#     df_table, updates_done = _update_table_on(
#         df_table      = df_table,
#         df_hgnc       = df_hgnc,
#         cols_table_on = [dt_cols["symbol"], "Chromosome"],
#         cols_hgnc_on  = ["symbol", "chromosome"],
#         updates_done  = updates_done
#     )
# 
#     #### try with alias symbol 
#     df_table, updates_done = _update_table_on(
#         df_table      = df_table,
#         df_hgnc       = df_hgnc,
#         cols_table_on = [dt_cols["symbol"], "Chromosome"],
#         cols_hgnc_on  = ["alias_symbol", "chromosome"],
#         updates_done  = updates_done
#     )
# 
#     #### try with prev symbol
#     df_hgnc_explode = explode(
#         df = df_hgnc,
#         cols = ["prev_symbol"],
#         sep = "|"
#     )
# 
#     df_table, updates_done = _update_table_on(
#         df_table      = df_table,
#         df_hgnc       = df_hgnc_explode,
#         cols_table_on = [dt_cols["symbol"], "Chromosome"],
#         cols_hgnc_on  = ["prev_symbol", "chromosome"],
#         updates_done  = updates_done
#     )
# 
#     #### correct incorrect matches
#     #### if chromosome number is different drop update
#     df_table["Chromosome_HGNC_updated"] = df_table["location_HGNC_updated"].apply(extract_chr)
# 
#     mask_incorrect = df_table["Chromosome"] != df_table["Chromosome_HGNC_updated"]
#     mask_incorrect = mask_incorrect & (df_table[["Chromosome", "Chromosome_HGNC_updated"]].isnull().sum(axis=1) == 0)
#     print("%d/%d incorrect matches" % (sum(mask_incorrect), len(mask_incorrect)))
#     df_table.loc[mask_incorrect, "HGNC_ID_updated"] = np.nan
#     df_table.loc[mask_incorrect, "HGNC_Symbol_updated"] = np.nan
#     df_table.loc[mask_incorrect, "Gene_HGNC_updated"] = np.nan
#     df_table.loc[mask_incorrect, "location_HGNC_updated"] = np.nan
#     df_table.loc[mask_incorrect, "Chromosome_HGNC_updated"] = np.nan
# 
#     mask_not_found = df_table["HGNC_ID_updated"].isnull()
#     print("%d/%d rows not updated" % (sum(mask_not_found), len(mask_not_found)))
# 
#     return df_table
# 
# def get_table_update_symbols(df_maf: DataFrame, col_gene_symbol: str) -> DataFrame:
#     df_table_original = get_table_gene_identifiers(
#         df_maf     = df_maf,
#         fields     = ["Gene", "Hugo_Symbol", "Chromosome", "SYMBOL_SOURCE", "HGNC_ID", "Entrez_Gene_Id"],
#         fields_dup = None
#     )
#     cols_original = df_table_original.columns
# 
#     #### drop duplicates i.e genes with same Gene - Hugo_Symbol - Chromosome combination
#     #### but different HGNC_ID (one is NaN for instance) or different Entrez_Gene_Id (one is 0 for instance).
#     #### keep the duplicate with the most info
#     info_hgnc = (~df_table_original.HGNC_ID.isnull()).astype(int)
#     info_entrez = (df_table_original.Entrez_Gene_Id != 0).astype(int)
#     df_table_original["info"] = info_hgnc + info_entrez
#     df_table_original = df_table_original.sort_values(by="info")
#     df_table_original = df_table_original.drop_duplicates(subset=["Gene", "Hugo_Symbol", "Chromosome"], keep="last")
#     del df_table_original["info"]
# 
#     cols_unique = ["Gene", "Hugo_Symbol", "Chromosome"]
#     df_table_original["Row"] = df_table_original[cols_unique].fillna("NaN").agg(lambda x:"/".join(x), axis=1)
#     df_table_original = df_table_original.sort_values("Row").reset_index(drop=True)
# 
#     #### add prefix HGNC: for matching ids from other tables
#     df_table_original["HGNC_ID"] = df_table_original["HGNC_ID"].apply(lambda x: x if np.isnan(x) else "HGNC:%d" % x)
# 
#     #### update hgnc id and hugo symbols using biomaRt from Ensembl
#     df_table_biomart = _update_table_biomart(df_table_original)
#     df_table_biomart = df_table_biomart.sort_values("Row").reset_index(drop=True)
# 
#     #### update HGNC symbols using the latest verions of hgnc_complete_set from EBI's database
#     #### update biomart HGNC symbols and original symbols in case the gene was not found in biomart
#     df_table_biomart_hgnc = _update_table_hgnc(
#         df_table = df_table_biomart,
#         dt_cols  = {"hgnc_id": ["HGNC_ID_biomart", "HGNC_ID"], "ensembl_id": "Gene", "symbol": "Hugo_Symbol"}
#     )
#     df_table_biomart_hgnc = df_table_biomart_hgnc.sort_values("Row").reset_index(drop=True)
# 
#     mask = df_table_biomart_hgnc["HGNC_Symbol_updated"].isnull() & (df_table_biomart_hgnc["SYMBOL_SOURCE"] == "HGNC")
#     if sum(mask) > 0:
#         print("\nWARNING: %d genes with SYMBOL_SOURCE HGNC could not be found anywhere" % sum(mask))
#         print(df_table_biomart_hgnc[mask].iloc[:,:6])
# 
#     #### assemble new data 
#     df_table_original["Hugo_Symbol_new"] = np.nan
#     df_table_original["HGNC_ID_new"] = np.nan
#     df_table_original["Gene_new"] = np.nan
#     df_table_original["Band_new"] = np.nan
# 
#     mask_new = ~df_table_biomart_hgnc.HGNC_Symbol_updated.isnull()
# 
#     df_table_original.loc[mask_new, "Hugo_Symbol_new"] = df_table_biomart_hgnc.loc[mask_new, "HGNC_Symbol_updated"]
#     df_table_original.loc[mask_new, "HGNC_ID_new"] = df_table_biomart_hgnc.loc[mask_new, "HGNC_ID_updated"]
#     df_table_original.loc[mask_new, "Gene_new"] = df_table_biomart_hgnc.loc[mask_new, "Gene_HGNC_updated"]
#     df_table_original.loc[mask_new, "Band_new"] = df_table_biomart_hgnc.loc[mask_new, "location_HGNC_updated"]
# 
#     df_table_original.loc[~mask_new, "Hugo_Symbol_new"] = df_table_biomart_hgnc.loc[~mask_new, "Hugo_Symbol"]
#     df_table_original.loc[~mask_new, "HGNC_ID_new"] = df_table_biomart_hgnc.loc[~mask_new, "HGNC_ID"]
#     df_table_original.loc[~mask_new, "Gene_new"] = df_table_biomart_hgnc.loc[~mask_new, "Gene"]
# 
#     return df_table_original
# 
# def update_gene_symbols(df_maf: DataFrame, col_gene_symbol: str="Hugo_Symbol") -> tuple:
#     df_table_update = get_table_update_symbols(
#         df_maf          = df_maf,
#         col_gene_symbol = col_gene_symbol
#     )
# 
#     cols_unique = ["Gene", "Hugo_Symbol", "Chromosome"]
#     df_maf["Row"] = df_maf[cols_unique].fillna("NaN").agg(lambda x:"/".join(x), axis=1)
# 
#     #### merge new Hugo_Symbol with df_maf
#     df_maf = df_maf.merge(
#         right = df_table_update[["Row", "Hugo_Symbol_new", "HGNC_ID_new", "Gene_new", "Band_new"]],
#         how   = "left",
#         on    = "Row"
#     )
# 
#     del df_maf["Row"]
# 
#     return df_maf, df_table_update

