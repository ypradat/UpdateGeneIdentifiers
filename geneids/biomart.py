# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 06 2021
@created: Jan 06 2021

Functions for laoding and updating tables using the biomaRt R package.

Functions
    _get_biomart_table
"""

import pandas as pd
import numpy as np

from rpy2.robjects.packages   import importr
from rpy2.robjects.pandas2ri  import rpy2py_dataframe
import rpy2.robjects          as ro
biomaRt = importr("biomaRt")

DataFrame = pd.core.frame.DataFrame

class Biomart(object):
    def __init__(self, biomart="ensembl", dataset="hsapiens_gene_ensembl", verbose=True):
        """
        Load the selected BioMart database by connecting direclty to an ensembl database withouth specifying the url
        pecify.
        Default mart dataset is the he human genome reference GRCh38.p13.

        Parameters
        ----------
        biomart: `biomart` argument of `biomaRt.useEnsembl`, default is 'ensembl'
        dataset: `dataset` argument of `biomaRt.useEnsembl`, default is 'hsapiens_gene_ensembl'
        verbose:  Set to True for intermediate messages
        """
        self.biomart = biomart
        self.dataset = dataset
        self.verbose = verbose

    def _get_mart(self):
        if self.verbose:
            print("loading biomaRt database...", flush=True)

        mart = biomaRt.useEnsembl(
            biomart = self.biomart,
            dataset = self.dataset,
        )

        if self.verbose:
            print("done!", flush=True)
        return mart

    def _get_bm_table(attributes: list, filters: str, values: list, mart: ro.methods.RS4) -> DataFrame:
        if all([type(e) == int for e in values]):
            r_values = ro.IntVector(values)
        elif all([type(e) == float for e in values]):
            r_values = ro.FloatVector(values)
        else:
            r_values = ro.StrVector(values)

        biomart_table = biomaRt.getBM(
            attributes = ro.StrVector(attributes),
            filters    = filters,
            values     = r_values,
            mart       = mart
        )

        df_biomart_table = rpy2py_dataframe(biomart_table)
        df_biomart_table = df_biomart_table.replace("", np.nan)

        return df_biomart_table

    def _discard_non_standard_chr(df: DataFrame) -> DataFrame:
        standard_chr = ["%d" % i for i in np.arange(1, 23)] + ["X", "Y"]
        mask_keep = df["chromosome_name"].isin(chromosomes)
        df = df.loc[mask_keep]
        if self.verbose:
            print("dropped %d/%d entries where 'chromosome_name' is non standard" % (sum(~mask_keep), len(mask_keep)))
        return df

    def _discard_duplicates_exceptional_rule(self, df, field_b, field_a):
        if field_b == "hgnc_symbol" and field_a == "entrezgene_id":
            # there are cases of entrezgene_id duplicates where one occurrence has no symbol or a symbol with 
            # an hyphen (readthrough transcript)
            mask_dup = df.duplicated(subset=[field_a], keep=False)
            mask_sym = df[field_b].apply(lambda x: True if type(x)==float or "-" in x else False)
            mask = mask_dup & mask_sym
            if self.verbose:
                print("applied entrezgened_id discard duplicates exceptional rule")
                print("dropped %d/%d duplicated entrezgene_id with choice" % (sum(mask), len(mask)))
            return df[~mask]
        else:
            if self.verbose:
                print("no exceptional rule for discarding duplicated")
            return df

    def _discard_duplicates(self, df, field, keep="first"):
        df_nodup = df.drop_duplicates(subset=[field], keep=keep)
        if self.verbose:
            print("dropped %d/%d duplicated %s" % (df.shape[0] - df_nodup.shape[0], df.shape[0], field))

    def _discard_updates_exceptional_rules(self, df, rule):
        if rule == "chromosome_name":
            fields_biomart = [x for x in df.columns if x.endswith("_biomart")]
            mask_chg_chr = df["chromosome_name"] != df["chromosome_name"]
            df.loc[mask_chg_chr, fields_biomart] = np.nan
            if self.verbose:
                print("applied chr exceptional discard updates rule")
                print("dropped %d/%d updates where chr was changed" % (sum(mask_chg_chr), len(mask_chg_chr)))
            return df
        else:
            if self.verbose:
                print("no exceptional rule for discarding updates")
            return df

    def _update_field_b_using_field_a(df, field_b, field_a, mart) ->
    tuple:
        """
        Match input df table with biomart table on field_a in order to update field_b.
        """
        attributes = [field_b, field_a, "chromosome_name", "start_position", "end_position", "band"]
        values_a_in_df = df[field_a].dropna().unique().tolist()
        if self.verbose:
            print("number of unique values of %s in df table: %d", (field_a, len(values_a_in_df)))

        if self.verbose:
            print("\nretrieving %s attributes ..." %s field_a, flush=True)

        attributes = fieldbiomart2attributes[field_a]
        df_biomart = _get_biomart_table(
            attributes = attributes,
            filters    = field_a,
            values     = values_a_in_df,
            mart       = mart,
        )
        if self.verbose:
            print("done!", flush=True)

        df_biomart = self._discard_non_standard_chr(df_biomart)
        df_biomart = self._discard_duplicates_exceptional_rule(df_biomart, field_b, field_a)
        df_biomart = self._discard_duplicates(df_biomart, field=field_a)

        field_a_new = "%s_biomart" % field_a
        field_b_new = "%s_biomart" % field_b
        df_biomart = df_biomart.rename(columns={field_a: field_a_new})

        mask_found = df[field_a].isin(df_biomart[field_a_new])
        df_found = df.loc[mask_found].copy()
        df_not_found = df.loc[~mask_found]
        df_found = df_found.merge(
            right     = df_biomart[[field_a_new, field_b_new]],
            how       = "left",
            left_on   = field_a,
            right_on  = field_a_new
        )
        assert df_found.shape[0] = len(mask_found)


        if self.verbose
            mask_new = df_found[field_b] != df_found[field_b_new]
            svars = (sum(mask_new), len(mask_new), field_b, field_a)
            print("updated %d/%d in %s where %s was found in biomart" % svars)

        return df_found, df_not_found

    def _update_field_using_biomart(self, df, field_to_update, fields_on_biomart=["hgnc_id", "entrezgene_id"]):
        mart = self._get_mart()

        dfs_ok = []
        df_on = df.copy()
        for field_on_biomart in fields_on_biomart:
            df_ok, df_on = self._update_field_b_using_field_a(df_on, field_b=field_to_update, field_a=field_on_biomart,
                                                              mart=mart)
            dfs_ok.append(dfs_ok)
        if df_on.shape[0] > 0:
            dfs_ok.append(df_on)

        df_ok = pd.concat(df_oks, axis=0)

        if self.verbose
            field_new = "%s_biomart" % field_to_update
            mask_new = df_ok[field_to_update] != df_found[field_new]
            print("updated %d/%d values of %s using biomart" % (sum(mask_new), len(mask_new), field_to_update))

        return df_ok

    def update_field_using_biomart(self, df, fields_id_df=["Gene", "Hugo_Symbol", "Chromosome"],
                                   fields_df2biomart = {"Hugo_Symbol": "hgnc_symbol", "HGNC_ID": "hgnc_id",
                                                        "Entrez_Gene_Id": "entrezgene_id",
                                                        "Chromosome": "chromosome_name"},
                                   field_to_update="hgnc_symbol", fields_on_biomart=["hgnc_id", "entrezgene_id"],
                                   fields_extra_after_update=[], rules_discard_updates=["chromosome_name"],
                                   update_mode="add"):

        # define unique row identifier
        df_new = df.copy()
        df_new["Row"] = df_new[fields_id_df].fillna("NaN").agg(lambda x:"/".join(x), axis=1)
        df_new["Row"] = df_new["Row"].sort_values("Row").reset_index(drop=True)

        # match colum names to that of biomart
        df_new = df_new.rename(columns=fields_df2biomart)
        fields_biomart2df = {v:k for k,v in fields_df2biomart.items()}

        # update fields
        df_new = self._update_field_using_biomart(df_new, field_to_update, fields_on_biomart)
        df_new = df_new.sort_values("Row").reset_index(drop=True)

        for rule_discard_update in rules_discard_updates:
            df_new = self._discard_updates_exceptional_rules(df, rule)

        # add fields to original table
        fields_biomart = [x for x in df_new.columns if x.endswith("_biomart")]
        fields_updated = [field_to_update] + fields_extra_after_update
        mask_where_updated = ~df_new["%s_biomart" % field_to_update].isnull()

        df_ori["%s_new" % field] = np.nan
        for field_extra in fields_extra:
            df_ori["%s_new" % field_extra] = np.nan

        # update status field
        mask_updated = ~df_new["%s_biomart" % field_to_update].isnull()
        field_updated = fields_biomart2df[field_to_update]
        df_new.loc[:, "%s_Was_Updated" % field_updated] = False
        df_new.loc[mask_where_updated, "%s_Updated" % field_updated] = True

        if update_mode == "inplace":
            # update inplace
            field_updated_bio = "%s_biomart" % field_updated
            df_new.loc[mask_updated,field_updated] = df_new.loc[mask_updated, field_updated_bio]

        elif update_mode == "add":
            # add updated field
            field_updated_add = "%s_Updated" % field_updated
            field_updated_bio = "%s_biomart" % field_updated

            df_new.loc[:, field_updated_add] = np.nan
            df_new.loc[mask_updated, field_updated_add] = df_new.loc[mask_updated, field_updated_bio]
            df_new.loc[~mask_updated, field_updated_add] = df_new.loc[~mask_updated, field_updated]
        else:
            raise ValueError("Unsupported value of %s for update_mode. Choose one of 'add' or 'replace'" % update_mode)

        # add extra fields (it should be in the attributes added by the biomart table)
        for field_extra in fields_extra_after_update:
            df_new.loc[:,field_extra] = df_new.loc[:, "%s_biomart" % field_extra]

        # delete _biomart fields
        for field_biomart in fields_biomart:
            del df_new[field_biomart]

        return df_new

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
