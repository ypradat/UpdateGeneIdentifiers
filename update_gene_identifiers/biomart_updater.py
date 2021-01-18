# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat
@modified: Jan 18 2021
@created: Jan 06 2021

Class for updating tables using the biomaRt R package.
"""

import pandas as pd
import numpy as np

from rpy2.robjects.packages   import importr
from rpy2.robjects.pandas2ri  import rpy2py_dataframe
import rpy2.robjects          as ro
biomaRt = importr("biomaRt")

DataFrame = pd.core.frame.DataFrame

class BiomartUpdater(object):
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
            mirror  = "www",
            verbose = self.verbose,
        )

        if self.verbose:
            print("done!", flush=True)
        return mart

    def _get_bm_table(self, attributes: list, filters: str, values: list, mart: ro.methods.RS4) -> DataFrame:
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

    def _discard_non_standard_chr(self, df: DataFrame) -> DataFrame:
        standard_chr = ["%d" % i for i in np.arange(1, 23)] + ["X", "Y"]
        mask_keep = df["chromosome_name"].isin(standard_chr)
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
        return df_nodup

    def _discard_updates_exceptional_rules(self, df, field_to_update, field_on_biomart, rule):
        suffix = "_biomart_from_%s" % field_on_biomart
        if rule == "chromosome_name":
            fields_biomart = [x for x in df.columns if x.endswith(suffix)]
            mask_null_update = df["%s%s" % (field_to_update, suffix)].isnull()
            mask_chg_chr = (df["chromosome_name"] != df["chromosome_name%s" % suffix]) & ~mask_null_update
            df.loc[mask_chg_chr, fields_biomart] = np.nan
            if self.verbose:
                print("applied chr exceptional discard updates rule")
                print("dropped %d/%d updates where chr was changed" % (sum(mask_chg_chr), len(mask_chg_chr)))
            return df
        else:
            if self.verbose:
                print("no exceptional rule for discarding updates")
            return df

    def _update_field(self, df, field_to_update, field_on_biomart, fields_extra, mart):
        attributes = [field_to_update, field_on_biomart, "chromosome_name", "start_position", "end_position", "band"]
        values_on_biomart_in_df = df[field_on_biomart].dropna().unique().tolist()
        if self.verbose:
            print("there are %d unique values of %s in df table" % (len(values_on_biomart_in_df), field_on_biomart))

        if self.verbose:
            print("retrieving %s attributes from biomart db ..." % field_on_biomart, flush=True)

        df_biomart = self._get_bm_table(
            attributes = attributes,
            filters    = field_on_biomart,
            values     = values_on_biomart_in_df,
            mart       = mart,
        )

        if self.verbose:
            print("retrieving %s attributes from biomart db ... done !" % field_on_biomart, flush=True)

        df_biomart = self._discard_non_standard_chr(df_biomart)
        df_biomart = self._discard_duplicates_exceptional_rule(df_biomart, field_to_update, field_on_biomart)
        df_biomart = self._discard_duplicates(df_biomart, field=field_on_biomart)

        field_on_biomart_new = "%s_biomart" % field_on_biomart
        field_to_update_new = "%s_biomart_from_%s" % (field_to_update, field_on_biomart)
        df_biomart = df_biomart.rename(columns={field_on_biomart: field_on_biomart_new, field_to_update:
                                                field_to_update_new})

        fields_extra_new = ["%s_biomart_from_%s" % (field, field_on_biomart) for field in fields_extra]
        df_biomart = df_biomart.rename(columns={k:v for k,v in zip(fields_extra, fields_extra_new)})

        mask_found = df[field_on_biomart].isin(df_biomart[field_on_biomart_new])
        df_found = df.loc[mask_found].copy()
        df_not_found = df.loc[~mask_found]
        df_found = df_found.merge(
            right     = df_biomart[[field_on_biomart_new, field_to_update_new] + fields_extra_new],
            how       = "left",
            left_on   = field_on_biomart,
            right_on  = field_on_biomart_new
        )

        if not df_found.shape[0] == sum(mask_found):
            raise Exception("The merge of the input table and biomart table resulted in duplicates")

        if self.verbose:
            print("found %d/%d values of %s in biomart" % (sum(mask_found), len(mask_found), field_on_biomart))
            mask_new = df_found[field_to_update] != df_found[field_to_update_new]
            svars = (sum(mask_new), len(mask_new), field_to_update, field_on_biomart)
            print("updated %d/%d values of %s where %s was found in biomart" % svars)

        df_updated = pd.concat((df_found, df_not_found), axis=0)
        df_updated = df_updated.sort_values("Index")
        del df_updated[field_on_biomart_new]

        return df_updated

    def _check_hgnc_id_format(self, df, field):
        if not all(df[field].apply(lambda x: (type(x)==float and np.isnan(x)) or (type(x)==str and x.startswith("HGNC:")))):
            raise ValueError("Values for the field %s should be null or start with 'HGNC:'" % field)


    def _check_id_formats(self, df, df2bio):
        for id_df, id_bio in df2bio.items():
            if id_bio == "hgnc_id":
                self._check_hgnc_id_format(df, id_df)

    def update_field(self, df, df2bio={"Hugo_Symbol": "hgnc_symbol", "HGNC_ID": "hgnc_id",
                                       "Entrez_Gene_Id": "entrezgene_id", "Chromosome": "chromosome_name"},
                     field_to_update="hgnc_symbol", fields_on_biomart=["hgnc_id", "entrezgene_id"],
                     fields_extra=[], rules_discard_updates={"chromosome_name": "chromosome_name"},
                     update_mode="add", suffix="Updated_Biomart"):

        # connect to biomart db
        mart = self._get_mart()

        # define unique row identifier
        df_new = df.copy()
        df_new["Index"] = np.arange(df_new.shape[0])
        df_new = df_new.sort_values("Index")

        # match colum names to that of biomart
        df_new = df_new.rename(columns=df2bio)
        bio2df = {v:k for k,v in df2bio.items()}

        # check format of some identifiers
        self._check_id_formats(df, df2bio)

        # update fields
        fields_discard_rules = list(set(rules_discard_updates.keys()))
        for field_on_biomart in fields_on_biomart:
            if self.verbose:
                print("="*40)
            df_new = self._update_field(df_new, field_to_update, field_on_biomart, fields_discard_rules, mart)
            if self.verbose:
                print("="*40)

            # discard update were required
            for rule_discard_update in rules_discard_updates.values():
                df_new = self._discard_updates_exceptional_rules(df_new, field_to_update, field_on_biomart,
                                                                 rule_discard_update)

        # merge updates
        def _merge_updates(x):
            u = list(set(x.tolist()))
            if np.nan in u:
                u.remove(np.nan)
            if len(u) == 0:
                return np.nan
            elif len(u) == 1:
                return u[0]
            else:
                return "CONFLICTING UPDATES"

        field_to_update_news = ["%s_biomart_from_%s" % (field_to_update, x) for x in fields_on_biomart]
        field_to_update_new = "%s_biomart_updated" % field_to_update
        df_new.loc[:, field_to_update_new] = df_new[field_to_update_news].apply(_merge_updates, axis=1)

        # add fields to original table
        fields_biomart = [x for x in df_new.columns if "_biomart" in x]
        mask_updated = (~df_new[field_to_update_new].isnull()) & (df_new[field_to_update]!=df_new[field_to_update_new])

        # update status field
        field_updated = bio2df[field_to_update]

        if update_mode == "inplace":
            # update inplace
            df_new.loc[mask_updated,field_updated] = df_new.loc[mask_updated, field_to_update_new]

        elif update_mode == "add":
            # add updated field
            field_updated_add = "%s_%s" % (field_updated, suffix)
            df_new.loc[:, field_updated_add] = np.nan
            df_new.loc[mask_updated, field_updated_add] = df_new.loc[mask_updated, field_to_update_new]
            df_new.loc[~mask_updated, field_updated_add] = df_new.loc[~mask_updated, field_to_update]
        else:
            raise ValueError("Unsupported value of %s for update_mode. Choose one of 'add' or 'replace'" % update_mode)

        # add extra fields (it should be in the attributes added by the biomart table)
        for field_extra in fields_extra:
            field_extra_news = ["%s_biomart_from_%s" % (field_extra, x) for x in fields_on_biomart]
            field_extra_new = "%s_Biomart" % bio2df[field_extra]
            df_new.loc[:, field_extra_new] = df_new[field_extra_news].apply(_merge_updates, axis=1)

        # delete _biomart fields
        for field_biomart in fields_biomart:
            del df_new[field_biomart]

        del df_new["Index"]
        df_new = df_new.rename(columns=bio2df)

        return df_new
