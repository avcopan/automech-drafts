"""Primary mechanism processing routines
"""

import automol
import pandas
from tqdm.auto import tqdm

from mecha import data, schema, species
from mecha.schema import Reactions, Species
from mecha.util import df_

tqdm.pandas()


def classify_reactions(
    rxn_df: pandas.DataFrame, spc_df: pandas.DataFrame
) -> pandas.DataFrame:
    """Classify the reactions in a mechanism

    :param rxn_df: The reactions dataframe
    :param spc_df: The species dataframe
    :return: The reactions dataframe, with reaction objects for classified reactions
    """
    rxn_df = schema.validate_reactions(rxn_df)
    spc_df = schema.validate_species(spc_df)

    rxn_df = combine_duplicates(rxn_df, first=False)

    # Do the classification with a progress bar, since it may take a while
    chi_dct = df_.lookup_dict(spc_df, Species.name, Species.chi)

    def objs_(eq):
        rchis, pchis = data.reac.parse_equation(eq, trans_dct=chi_dct)
        objs = automol.reac.from_amchis(rchis, pchis, stereo=False)
        return objs if objs else pandas.NA

    rxn_df[Reactions.obj] = rxn_df[Reactions.eq].progress_apply(objs_)

    # Separate out the unclassified reactions
    unc_rxn_df = rxn_df[rxn_df[Reactions.obj].isna()].drop(columns=[Reactions.obj])
    rxn_df = rxn_df[rxn_df[Reactions.obj].notna()].drop(columns=[Reactions.rate])

    # Expand reactions with multiple possible mechanisms
    def amchi_(obj):
        return automol.reac.ts_amchi(obj)

    rxn_df = expand_duplicates(rxn_df)
    rxn_df[Reactions.chi] = rxn_df[Reactions.obj].progress_apply(amchi_)

    # Expand duplicates among the unclassified reactions again
    unc_rxn_df = expand_duplicates(unc_rxn_df)

    return schema.validate_reactions(rxn_df), schema.validate_reactions(unc_rxn_df)


def expand_stereo(
    rxn_df: pandas.DataFrame,
    spc_df: pandas.DataFrame,
    enant: bool = True,
    expand_species: bool = True,
) -> pandas.DataFrame:
    """Classify the reactions in a mechanism

    :param rxn_df: The reactions dataframe
    :param spc_df: The species dataframe
    :param enant: Distinguish between enantiomers?, defaults to True
    :param expand_species: Expand the species dataframe?
        If set to False, the species dataframe must already be expanded
    :return: The reactions dataframe, with reaction objects for classified reactions
    """

    # Expand species, if requested
    if expand_species:
        spc_df = species.expand_stereo(spc_df, enant=enant)

    if not enant:
        raise NotImplementedError("Reduced expansion not yet implemented")

    rxn_df = schema.validate_reactions(rxn_df)
    rxn_df = rename_with_original_columns(rxn_df)

    # Expand reaction objects
    def objs_(obj):
        if isinstance(obj, str):
            obj = automol.reac.from_string(obj)
        return automol.reac.expand_stereo(obj, enant=enant)

    rxn_df[Reactions.obj] = rxn_df[Reactions.obj].progress_apply(objs_)
    rxn_df = rxn_df.explode(Reactions.obj)

    # Get stereo-resolved reaction equations
    name_dct = df_.lookup_dict(spc_df, [Species.orig_name, Species.chi], Species.name)

    def eq_(row):
        rname0s, pname0s = data.reac.parse_equation(row[Reactions.orig_eq])
        rchis, pchis = automol.reac.amchis(row[Reactions.obj])
        rnames = tuple(map(name_dct.get, zip(rname0s, rchis)))
        pnames = tuple(map(name_dct.get, zip(pname0s, pchis)))
        # Make sure we print the ChIs if something went wrong
        assert all(isinstance(n, str) for n in rnames + pnames), f"{rchis} = {pchis}"
        return data.reac.form_equation(rnames, pnames)

    rxn_df[Reactions.eq] = rxn_df.progress_apply(eq_, axis=1)

    return schema.validate_reactions(rxn_df)


# Helpers
def combine_duplicates(
    rxn_df: pandas.DataFrame, first: bool = False
) -> pandas.DataFrame:
    """Combine duplicate reactions, so each reaction only appears once

    :param rxn_df: The reactions dataframe
    :param first: Keep only the first column value? Otherwise, the values will be
        grouped into lists
    :return: The reactions dataframe
    """
    agg_ = "first" if first else list
    agg_dct = {c: agg_ for c in schema.DUP_DIFF_COLS if c in rxn_df}
    agg_dct.update(
        {c: "first" for c in rxn_df.columns if c not in schema.DUP_DIFF_COLS}
    )
    return rxn_df.groupby(Reactions.eq, as_index=False).agg(agg_dct)


def expand_duplicates(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    """Expand duplicate reactions, so they appear multiple times

    :param rxn_df: The reactions dataframe
    :return: The reactions dataframe
    """
    exp_cols = [c for c in schema.DUP_DIFF_COLS if c in rxn_df]
    rxn_df[exp_cols] = rxn_df[exp_cols].map(list)
    return rxn_df.explode(exp_cols)


def rename_with_original_columns(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    """Rename the reactions dataframe with original column names
    (orig_reactants, orig_products, orig_chi)

    Used when expanding stereo

    :param rxn_df: The reactions dataframe
    :return: The reactions dataframe
    """
    return rxn_df.rename(columns=dict(zip(schema.R_CURR_COLS, schema.R_ORIG_COLS)))
