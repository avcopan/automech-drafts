"""Primary mechanism processing routines
"""

import itertools

import automol
import pandas
from tqdm import tqdm

from mecha import schema


def combine_duplicates(
    rxn_df: pandas.DataFrame, first: bool = False
) -> pandas.DataFrame:
    """Combine duplicate reactions, so each reaction only appears once

    :param rxn_df: The reactions dataframe
    :param first: Keep only the first column value? Otherwise, the values will be
        grouped into lists
    :return: The reactions dataframe
    """
    agg_ = "first" if first else tuple
    rxn_df = rxn_df.groupby(["reactants", "products"], as_index=False).agg(agg_)

    if not first:
        rxn_df = rename_with_plural_columns(rxn_df)

    return schema.validate_reactions(rxn_df)


def expand_duplicates(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    """Expand duplicate reactions, so they appear multiple times

    :param rxn_df: The reactions dataframe
    :return: The reactions dataframe
    """
    all_exp_cols = [schema.Reactions.rates, schema.Reactions.steps]
    exp_cols = [c for c in all_exp_cols if c in rxn_df]
    rxn_df = rxn_df.explode(exp_cols)
    return rename_with_singular_columns(rxn_df)


def classify_reactions(
    rxn_df: pandas.DataFrame, spc_df: pandas.DataFrame, expand: bool = True
) -> pandas.DataFrame:
    """Classify the reactions in a mechanism

    :param rxn_df: The reactions dataframe
    :param spc_df: The species dataframe
    :param expand: Expand reactions with multiple possible steps?
    :return: The reactions dataframe, with reaction objects for classified reactions
    """
    rxn_df = schema.validate_reactions(rxn_df)
    spc_df = schema.validate_species(spc_df)

    rxn_df = combine_duplicates(rxn_df, first=False)
    rxn_df.drop(columns=[schema.Reactions.rates], inplace=True)

    name2chi = dict(zip(spc_df["name"], spc_df["chi"]))

    def steps_(rcts, prds):
        rchis = list(map(name2chi.get, rcts))
        pchis = list(map(name2chi.get, prds))
        try:
            return automol.reac.from_amchis(rchis, pchis, stereo=False)
        except Exception as exc:
            return exc

    # Do the classification with a progress bar, since it may take a while
    rxns_lst = list(
        zip(rxn_df[schema.Reactions.reactants], rxn_df[schema.Reactions.products])
    )
    rxn_df[schema.Reactions.steps] = list(itertools.starmap(steps_, tqdm(rxns_lst)))

    # Expand reactions with multiple possible steps, if requested
    if expand:
        rxn_df = expand_duplicates(rxn_df)

    return schema.validate_reactions(rxn_df)


# Helpers
def rename_with_plural_columns(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    return rxn_df.rename(
        columns={
            schema.Reactions.rate: schema.Reactions.rates,
            schema.Reactions.step: schema.Reactions.steps,
        }
    )


def rename_with_singular_columns(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    return rxn_df.rename(
        columns={
            schema.Reactions.rates: schema.Reactions.rate,
            schema.Reactions.steps: schema.Reactions.step,
        }
    )
