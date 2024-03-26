"""Primary mechanism processing routines
"""

from typing import Tuple

import automol
import pandas
from tqdm.auto import tqdm

from mecha import data, schema, species
from mecha.schema import Reactions, Species
from mecha.util import df_

tqdm.pandas()


def display_reactions(
    rxn_df: pandas.DataFrame,
    spc_df: pandas.DataFrame,
    keys: Tuple[str, ...] = (Reactions.eq,),
    stereo: bool = True,
):
    """Display the reactions in a mechanism

    :param rxn_df: The reactions dataframe
    :param spc_df: The species dataframe
    :param keys: Keys of extra columns to print
    :param stereo: Display with stereochemistry?
    """
    rxn_df = schema.validate_reactions(rxn_df)
    spc_df = schema.validate_species(spc_df)

    # Do the classification with a progress bar, since it may take a while
    chi_dct = df_.lookup_dict(spc_df, Species.name, Species.chi)

    def display_(row):
        # Print the requested information
        for key in keys:
            val = row[key] if key in row else "No such column"
            print(f"{key}: {val}")

        # Display the reaction
        eq = row[Reactions.eq]
        rchis, pchis, _ = data.reac.read_chemkin_equation(eq, trans_dct=chi_dct)
        if not all(isinstance(n, str) for n in rchis + pchis):
            print(f"Some ChIs missing from species table: {rchis} = {pchis}")
        else:
            automol.amchi.display_reaction(rchis, pchis, stereo=stereo)

    rxn_df.apply(display_, axis=1)


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
        rchis, pchis, _ = data.reac.read_chemkin_equation(eq, trans_dct=chi_dct)
        objs = automol.reac.from_amchis(rchis, pchis, stereo=False)
        return objs if objs else pandas.NA

    rxn_df[Reactions.obj] = rxn_df[Reactions.eq].progress_apply(objs_)

    # Separate out the unclassified reactions
    err_df = rxn_df[rxn_df[Reactions.obj].isna()].drop(columns=[Reactions.obj])
    rxn_df = rxn_df[rxn_df[Reactions.obj].notna()].drop(columns=[Reactions.rate])

    # Expand reactions with multiple possible mechanisms
    def amchi_(obj):
        return automol.reac.ts_amchi(obj)

    rxn_df = expand_duplicates(rxn_df)
    rxn_df[Reactions.chi] = rxn_df[Reactions.obj].progress_apply(amchi_)

    # Expand duplicates among the unclassified reactions again
    err_df = expand_duplicates(err_df)

    return schema.validate_reactions(rxn_df), schema.validate_reactions(err_df)


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

    # Get stereo-resolved reaction equations (initially as lists)
    name_dct = df_.lookup_dict(spc_df, [Species.orig_name, Species.chi], Species.name)

    def eq_(row):
        eq0 = row[Reactions.orig_eq]
        objs = row[Reactions.obj]
        eqs = []
        for obj in objs:
            rname0s, pname0s, _ = data.reac.read_chemkin_equation(eq0)
            rchis, pchis = automol.reac.amchis(obj)
            rnames = tuple(map(name_dct.get, zip(rname0s, rchis)))
            pnames = tuple(map(name_dct.get, zip(pname0s, pchis)))
            if not all(isinstance(n, str) for n in rnames + pnames):
                return pandas.NA
            eqs.append(data.reac.write_chemkin_equation(rnames, pnames))
        return tuple(eqs)

    rxn_df[Reactions.eq] = rxn_df.progress_apply(eq_, axis=1)

    # Separate out the reactions that had errors
    err_df = rxn_df[rxn_df[Reactions.eq].isna()]
    rxn_df = rxn_df[rxn_df[Reactions.eq].notna()]

    err_df = err_df.drop(columns=[Reactions.eq]).rename(
        columns={Reactions.orig_eq: Reactions.eq}
    )
    rxn_df = rxn_df.explode([Reactions.eq, Reactions.obj])

    return schema.validate_reactions(rxn_df), schema.validate_reactions(err_df)


def to_mechanalyzer(
    rxn_df: pandas.DataFrame, spc_df: pandas.DataFrame, drop: bool = False
) -> Tuple[str, str]:
    """Get old reaction and species dictionaries for a mechanism

    :param rxn_df: The reactions dataframe
    :param spc_df: The species dataframe
    :param drop: Drop species which have no reactions?, defaults to False
    :return: The mechanism and species list, in CHEMKIN and CSV formats, respectively
    """
    import chemkin_io

    from mechanalyzer.parser import spc as ma_species_io_

    rate0 = data.rate.SimpleRate()
    headers0 = ["name", "inchi", "smiles", "charge", "mult"]

    rxn_df = schema.validate_reactions(rxn_df)
    spc_df = schema.validate_species(spc_df)

    # 1. Form the reactions dictionary
    def reaction_dict_item_(row):
        eq = row[Reactions.eq]
        rate = row[Reactions.rate] if Reactions.rate in row else rate0
        rxn = data.reac.from_equation(eq, rate)
        key = data.reac.chemkin_reagents(rxn, tuple_coll=True)
        val = data.rate.to_old_object(data.reac.rate(rxn))
        return (key, val)

    rxn_dct = dict(rxn_df.progress_apply(reaction_dict_item_, axis=1).to_list())

    # 3. Drop unused species, if requested
    if drop:
        spc_set = {s for r, p, _ in rxn_dct for s in r + p}
        spc_df = spc_df[spc_df[Species.name].isin(spc_set)]
        assert (
            set(spc_df[Species.name]) == spc_set
        ), f'Missing species: {set(spc_df["name"]) - spc_set}'

    # 2. Re-format the species dataframe
    spc_df = schema.validate_species(spc_df, smi=True)
    spc_df[Species.chi] = spc_df[Species.chi].progress_apply(automol.chi.inchi_to_amchi)
    spc_df = spc_df.rename(columns={Species.chi: "inchi", Species.smi: "smiles"})

    # 4. build the mechanism string
    spc_df["fml"] = spc_df["inchi"].progress_apply(automol.amchi.formula)
    spc_df = spc_df[~spc_df["name"].duplicated(keep="first")]
    spc_df = spc_df.set_index("name")
    spc_dct = spc_df.to_dict("index")
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        mech_spc_dct=spc_dct, rxn_param_dct=rxn_dct
    )

    # 5. build the species string
    spc_df = spc_df.drop(columns=spc_df.columns.difference(headers0))
    spc_str = ma_species_io_.csv_string(spc_df.to_dict("index"), spc_df.columns)
    return mech_str, spc_str


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
