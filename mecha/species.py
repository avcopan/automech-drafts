"""Functions concerning only the species in the mechanism
"""

import automol
import pandas
from tqdm.auto import tqdm

from mecha import chemkin, schema
from mecha.schema import Species

tqdm.pandas()


def expand_stereo(spc_df: pandas.DataFrame, enant: bool = True) -> pandas.DataFrame:
    """Expand stereochemistry for a list of species

    :param spc_df: The species dataframe
    :param enant: Distinguish between enantiomers?, defaults to True
    :return: The stereo-expanded species dataframe
    """

    def expand_amchi_(chi):
        return automol.amchi.expand_stereo(chi, enant=enant)

    def name_(row):
        name = row[Species.orig_name]
        chi = row[Species.chi]
        return chemkin.name.with_stereo_suffix(name, chi, racem=not enant)

    spc_df = schema.validate_species(spc_df)
    spc_df = rename_with_original_columns(spc_df)
    spc_df[Species.chi] = spc_df[Species.orig_chi].progress_apply(expand_amchi_)
    spc_df = spc_df.explode(Species.chi)
    spc_df[Species.name] = spc_df.apply(name_, axis=1)
    return schema.validate_species(spc_df)


# Helpers
def rename_with_original_columns(spc_df: pandas.DataFrame) -> pandas.DataFrame:
    """Rename the species dataframe with original column names
    (orig_name, orig_chi, orig_smi)

    Used when expanding stereo

    :param spc_df: The species dataframe
    :return: The species dataframe
    """
    return spc_df.rename(columns=dict(zip(schema.S_CURR_COLS, schema.S_ORIG_COLS)))
