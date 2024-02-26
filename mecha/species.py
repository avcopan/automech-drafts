"""Functions concerning only the species in the mechanism
"""

import automol
import pandas

from mecha import chemkin, schema


def expand_stereo(spc_df: pandas.DataFrame, enant: bool = True) -> pandas.DataFrame:
    """Expand stereochemistry for a list of species

    :param spc_df: The species dataframe
    :param enant: Distinguish between enantiomers?, defaults to True
    :return: The stereo-expanded species dataframe
    """

    def expand_amchi(chi):
        return automol.amchi.expand_stereo(chi, enant=enant)

    def update_name(row):
        orig_name = row[schema.Species.orig_name]
        chi = row[schema.Species.chi]
        return chemkin.name.with_stereo_suffix(orig_name, chi, racem=not enant)

    spc_df = schema.validate_species(spc_df)
    spc_df = spc_df.rename(
        columns={
            schema.Species.name: schema.Species.orig_name,
            schema.Species.chi: schema.Species.orig_chi,
            schema.Species.smi: schema.Species.orig_smi,
        }
    )
    spc_df[schema.Species.chi] = spc_df[schema.Species.orig_chi].apply(expand_amchi)
    spc_df = spc_df.explode(schema.Species.chi)
    spc_df[schema.Species.name] = spc_df.apply(update_name, axis=1)
    return schema.validate_species(spc_df)
