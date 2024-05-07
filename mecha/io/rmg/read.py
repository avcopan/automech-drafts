""" RMG parsers
"""

import os
from typing import Any, Dict, Optional

import automol
import pandas
import pyparsing as pp
from automol.graph import RMG_ADJACENCY_LIST
from pyparsing import pyparsing_common as ppc
from tqdm.auto import tqdm

from mecha import schema
from mecha.data.reac import SPECIES_NAME

MULTIPLICITY = pp.CaselessLiteral("multiplicity") + ppc.integer("mult")
SPECIES_ENTRY = (
    SPECIES_NAME("species") + pp.Opt(MULTIPLICITY) + RMG_ADJACENCY_LIST("adj_list")
)
SPECIES_DICT = pp.OneOrMore(pp.Group(SPECIES_ENTRY))("dict")


def species_dictionary(inp: str, out: Optional[str] = None) -> Dict[str, Any]:
    """Parse a species dictionary string

    :param inp: An RMG species dictionary, as a file path or string
    :param out: Optionally, write the output to this file path
    :return: A dictionary mapping CHEMKIN names onto automol graphs
    """
    inp = open(inp).read() if os.path.exists(inp) else inp

    spc_par_rets = SPECIES_DICT.parseString(inp).asDict()["dict"]

    names = []
    mults = []
    smis = []
    chis = []
    for spc_par_ret in tqdm(spc_par_rets):
        adj_par_ret = spc_par_ret["adj_list"]
        gra = automol.graph.from_parsed_rmg_adjacency_list(adj_par_ret)

        names.append(spc_par_ret["species"])
        mults.append(spc_par_ret.get("mult", 1))
        chis.append(automol.graph.inchi(gra))
        smis.append(automol.graph.smiles(gra))

    spc_df = pandas.DataFrame(
        {
            schema.Species.name: names,
            schema.Species.mult: mults,
            schema.Species.chi: chis,
            schema.Species.smi: smis,
        }
    )

    spc_df = schema.validate_species(spc_df)
    if out is not None:
        spc_df.to_csv(out)

    return spc_df
