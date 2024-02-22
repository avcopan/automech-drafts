""" RMG parsers
"""

from typing import Any, Dict

import automol
import pandas
import pyparsing as pp
from automol.graph import RMG_ADJACENCY_LIST
from pyparsing import pyparsing_common as ppc

from autoio import schema
from autoio.chemkin import SPECIE

MULTIPLICITY = pp.CaselessLiteral("multiplicity") + ppc.integer("mult")
SPECIES_ENTRY = (
    SPECIE("species") + pp.Opt(MULTIPLICITY) + RMG_ADJACENCY_LIST("adj_list")
)
SPECIES_DICT = pp.OneOrMore(pp.Group(SPECIES_ENTRY))("dict")


def species_dictionary(spc_dict_str) -> Dict[str, Any]:
    """Parse a species dictionary string

    :param spc_dict_str: A species dictionary string, consisting of alternating CHEMKIN
        names and adjacency lists
    :return: A dictionary mapping CHEMKIN names onto automol graphs
    """
    spc_par_rets = SPECIES_DICT.parseString(spc_dict_str).asDict()["dict"]
    names = []
    mults = []
    smis = []
    chis = []
    for spc_par_ret in spc_par_rets:
        adj_par_ret = spc_par_ret["adj_list"]
        gra = automol.graph.from_parsed_rmg_adjacency_list(adj_par_ret)

        names.append(spc_par_ret["species"])
        mults.append(spc_par_ret.get("mult", 1))
        smis.append(automol.graph.smiles(gra))
        chis.append(automol.graph.inchi(gra))

    spc_df = pandas.DataFrame(
        {
            schema.Species.name: names,
            schema.Species.mult: mults,
            schema.Species.smi: smis,
            schema.Species.chi: chis,
        }
    )
    return schema.validate(schema.Species, spc_df)
