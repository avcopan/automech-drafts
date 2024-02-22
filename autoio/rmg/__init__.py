""" RMG parsers
"""

from typing import Any, Dict

import automol
import pandas
import pyparsing as pp
from automol.graph import RMG_ADJACENCY_LIST
from pyparsing import pyparsing_common as ppc

from autoio import schema
from autoio.chemkin._read import SPECIE

MULTIPLICITY = pp.CaselessLiteral("multiplicity") + ppc.integer("mult")
SPECIES_ENTRY = (
    SPECIE("species") + pp.Opt(MULTIPLICITY) + RMG_ADJACENCY_LIST("adj_list")
)
SPECIES_DICT = pp.OneOrMore(pp.Group(SPECIES_ENTRY))("dict")


def species_dictionary(rmg_spc_str) -> Dict[str, Any]:
    """Parse a species dictionary string

    :param rmg_spc_str: An RMG species dictionary string
    :return: A dictionary mapping CHEMKIN names onto automol graphs
    """
    spc_par_rets = SPECIES_DICT.parseString(rmg_spc_str).asDict()["dict"]
    names = []
    mults = []
    smis = []
    chis = []
    for spc_par_ret in spc_par_rets:
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
    return schema.validate_species(spc_df, smi=True)
