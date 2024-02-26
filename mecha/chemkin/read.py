""" CHEMKIN parsers

TODO: Format the reactions and species output into pandas tables
"""

import re
from typing import Dict, List, Optional, Tuple

import pandas
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from mecha import data, schema

# generic
COMMENT_REGEX = re.compile(r"!.*$", flags=re.M)
COMMENT = pp.Suppress(pp.Literal("!")) + ... + pp.Suppress(pp.LineEnd())
COMMENTS = pp.ZeroOrMore(COMMENT)

# units
E_UNIT = pp.Opt(
    pp.CaselessKeyword("CAL/MOLE")
    ^ pp.CaselessKeyword("KCAL/MOLE")
    ^ pp.CaselessKeyword("JOULES/MOLE")
    ^ pp.CaselessKeyword("KJOULES/MOLE")
    ^ pp.CaselessKeyword("KELVINS")
)
A_UNIT = pp.Opt(pp.CaselessKeyword("MOLES") ^ pp.CaselessKeyword("MOLECULES"))

# species
SPECIE = pp.Combine(
    pp.WordStart(pp.alphas) + pp.Word(pp.printables, exclude_chars="+=<>!")
)

# reactions
ARROW = pp.Literal("=") ^ pp.Literal("=>") ^ pp.Literal("<=>")
FALLOFF = pp.Combine(
    pp.Literal("(") + pp.Literal("+") + pp.Literal("M") + pp.Literal(")"),
    adjacent=False,
)
DUP = pp.Opt(pp.CaselessKeyword("DUP") ^ pp.CaselessKeyword("DUPLICATE"))


# reactions
def reactions_block(mech_str: str, comments: bool = True) -> str:
    """Get the reactions block, starting with 'REACTIONS' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "REACTIONS", comments=comments)


def reaction_units(mech_str: str, default: bool = True) -> Tuple[str, str]:
    """Get the E and A units for reaction rate constants

    :param mech_str: A CHEMKIN mechanism string
    :param default: Return default values, if missing?
    :return: The units for E and A, respectively
    """
    e_default = "CAL/MOL" if default else None
    a_default = "MOLES" if default else None

    rxn_block_str = reactions_block(mech_str, comments=False)
    parser = E_UNIT("e_unit") + A_UNIT("a_unit")
    unit_dct = parser.parseString(rxn_block_str).as_dict()
    e_unit = unit_dct["e_unit"].upper() if "e_unit" in unit_dct else e_default
    a_unit = unit_dct["a_unit"].upper() if "a_unit" in unit_dct else a_default
    return e_unit, a_unit


def reaction_kinetics(mech_str: str) -> List[str]:
    """Get the reaction entries (reactions with kinetics and comments, if any)

    :param mech_str: A CHEMKIN mechanism string
    :return: The reaction entries
    """
    # Build the parser
    r_expr = pp.Group(
        pp.delimitedList(SPECIE, delim="+")("species") + pp.Opt(FALLOFF)("falloff")
    )
    eq_expr = r_expr("reactants") + ARROW("arrow") + r_expr("products")
    rxn_expr = (
        eq_expr
        + number_list_expr(3)("arrh")
        + pp.Opt(rate_params_expr("LOW", 3))("arrh0")
        + pp.Opt(rate_params_expr("TROE", 3, 4))("troe")
        + pp.Opt(pp.OneOrMore(pp.Group(rate_params_expr("PLOG", 4))))("plog")
        + DUP("dup")
    )
    parser = pp.Suppress(...) + pp.OneOrMore(pp.Group(rxn_expr))

    # Do the parsing
    rxn_block_str = reactions_block(mech_str, comments=False)
    rcts_lst = []
    prds_lst = []
    rate_lst = []
    for res in parser.parseString(rxn_block_str):
        dct = res.asDict()
        rxn = data.reac.from_chemkin(
            rcts=list(dct["reactants"]["species"]),
            prds=list(dct["products"]["species"]),
            arrow=dct["arrow"],
            plus_m=dct.get("falloff", ""),
            arrh=dct.get("arrh", None),
            arrh0=dct.get("arrh0", None),
            troe=dct.get("arrh0", None),
        )

        rcts_lst.append(data.reac.reactants(rxn))
        prds_lst.append(data.reac.products(rxn))
        rate_lst.append(data.reac.rate(rxn))

    rxn_df = pandas.DataFrame(
        {
            schema.Reactions.reactants: rcts_lst,
            schema.Reactions.products: prds_lst,
            schema.Reactions.rate: rate_lst,
        }
    )

    return schema.validate_reactions(rxn_df)


# species
def species_block(mech_str: str, comments: bool = True) -> str:
    """Get the species block, starting with 'SPECIES' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "SPECIES", comments=comments)


def species(mech_str: str) -> List[str]:
    """Get the list of species

    :param mech_str: A CHEMKIN mechanism string
    :return: The species
    """
    parser = pp.OneOrMore(SPECIE)
    spc_block_str = species_block(mech_str, comments=False)
    return parser.parseString(spc_block_str).asList()


def species_with_comments(mech_str: str) -> Dict[str, List[str]]:
    """Get the list of species, along with their comments

    :param mech_str: A CHEMKIN mechanism string
    :return: A dictionary mapping species onto their comments
    """
    parser = pp.Suppress(...) + pp.OneOrMore(pp.Group(SPECIE + pp.Group(COMMENTS)))
    spc_block_str = species_block(mech_str, comments=True)
    return dict(parser.parseString(spc_block_str).asList())


# therm
def therm_block(mech_str: str) -> str:
    """Get the therm block, starting with 'REACTIONS' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "THERM")


# generic
def block(mech_str: str, key: str, comments: bool = False) -> str:
    """Get a keyword block, starting with a key and ending in 'END'

    :param mech_str: The mechanism string
    :param key: The key that the block starts with
    :param comments: Include comments?
    :return: The block
    """
    block_par = pp.Suppress(...) + pp.QuotedString(
        key, end_quote_char="END", multiline=True
    )
    (block_str,) = block_par.parseString(mech_str).asList()
    # Remove comments, if requested
    if not comments:
        block_str = without_comments(block_str)
    return block_str


def without_comments(mech_str: str) -> str:
    """Get a CHEMKIN string or substring with comments removed

    :param mech_str: The mechanism string, or any substring of it
    :return: The string, without comments
    """
    return re.sub(COMMENT_REGEX, "", mech_str)


def all_comments(mech_str: str) -> List[str]:
    """Get all comments from a CHEMKIN string or substring

    :param mech_str: The mechanism string, or any substring of it
    :return: The comments
    """
    return re.findall(COMMENT_REGEX, mech_str)


# helpers
def number_list_expr(
    nmin: int, nmax: Optional[int] = None, delim: str = ""
) -> pp.core.ParseExpression:
    """A parse expression for a list of numbers

    :param nmin: The minimum list length
    :param nmax: The maximum list length (defaults to `nmin` if `None`)
    :param delim: The delimiter between numbers, defaults to ""
    :return: The parse expression
    """
    nmax = nmin if nmax is None else nmax
    return pp.delimitedList(ppc.number, delim=delim, min=nmin, max=nmax)


def rate_params_expr(
    key: str, nmin: int, nmax: Optional[int] = None
) -> pp.core.ParseExpression:
    """A parse expression for rate parameters after a CHEMKIN reaction

    :param key: The keyword for these rate parameters
    :param nmin: The minimum parameter list length
    :param nmax: The maximum parameter list length (defaults to `nmin` if `None`)
    :return: The parse expression
    """
    keyword = pp.Suppress(pp.CaselessLiteral(key))
    slash = pp.Suppress(pp.Literal("/"))
    params = number_list_expr(nmin, nmax=nmax, delim="")
    return keyword + slash + params + slash
