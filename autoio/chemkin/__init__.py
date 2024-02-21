""" CHEMKIN parsers
"""

import re
from typing import Dict, List, Optional, Tuple

import pyparsing as pp
from pyparsing import pyparsing_common as ppc


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


# generic
COMMENT = pp.Suppress(pp.Literal("!")) + ... + pp.Suppress(pp.LineEnd())
COMMENTS = pp.delimitedList(COMMENT, delim=pp.LineEnd())

# units
E_UNIT = (
    pp.CaselessKeyword("CAL/MOLE")
    ^ pp.CaselessKeyword("KCAL/MOLE")
    ^ pp.CaselessKeyword("JOULES/MOLE")
    ^ pp.CaselessKeyword("KJOULES/MOLE")
    ^ pp.CaselessKeyword("KELVINS")
)
A_UNIT = pp.CaselessKeyword("MOLES") ^ pp.CaselessKeyword("MOLECULES")
REACTION_UNITS = pp.Opt(E_UNIT)("e_unit") + pp.Opt(A_UNIT)("a_unit")

# species
SPECIE = pp.Combine(
    pp.WordStart(pp.alphas) + pp.Word(pp.printables, exclude_chars="+=<>!")
)
SPECIES = pp.OneOrMore(SPECIE("species"))
SPECIES_WITH_COMMENTS = pp.OneOrMore(
    pp.Group(SPECIE("species") + pp.Opt(COMMENTS)("comments"))
)

# reactions
ARROW = pp.Literal("=") ^ pp.Literal("<=") ^ pp.Literal("=>") ^ pp.Literal("<=>")
FALLOFF = pp.Combine(
    pp.Literal("(") + pp.Literal("+") + pp.Literal("M") + pp.Literal(")"),
    adjacent=False,
)
REAGENTS = pp.Group(
    pp.delimitedList(SPECIE, delim="+")("species") + pp.Opt(FALLOFF)("falloff")
)
REACTION = REAGENTS("reactants") + ARROW("arrow") + REAGENTS("products")
RATE = (
    number_list_expr(3)("arrh")
    + pp.Opt(rate_params_expr("LOW", 3))("lind")
    + pp.Opt(rate_params_expr("TROE", 3, 4))("troe")
    + pp.ZeroOrMore(rate_params_expr("PLOG", 4))("plog")
)
REACTION_ENTRY = pp.Group(REACTION)("reaction") + pp.Group(RATE)("rate")


# reactions
def reactions_block(mech_str: str) -> str:
    """Get the reactions block, starting with 'REACTIONS' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "REACTIONS")


def reaction_units(mech_str: str, default: bool = True) -> Tuple[str, str]:
    """Get the E and A units for reaction rate constants

    :param mech_str: A CHEMKIN mechanism string
    :param default: Return default values, if missing?
    :return: The units for E and A, respectively
    """
    e_default = "CAL/MOL" if default else None
    a_default = "MOLES" if default else None

    reac_str = reactions_block(mech_str)
    unit_expr = REACTION_UNITS
    unit_dct = unit_expr.parseString(reac_str).as_dict()
    e_unit = unit_dct["e_unit"].upper() if "e_unit" in unit_dct else e_default
    a_unit = unit_dct["a_unit"].upper() if "a_unit" in unit_dct else a_default
    return e_unit, a_unit


# species
def species_block(mech_str: str) -> str:
    """Get the species block, starting with 'SPECIES' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "SPECIES")


def species(mech_str: str) -> List[str]:
    """Get the list of species

    :param mech_str: A CHEMKIN mechanism string
    :return: The species
    """
    spc_block_str = without_comments(species_block(mech_str))
    return SPECIES.parseString(spc_block_str).asList()


def species_with_comments(mech_str: str) -> Dict[str, List[str]]:
    """Get the list of species, along with their comments

    :param mech_str: A CHEMKIN mechanism string
    :return: A dictionary mapping species onto their comments
    """
    spc_block_str = species_block(mech_str)
    dct = SPECIES_WITH_COMMENTS("_").parseString(spc_block_str).asDict()
    print(dct)


# therm
def therm_block(mech_str: str) -> str:
    """Get the therm block, starting with 'REACTIONS' and ending in 'END'

    :param mech_str: A CHEMKIN mechanism string
    :return: The block
    """
    return block(mech_str, "THERM")


# generic
def block(mech_str: str, key: str) -> str:
    """Get a keyword block, starting with a key and ending in 'END'

    :param mech_str: The mechanism string
    :param key: The key that the block starts with
    :return: The block
    """
    block_par = pp.Suppress(...) + pp.QuotedString(
        key, end_quote_char="END", multiline=True
    )
    (block_str,) = block_par.parseString(mech_str).asList()
    return block_str


COMMENT_REGEX = re.compile(r"!.*$", flags=re.M)


def without_comments(mech_str: str) -> str:
    """Get a CHEMKIN string with comments removed

    :param mech_str: The mechanism string, or any substring of it
    :return: The string, without comments
    """
    return re.sub(COMMENT_REGEX, "", mech_str)
