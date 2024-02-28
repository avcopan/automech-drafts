"""Reaction dataclasses
"""

import dataclasses
from typing import Dict, Optional, Tuple

import pyparsing as pp

from mecha.data import rate as rate_


@dataclasses.dataclass
class Reaction:
    """A reaction

    Stores a tuple of rates to handle duplicate reactions

    :param reactants: Names for the reactants
    :param products: Names for the products
    :param rates: The reaction rates
    """

    reactants: Tuple[str, ...]
    products: Tuple[str, ...]
    rates: Tuple[rate_.Rate, ...] = ()


# constructors
def from_chemkin(
    rcts: Tuple[str, ...],
    prds: Tuple[str, ...],
    arrow: str = "=",
    plus_m: str = "",
    arrh: Optional[rate_.Params3] = None,
    arrh0: Optional[rate_.Params4] = None,
    troe: Optional[rate_.Params3or4] = None,
    plog: Optional[Tuple[rate_.Params4, ...]] = None,
) -> Reaction:
    """Build a single-rate Reaction object from data

    :param rcts: The CHEMKIN reactants
    :param prds: The CHEMKIN products
    :param arrow: The CHEMKIN arrow, indicating whether or not the reaction is reversible
    :param plus_m: The CHEMKIN M collider, 'M' or '(+M)', indicating the type of
        pressure dependence for simple reactions
    :param arrh: The high-pressure Arrhenius parameters, defaults to None
    :param arrh0: The low-pressure Arrhenius parameters, defaults to None
    :param troe: The Troe parameters, defaults to None
    :param plog: The Plog parameters, defaults to None
    :return: The reaction object
    """
    colliders = ("M", "(M)")

    rcts = [s.replace(" ", "").replace("+", "") for s in rcts]
    prds = [s.replace(" ", "").replace("+", "") for s in prds]

    plus_mr = next((s for s in rcts if s in colliders), plus_m)
    plus_mp = next((s for s in prds if s in colliders), plus_m)
    assert (
        plus_mr == plus_mp
    ), f"Inconsistent colliders, {plus_mr} != {plus_mp}, for reaction {rcts} = {prds}"

    rcts = [s for s in rcts if s not in colliders]
    prds = [s for s in prds if s not in colliders]

    rate = rate_.from_chemkin(
        arrow=arrow, plus_m=plus_mr, arrh=arrh, arrh0=arrh0, troe=troe, plog=plog
    )
    return Reaction(reactants=tuple(rcts), products=tuple(prds), rates=(rate,))


# getters
def reactants(rxn: Reaction) -> Tuple[str, ...]:
    """The list of reactants

    :param rxn: A reaction object
    :return: The CHEMKIN names of the reactants
    """
    return rxn.reactants


def products(rxn: Reaction) -> Tuple[str, ...]:
    """The list of products

    :param rxn: A reaction object
    :return: The CHEMKIN names of the products
    """
    return rxn.products


def rates(rxn: Reaction) -> Tuple[rate_.Rate, ...]:
    """The rate constants

    :param rxn: A reaction object
    :return: The rate objects
    """
    return rxn.rates


def rate(rxn: Reaction) -> Tuple[rate_.Rate]:
    """The rate constant or, if multiple, the (arbitrary) first one in the list

    :param rxn: A reaction object
    :return: The rate object
    """
    return rxn.rates[0]


# properties
def equation(rxn: Reaction) -> str:
    """Get the CHEMKIN equation of a reaction

    :param rxn: A reaction object
    :return: The reaction CHEMKIN equation
    """
    return form_equation(reactants(rxn), products(rxn))


# helpers
SPECIE = pp.Combine(
    pp.WordStart(pp.alphas) + pp.Word(pp.printables, exclude_chars="+=<>!")
)
ARROW = pp.Literal("=") ^ pp.Literal("=>") ^ pp.Literal("<=>")
FALLOFF = pp.Combine(
    pp.Literal("(") + pp.Literal("+") + pp.Literal("M") + pp.Literal(")"),
    adjacent=False,
)


def parse_equation(
    eq: str, trans_dct: Optional[Dict[str, str]] = None
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """Parse the CHEMKIN equation of a reaction into reactants and products

    :param eq: The reaction CHEMKIN equation
    :return: The reactants and products
    """

    def trans_(name):
        return name if trans_dct is None else trans_dct.get(name)

    r_expr = pp.Group(
        pp.delimitedList(SPECIE, delim="+")("species") + pp.Opt(FALLOFF)("falloff")
    )
    parser = r_expr("reactants") + ARROW("arrow") + r_expr("products")
    dct = parser.parseString(eq).asDict()
    rcts = tuple(map(trans_, dct["reactants"]["species"]))
    prds = tuple(map(trans_, dct["products"]["species"]))
    return (rcts, prds)


def equation_reactants(eq: str) -> Tuple[str, ...]:
    """Get the reactants of a CHEMKIN equation

    :param eq: The reaction CHEMKIN equation
    :return: The reactants
    """
    rcts, _ = parse_equation(eq)
    return rcts


def equation_products(eq: str) -> Tuple[str, ...]:
    """Get the products of a CHEMKIN equation

    :param eq: The reaction CHEMKIN equation
    :return: The products
    """
    _, prds = parse_equation(eq)
    return prds


def equation_reagents(eq: str, prod: bool = False) -> Tuple[str, ...]:
    """Get the reagents (reactants or products) of a CHEMKIN equation

    :param eq: The reaction CHEMKIN equation
    :param prod: Get the products, instead of the reactants?, defaults to False
    :return: The reagents
    """
    return equation_products(eq) if prod else equation_reactants(eq)


def form_equation(rcts: Tuple[str, ...], prds: Tuple[str, ...]) -> str:
    """Form the CHEMKIN equation of a reaction from reactants and products

    :param rcts: The reactant names
    :param prds: The product names
    :return: The reaction CHEMKIN equation
    """
    assert all(isinstance(n, str) for n in rcts), f"Invalid reactants: {rcts}"
    assert all(isinstance(n, str) for n in prds), f"Invalid products: {prds}"

    rcts_str = " + ".join(rcts)
    prds_str = " + ".join(prds)
    return " = ".join([rcts_str, prds_str])
