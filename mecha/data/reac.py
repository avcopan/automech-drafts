"""Reaction dataclasses
"""

import dataclasses
from typing import Optional, Tuple

from mecha.data import rate as rate_


@dataclasses.dataclass
class Reaction:
    """A reaction

    Stores a tuple of rates to handle duplicate reactions

    :param reactants: Names for the reactants
    :param products: Names for the products
    :param rates: The reaction rates
    """

    reactants: Tuple[str]
    products: Tuple[str]
    rates: Tuple[rate_.Rate] = ()


# constructors
def from_chemkin(
    rcts: Tuple[str],
    prds: Tuple[str],
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
def reactants(rxn: Reaction) -> Tuple[str]:
    """The list of reactants

    :param rxn: A reaction object
    :return: The CHEMKIN names of the reactants
    """
    return rxn.reactants


def products(rxn: Reaction) -> Tuple[str]:
    """The list of products

    :param rxn: A reaction object
    :return: The CHEMKIN names of the products
    """
    return rxn.products


def rates(rxn: Reaction) -> Tuple[rate_.Rate]:
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
    """Get the chemical equation for a reaction, as a CHEMKIN string

    :param rxn: A reaction object
    :return: The chemical equation
    """
