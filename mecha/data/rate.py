"""Kinetics dataclasses

Eventually, we will want to add a "Raw" rate type with k(T,P) values stored in an xarray
"""

import abc
import dataclasses
import enum
from typing import Optional, Tuple, Union

Params3 = Tuple[float, float, float]
Params4 = Tuple[float, float, float, float]
Params3or4 = Union[Params3, Params4]


class RateType(str, enum.Enum):
    """The type of reaction rate (type of pressure dependence)"""

    CONSTANT = "Constant"
    FALLOFF = "Falloff"
    ACTIVATED = "Activated"
    PLOG = "Plog"
    CHEB = "Chebyshev"


class BlendType(str, enum.Enum):
    """The type of blending function for high and low-pressure rates"""

    LIND = "Lind"
    TROE = "Troe"


@dataclasses.dataclass
class ArrheniusFunction:
    """An Arrhenius or Landau-Teller function (see cantera.Arrhenius)

    :param A: The pre-exponential factor [(m^3/kmol)**o/s]
    :param b: The temperature exponent
    :param E: The activation energy E [J/kmol]
    :param B: The Landau-Teller B-factor
    :param C: The Landau-Teller C-factor
    """

    A: float = 1.0
    b: float = 0.0
    E: float = 0.0
    B: float = 0.0
    C: float = 0.0


@dataclasses.dataclass
class BlendingFunction:
    """A blending function for high and low-pressure rates (see cantera.Falloff)

    Types:
        Lind   - coeffs: (None)
        Troe   - coeffs: a, T, T***, T*, (T**)

    :param coeffs: A list of coefficients for the parametrization
    :param type_: The type of parametrization: "Lind", "Troe"
    """

    coeffs: Optional[Tuple[float, ...]] = None
    type_: BlendType = BlendType.LIND

    def __post_init__(self):
        self.type_ = BlendType(self.type_)


class Rate(abc.ABC):
    """Base class for reaction rates"""

    @property
    @abc.abstractproperty
    def type_(self):
        pass

    @property
    @abc.abstractproperty
    def is_rev(self):
        pass


@dataclasses.dataclass
class SimpleRate(Rate):
    """Simple reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    Types:
        Constant    - k: The rate coefficient
        Falloff     - k: The high-pressure rate coefficient (M-independent)
                    - k0: The low-pressure rate coefficient
                    - f: The blending function, F(T, P_r)
        Activated   - k: The high-pressure rate coefficient
                    - k0: The low-pressure rate coefficient (M-independent)
                    - f: The blending function, F(T, P_r)

    :param k: Rate coefficient (or high-pressure limit) for the reaction
    :param k0: Rate coefficient for the low-pressure limit
    :param f: Falloff function for blending the high- and low-pressure rate coefficients
    :param is_rev: Is this a reversible reaction?
    :param type_: The type of reaction: "Constant", "Falloff", "Activated"
    """

    k: ArrheniusFunction = ArrheniusFunction()
    k0: Optional[ArrheniusFunction] = None
    f: Optional[BlendingFunction] = None
    is_rev: bool = True
    type_: RateType = RateType.CONSTANT

    def __post_init__(self):
        self.type_ = RateType(self.type_)
        assert self.type_ in (RateType.CONSTANT, RateType.FALLOFF, RateType.ACTIVATED)

        if self.type_ != RateType.CONSTANT:
            self.f = BlendingFunction() if self.f is None else self.f


@dataclasses.dataclass
class PlogRate(Rate):
    """P-Log reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    :param ks: Rate coefficients at specific pressures, k_P_, k_P2, ...
    :param Ps: An array of pressures, P_, P2, ... [Pa]
    :param k: Optional high-pressure rate
    :param is_rev: Is this a reversible reaction?
    """

    ks: Tuple[ArrheniusFunction, ...]
    Ps: Tuple[float, ...]
    k: Optional[ArrheniusFunction] = None
    is_rev: bool = True
    type_: RateType = RateType.PLOG

    def __post_init__(self):
        self.type_ = RateType(self.type_)
        assert self.type_ == RateType.PLOG


def from_chemkin(
    arrow: str = "=",
    plus_m: str = "",
    arrh: Optional[Params3] = None,
    arrh0: Optional[Params4] = None,
    troe: Optional[Params3or4] = None,
    plog: Optional[Tuple[Params4, ...]] = None,
) -> Rate:
    """Create a rate object from CHEMKIN data

    Ignores Chebyshev for now...

    :param arrow: The CHEMKIN arrow, indicating whether or not the reaction is reversible
    :param plus_m: The CHEMKIN M collider, 'M' or '(+M)', indicating the type of
        pressure dependence for simple reactions
    :param arrh: The high-pressure Arrhenius parameters, defaults to None
    :param arrh0: The low-pressure Arrhenius parameters, defaults to None
    :param troe: The Troe parameters, defaults to None
    :param plog: The Plog parameters, defaults to None
    :return: The rate object
    """
    # Assess reversibility based on the arrow
    arrow = arrow.strip()
    assert arrow in ("=", "<=>", "=>"), f"Invalid CHEMKIN arrow: {arrow}"
    is_rev = arrow in ("=", "<=>")
    # Determine the high and low-pressure Arrhenius constants, if present
    k = None if arrh is None else ArrheniusFunction(*arrh)
    k0 = None if arrh0 is None else ArrheniusFunction(*arrh0)

    # If this is a Plog rate, return early
    # Chebyshev rates are could be handled similarly, if needed
    if plog is not None:
        ks = [c[1:] for c in plog]
        Ps = ([c[0] for c in plog],)
        return PlogRate(ks=ks, Ps=Ps, k=k, is_rev=is_rev, type_=RateType.PLOG)

    # Otherwise, this is a simple case
    # Determine the blending function if Troe coefficients are present
    f = None if troe is None else BlendingFunction(troe, type_=BlendType.TROE)

    # Determine the pressure dependency type from the M collider
    m2t = {"": RateType.CONSTANT, "M": RateType.ACTIVATED, "(M)": RateType.FALLOFF}
    type_ = m2t[plus_m.replace(" ", "").replace("+", "")]
    assert f is None or type_ in (
        RateType.ACTIVATED,
        RateType.FALLOFF,
    ), "Troe coefficients without +M"

    return SimpleRate(k=k, k0=k0, f=f, is_rev=is_rev, type_=type_)


@dataclasses.dataclass
class ChebRate(Rate):
    """Chebyshev reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    :param T0: The minimum temperature [K] for the Chebyshev fit
    :param T_: The minimum temperature [K] for the Chebyshev fit
    :param P0: The minimum pressure [K] for the Chebyshev fit
    :param P_: The minimum pressure [K] for the Chebyshev fit
    :param coeffs: The Chebyshev expansion coefficients
    :param is_rev: Is this a reversible reaction?
    """

    T0: float
    T_: float
    P0: float
    P_: float
    coeffs: Tuple[Tuple[float, ...], ...]
    is_rev: bool = True
    type_: str = RateType.CHEB

    def __post_init__(self):
        self.type_ = RateType(self.type_)
        assert self.type_ == RateType.CHEB
