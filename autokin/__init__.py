"""Handling of data related to thermodynamics and kinetics
"""

import abc
import dataclasses
from typing import Optional, Tuple


@dataclasses.dataclass
class SpeciesThermo:
    """Species thermodynamic properties parametrization (see cantera.SpeciesThermo)

    Types:
        NASA7   - coeffs: a0, a1, a2, a3, a4, a5, a6
        NASA9   - coeffs: a0, a1, a2, a3, a4, a5, a6, a7, a8
        Shomate - coeffs: A, B, C, D, E, F, G

    :param T0: The minimum temperature [K] at which the parametrization is valid
    :param T_: The maximum temperature [K] at which the parametrization is valid
    :param coeffs: A list of coefficients for the parametrization
    :param type_: The type of parametrization: "NASA7", "NASA9", "Shomate"
    :param P_ref: The reference pressure [Pa] for the parametrization
    """

    T0: float
    T_: float
    coeffs: Tuple[float, ...]
    type_: str = "NASA7"
    P_ref: float = 101325  # 1 atm


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
class FalloffFunction:
    """A Lindemann, Troe, or SRI function (see cantera.Falloff)

    Types:
        Lindemann   - coeffs: (None)
        Troe        - coeffs: a, T, T***, T*, (T**)
        SRI         - coeffs: a, b, c, (d), (e)

    :param coeffs: A list of coefficients for the parametrization
    :param type_: The type of parametrization: "Lindemann", "Troe", "SRI"
    """

    coeffs: Optional[Tuple[float, ...]] = None
    type_: str = "Lindemann"


class ReactionRate(abc.ABC):
    """Base class for reaction rates"""

    @property
    @abc.abstractproperty
    def type_(self):
        pass


@dataclasses.dataclass
class ReactionRateSimple(ReactionRate):
    """Simple reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    Types:
        Elementary  - k: The rate coefficient
        Falloff     - k: The high-pressure rate coefficient (M-independent)
                    - k0: The low-pressure rate coefficient
                    - f: The blending function, F(T, P_r)
        Activated   - k: The high-pressure rate coefficient
                    - k0: The low-pressure rate coefficient (M-independent)
                    - f: The blending function, F(T, P_r)

    :param k: Rate coefficient (or high-pressure limit) for the reaction
    :param k0: Rate coefficient for the low-pressure limit
    :param f: Falloff function for blending the high- and low-pressure rate coefficients
    :param type_: The type of parametrization: "Elementary", "Falloff", "Activated"
    """

    k: ArrheniusFunction
    k0: Optional[ArrheniusFunction] = None
    f: Optional[FalloffFunction] = None
    type_: str = "Elementary"

    def __post_init__(self):
        assert self.type_ in ("Elementary", "Falloff", "Activated")

        if self.type_ != "Elementary":
            self.f = FalloffFunction() if self.f is None else self.f


@dataclasses.dataclass
class ReactionRatePlog(ReactionRate):
    """P-Log reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    :param ks: Rate coefficients at specific pressures, k_P_, k_P2, ...
    :param Ps: An array of pressures, P_, P2, ... [Pa]
    """

    ks: Tuple[ArrheniusFunction, ...]
    Ps: Tuple[float, ...]
    type_: str = "Plog"

    def __post_init__(self):
        assert self.type_ == "Plog"


@dataclasses.dataclass
class ReactionRateCheb(ReactionRate):
    """Chebyshev reaction rate, k(T,P) parametrization (see cantera.ReactionRate)

    :param T0: The minimum temperature [K] for the Chebyshev fit
    :param T_: The minimum temperature [K] for the Chebyshev fit
    :param P0: The minimum pressure [K] for the Chebyshev fit
    :param P_: The minimum pressure [K] for the Chebyshev fit
    :param ks: Arrhenius rate expressions in the parametrization, if any
    :param coeffs: Pressure-dependence coefficients in the parametrization, if any
    :param Ps: Pressure range P0 -> P (Chebyshev) or array, P_, P2, ... (Plog) [Pa]
    :param Ts: Temperature range T0 -> T (Chebyshev) [K]
    """

    T0: float
    T_: float
    P0: float
    P_: float
    coeffs: Tuple[Tuple[float, ...], ...]
    type_: str = "Chebyshev"

    def __post_init__(self):
        assert self.type_ == "Chebyshev"
