"""DataFrame utilities
"""

from typing import Dict, Tuple, Union

import pandas

Key = str
Keys = Tuple[str, ...]
Key_ = Union[Key, Keys]
Value = object
Values = Tuple[object, ...]
Value_ = Union[Value, Values]


def lookup_dict(df: pandas.DataFrame, in_: Key_, out_: Key_) -> Dict[Value_, Value_]:
    """Form a lookup dictionary mapping one column onto another in a dataframe

    Allows mappings between sets of columns, using a tuple of column keys

    :param df: The dataframe
    :param in_: The input key or keys
    :param out_: The output key or keys
    :return: The dictionary mapping input values to output values
    """

    def check_(key_):
        return key_ in df if isinstance(key_, str) else all(k in df for k in key_)

    def values_(key_):
        return df[key_] if isinstance(key_, str) else zip(*(df[k] for k in key_))

    assert check_(in_), f"{in_} not in {df}"
    assert check_(out_), f"{out_} not in {df}"

    return dict(zip(values_(in_), values_(out_)))
