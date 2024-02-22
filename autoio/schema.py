from typing import Optional

import pandas
import pandera as pa
from pandera.typing import Series


class Species(pa.DataFrameModel):
    name: Series[str] = pa.Field(coerce=True)
    mult: Series[int] = pa.Field(coerce=True)
    smi: Optional[Series[str]]
    chi: Optional[Series[str]]
    charge: Series[int] = pa.Field(coerce=True, default=0)


def validate(model: pa.DataFrameModel, df: pandas.DataFrame) -> pandas.DataFrame:
    """Validate a pandas dataframe based on a model

    :param model: The model
    :param df: The DataFrame
    :return: The validated DataFrame
    """
    schema = model.to_schema()
    schema.add_missing_columns = True
    return schema.validate(df)
