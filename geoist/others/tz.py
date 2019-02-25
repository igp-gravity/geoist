from typing import Union
import numpy as np
import datetime
from pytz import UTC
from dateutil.parser import parse


def forceutc(t: Union[str, datetime.datetime, datetime.date, np.datetime64]) -> Union[datetime.datetime, datetime.date]:
    """
    Add UTC to datetime-naive and convert to UTC for datetime aware

    input: python datetime (naive, utc, non-utc) or Numpy datetime64  #FIXME add Pandas and AstroPy time classes
    output: utc datetime
    """
    # need to passthrough None for simpler external logic.
# %% polymorph to datetime
    if isinstance(t, str):
        t = parse(t)
    elif isinstance(t, np.datetime64):
        t = t.astype(datetime.datetime)
    elif isinstance(t, datetime.datetime):
        pass
    elif isinstance(t, datetime.date):
        return t
    elif isinstance(t, (np.ndarray, list, tuple)):
        return np.asarray([forceutc(T) for T in t])
    else:
        raise TypeError('datetime only input')
# %% enforce UTC on datetime
    if t.tzinfo is None:  # datetime-naive
        t = t.replace(tzinfo=UTC)
    else:  # datetime-aware
        t = t.astimezone(UTC)  # changes timezone, preserving absolute time. E.g. noon EST = 5PM UTC

    return t
