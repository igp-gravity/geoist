"""These functions separated due to time-consuming imports"""
from datetime import timedelta
import numpy as np
import xarray
from matplotlib.dates import DateFormatter
from matplotlib.dates import MinuteLocator, SecondLocator


def tickfix(t, fg, ax, tfmt: str='%H:%M:%S'):
    majtick, mintick = timeticks(t[-1] - t[0])
    if majtick:
        ax.xaxis.set_major_locator(majtick)
    if mintick:
        ax.xaxis.set_minor_locator(mintick)
    ax.xaxis.set_major_formatter(DateFormatter(tfmt))
    fg.autofmt_xdate()

    ax.autoscale(True, 'x', tight=True)

    # ax.tick_params(axis='both',which='both')
    # ax.grid(True,which='both')


def timeticks(tdiff):
    """
    NOTE do NOT use "interval" or ticks are misaligned!  use "bysecond" only!
    """
    if isinstance(tdiff, xarray.DataArray):  # len==1
        tdiff = timedelta(seconds=tdiff.values / np.timedelta64(1, 's'))

    assert isinstance(tdiff, timedelta), 'expecting datetime.timedelta'

    if tdiff > timedelta(hours=2):
        return None, None

    elif tdiff > timedelta(minutes=20):
        return MinuteLocator(byminute=range(0, 60, 5)), MinuteLocator(byminute=range(0, 60, 2))

    elif (timedelta(minutes=10) < tdiff) & (tdiff <= timedelta(minutes=20)):
        return MinuteLocator(byminute=range(0, 60, 2)), MinuteLocator(byminute=range(0, 60, 1))

    elif (timedelta(minutes=5) < tdiff) & (tdiff <= timedelta(minutes=10)):
        return MinuteLocator(byminute=range(0, 60, 1)), SecondLocator(bysecond=range(0, 60, 30))

    elif (timedelta(minutes=1) < tdiff) & (tdiff <= timedelta(minutes=5)):
        return SecondLocator(bysecond=range(0, 60, 30)), SecondLocator(bysecond=range(0, 60, 10))

    elif (timedelta(seconds=30) < tdiff) & (tdiff <= timedelta(minutes=1)):
        return SecondLocator(bysecond=range(0, 60, 10)), SecondLocator(bysecond=range(0, 60, 2))

    else:
        return SecondLocator(bysecond=range(0, 60, 2)),  SecondLocator(bysecond=range(0, 60, 1))
