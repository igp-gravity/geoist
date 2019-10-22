# -*- coding: utf-8 -*-
"""
Monkey-patch plt.figure() to support Ctrl+C for copying to clipboard as an image

@author: Josh Burnett
Modified from code found on Stack Exchange:
    https://stackoverflow.com/questions/31607458/how-to-add-clipboard-support-to-matplotlib-figures
    https://stackoverflow.com/questions/34322132/copy-image-to-clipboard-in-python3

"""


import matplotlib.pyplot as plt
from win32gui import GetWindowText, GetForegroundWindow
from PIL import Image
import win32clipboard
import io

__version__ = (2, 0, 0)
oldfig = plt.figure


def copyfig(fig=None):
    # store the image in a buffer using savefig(), this has the
    # advantage of applying all the default savefig parameters
    # such as background color; those would be ignored if you simply
    # grab the canvas
    if fig is None:
        # find the figure window that has UI focus right now (not necessarily the same as plt.gcf())
        fig_window_text = GetWindowText(GetForegroundWindow())
        for i in plt.get_fignums():
            if plt.figure(i).canvas.get_window_title() == fig_window_text:
                fig = plt.figure(i)
                break

    with io.BytesIO() as buf:
        fig.savefig(buf)
        im = Image.open(buf)

        with io.BytesIO() as output:
            im.convert("RGB").save(output, "BMP")
            data = output.getvalue()[14:]  # The file header off-set of BMP is 14 bytes

    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)  # DIB = device independent bitmap
    win32clipboard.CloseClipboard()


def newfig(*args, **kwargs):
    fig = oldfig(*args, **kwargs)

    def clipboard_handler(event):
        if event.key == 'ctrl+c':
            copyfig()

    fig.canvas.mpl_connect('key_press_event', clipboard_handler)
    return fig


plt.figure = newfig
