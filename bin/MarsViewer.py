#!/usr/bin/env python3
"""
The MarsViewer executable is a pure-Python implementation of a \
primitive image viewer.

It is designed to display an image preview in the terminal where accessing NAS through X11 is not possible or is slow. Supported formats are PNG, EPS, and PDF.

Adapted from:
https://github.com/Belval/pdf2image

The executable requires no arguments.

Third-party Requirements:
    * numpy
    * Image

List of Functions:
    * get_ansi_color_code - 
    * get_color - 
    * show_image - 
"""

# load generic Python modules
import sys          # system commands
import os           # access operating system functions
import numpy as np
from PIL import Image

# load amesCAP modules
from amescap.pdf2image import convert_from_bytes

# ======================================================
#                  DEFINITIONS
# ======================================================

def get_ansi_color_code(r, g, b):
    """
    -

    -

    Parameters
    ----------
    r : float
        -
    g : float
        -
    b : float
        -
    
    Raises
    ------

    Returns
    -------
    -
        -
    """
    if r == g and g == b:
        if r < 8:
            return 16
        elif r > 248:
            return 231
        else:
            return round(((r - 8) / 247) * 24) + 232
    else:
        return 16 + (36 * round(r / 255 * 5)) + (6 * round(g / 255 * 5)) + round(b / 255 * 5)

def get_color(r, g, b):
    """
    -

    -

    Parameters
    ----------
    r : float
        -
    g : float
        -
    b : float
        -
    
    Raises
    ------

    Returns
    -------
    -
        -
    """
    return f"\x1b[48;5;{int(get_ansi_color_code(r, g, b))}m \x1b[0m"

def show_image(img_path):
    """
    -

    -

    Parameters
    ----------
    img_path : str
        -
    
    Raises
    ------

    Returns
    -------
    -
        -
    """
    try:
        # If PDF, merge images together
        if img_path[-3:] == 'pdf':
            print('Extracting pdf pages...')
            images = convert_from_bytes(open(img_path, 'rb').read())
            widths, heights = zip(*(i.size for i in images))
            total_height = sum(heights)
            max_width = max(widths)

            img = Image.new('RGB', (max_width, total_height))

            y_offset = 0
            for im in images:
                img.paste(im, (0, y_offset))
                y_offset += im.size[1]
        else:
            img = Image.open(img_path)
    except Exception as e:
        exit(e)
    
    rows, columns = os.popen('stty size', 'r').read().split()
    mywidth = int(columns)
    scale_factor = 0.5  # account for height of line in terminal

    wpercent = (mywidth/float(img.size[0]))
    hsize = int(scale_factor*(float(img.size[1])*float(wpercent)))
    img = img.resize((mywidth, hsize), Image.ANTIALIAS)

    img_arr = np.asarray(img)
    h, w, c = img_arr.shape

    for x in range(h):
        for y in range(w):
            pix = 1.*img_arr[x][y]
            txt_col = get_color(pix[0], pix[1], pix[2])
            sys.stdout.write(txt_col)
        print('')

# ======================================================
#                  MAIN PROGRAM
# ======================================================
print(sys.argv[0])
print(sys.argv[1])

def main():
    if len(sys.argv) > 1:
        img_path = sys.argv[1]
        show_image(img_path)

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
