"""
"""
import sys
import os

from plt_exposure_irma import fig03_fig04
from plt_local_rp import fig05
from plt_analysis import fig06, fig07

DATA_DIR = os.path.dirname(__file__)
""" Input/output data folder relative path """

IBTRACS_DIR = os.path.join(DATA_DIR, 'tracks')
""" Tracks data in DATA_DIR """

FIG_DIR = DATA_DIR
""" Folder where the images are written """

def main(argv):
    print('Input/Output data folder: ', DATA_DIR)

    # FIG03 and FIG04
    fig03_fig04(IBTRACS_DIR, DATA_DIR, FIG_DIR) # 5min

    # FIG 05
    fig05(DATA_DIR, FIG_DIR)

    # FIG 06
    fig06(DATA_DIR, FIG_DIR)

    # FIG 07
    fig07(DATA_DIR, FIG_DIR)

if __name__ == "__main__":
   main(sys.argv[1:])
