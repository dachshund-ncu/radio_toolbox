#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# -- Imports --
from astropy.io import fits
import sys, os
# ----------------------------------

def main():
    if len(sys.argv) < 2:
        print(f"Usage: python3 {os.path.basename(__file__)} spectrum.fits")
        return
    # -- opening .fits file --
    hdul = fits.open(sys.argv[1])
    print(repr(hdul[1].header))
    # ------------------------

if __name__ == '__main__':
    main()


