#! /usr/bin/env python3
# -*- coding: utf-8 -*-

__doc__ = "This is an automated data rotator\n \
    It is meant only for rotating spectra, that have different\n\
    V_sys keyword in the .fits file header"

import numpy as np
from astropy.io import fits
from astropy.time import Time
import glob
import sys, os
import math

def display_welcome_message():
    '''
    Simply displays welcome message
    '''
    print("-------------------------------------")

def __get_flagged_filenames(directory: str) -> np.ndarray:
    '''
    searches for file, titled 'flagged_obs.dat', and 
    reads its contents
    '''
    if os.path.isfile(os.path.join(directory, 'flagged_obs.dat')):
        contents = np.loadtxt(os.path.join(directory, 'flagged_obs.dat'), dtype=str)
    else:
        return np.asarray([], dtype=str)
    try:
        len(contents)
        return contents
    except TypeError:
        return np.asarray([contents], dtype=str)

def __get_fits_filenames(directory: str, flagged_filenames: np.ndarray) -> np.ndarray:
    '''
    Returns non-flagged fits filenames in the given directory
    '''
    prelim_table = glob.glob(os.path.join(directory, '*.fits'))
    for i in range(len(prelim_table)-1, -1, -1):
        for flagged_name in flagged_filenames:
            if os.path.basename(prelim_table[i]) == flagged_name:
                prelim_table.pop(i)
                break
    return np.asarray(prelim_table, dtype=str)

def __load_epochs(filenames: np.ndarray) -> np.ndarray:
    '''
    Loads epochs of the given filenames
    '''
    mjds = []
    for f in filenames:
        hdul = fits.open(f)
        time = Time([hdul[1].header["TIME"]], format='isot', scale='utc')
        mjds.append(time.mjd[0])
    return np.asarray(mjds, dtype=float)

def __bubblesort(file_names: np.ndarray) -> np.ndarray:
    '''
    bubblesorts loaded files, require loading data from them
    '''
    epochs_table = __load_epochs(file_names)
    # -- main sorting loop --
    for i in range(len(epochs_table) - 1):
        for j in range(0, len(epochs_table) - i - 1):
            if epochs_table[j] > epochs_table[j+1]:
                epochs_table[j], epochs_table[j+1] = epochs_table[j+1], epochs_table[j]
                file_names[j], file_names[j+1] = file_names[j+1], file_names[j]
    return file_names

def __get_target_velocity(file_names: np.ndarray) -> float:
    '''
    Returns the target velocity from the first epoch
    amongst the loaded .fits files
    '''
    hdul = fits.open(file_names[0])[1].header
    return hdul["VSYS"]

def __rotate_with_interpolation(table: np.ndarray, no_of_channels: float, label: str)-> np.ndarray:
    '''
    rotates <<table>> with a floating point amount of channels (no_of_channels)
    '''
    # exceptions
    if no_of_channels >= 1.0:
        print(f"--> floating part of rotation cannot be higher, than 1!")
        return table # to not interrupt

    # getting over with rotation
    xT = [0, 1]
    rotated_spectrum_table = np.zeros(len(table))
    for i in range(len(table)):
        # -- getting a pair of points --
        if i+1 < len(table):
            yT = [table[i], table[i+1]]
        else:
            yT = [table[i], table[len(table) - (i+1)]]
        # -- fitting polynomial --
        p = np.polyfit(xT, yT, 1)
        rotated_spectrum_table[i] = np.polyval(p, -no_of_channels)
    return rotated_spectrum_table

def __rotate_fits_files(file_names: np.ndarray, target_velocity: float):
    '''
    Rotates fits files so the desired velocity lays in the middle 1024
    '''
    c = 299792.458 # km/s
    for f in file_names:
        # open fits file
        hdul = fits.open(f)
        # get spectral tables
        data = np.asarray(list(hdul[1].data))
        lhc = data[:,0][::-1]
        rhc = data[:,1][::-1]
        # calculate dChan
        # - read needed data -
        header = hdul[1].header
        restf = header["RESTFRQ"] / 1e6 # rest frequency (MHz)
        freq_rang = header["FRQ_RANG"] # bw (MHz)
        dv = header["VSYS"] - target_velocity
        # - calculate -
        nchans = len(lhc) # number of channels (IMPORTANT: assumes LHC and RHC are arrays of the same size)
        freqs_edges = np.asarray([restf - freq_rang/2.0, restf + freq_rang/2.0])
        vels = -c * ((np.asarray(freqs_edges) / restf) - 1.0)
        spectral_resolution = (vels[0] - vels[1]) / nchans
        freq_resolution = freq_rang / nchans

        # get the number of channels to rotate
        channels_to_rotate = dv / spectral_resolution
        df = freq_resolution * channels_to_rotate # for changing frequency

        # some exceptions
        if channels_to_rotate == 0.0:
            print(f"--> {os.path.basename(f)}: no rotation needed!")
            continue
        elif abs(channels_to_rotate) > nchans:
            print(f"--> \"{os.path.basename(f)}\": number of channels to rotate ({int(channels_to_rotate)}) exceeds the number of channels in the spectrum ({nchans}) - this is likely due to combining different sources in one direcrtory, which is not allowed.")
            break
        else:
            print(f"--> {os.path.basename(f)}: rotating by {round(channels_to_rotate,2)} channels")

        # - only positive rotations allowed - if it is negative, simply overflow the nChans
        # if channels_to_rotate < 0.0:
        #     channels_to_rotate = nchans + channels_to_rotate
        
        # split into integer and floating values
        # integer_rotate = int(channels_to_rotate) # integer part
        # floating_rotate = mod # floating-point part
        floating_rotate , integer_rotate = math.modf(channels_to_rotate)
        # rotation
        # - integer part -
        lhc = np.roll(lhc, int(integer_rotate))
        rhc = np.roll(rhc, int(integer_rotate))
        # - floating part -
        lhc = __rotate_with_interpolation(lhc, floating_rotate, f)
        rhc = __rotate_with_interpolation(rhc, floating_rotate, f)

        # -- this will be called in the final version
        # temporarily we save to new files
        columnPol1 = fits.Column(name='Pol 1', format='E', array=lhc[::-1])
        columnPol2 = fits.Column(name='Pol 2', format='E', array=rhc[::-1])
        pHeader = hdul[0]
        dHeader = fits.BinTableHDU.from_columns([columnPol1, columnPol2]) 
        dHeader.header = hdul[1].header
        dHeader.header["VSYS"] = target_velocity
        dHeader.header["FRQ_BEG"] = dHeader.header["FRQ_BEG"] + df
        dHeader.header["FRQ_MID"] = dHeader.header["FRQ_MID"] + df
        dHeader.header["FRQ_END"] = dHeader.header["FRQ_END"] + df
        hdul_n = fits.HDUList([pHeader, dHeader])
        hdul_n.writeto(f, overwrite=True)
def main():
    '''
    Body of this script
    '''
    flagged_filenames = __get_flagged_filenames(sys.argv[1])
    file_names = __get_fits_filenames(sys.argv[1], flagged_filenames)
    file_names = __bubblesort(file_names)
    target_velocity = __get_target_velocity(file_names)
    __rotate_fits_files(file_names, target_velocity)

if __name__ == '__main__':
    display_welcome_message()
    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]):
            main()
        else:
            print(f"--> {sys.argv[1]} is not a directory!")
    else:
        print(f"---> Usage: automatic_data_rotator.py directory_with_fits_files")
    print("-------------------------------------")