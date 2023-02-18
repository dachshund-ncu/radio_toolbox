# Spectral toolbox
Set of simple tools to manipulate/fix/edit spectral data. At least python 3.8 required!

## Installing required packages
```bash
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

## Automatic data rotator
Simply rotates all of the .fits spectra in the directory, given as a script argument.

Usage:
```bash
python3 automatic_data_rotator.py directory_with_fits_files
```

Requirements:
- numpy
- astropy

Methodology: it reads the V_lsr of the oldest observation made in the direcrtory and then shifts all of the .fits files so all have same velocities in same channels.
### Important
```Automatic data rotator``` overwrites old .fits files so make sure you have a backup

## Header viewer
Displays a header info of the .fits file in the terminal window

Usage:
```bash
python3 header_viewer.py spectral_file.fits
```

## Wise photometry reader
Simple set of classes for easy reading of the crucial parts of the WISE satellite photometry tables, that are normally accesed from the [WISE catalog page](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?submit=Select&projshort=WISE).

Usage:
```python
# neowise
wise_data = neowise_reader()
wise_data.read("neowise_data.tbl")

# allwise
wise_data = allwise_reader()
wise_data.read("allwise_data.tbl")
```