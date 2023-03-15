__doc__ = "This file contains simple set of classes that allow simple fits_files reading. There are two classes:\
    - Spectrum, which allows for reading single fits file\
    - setOfSpec, which reads a bunch of the .fits files\
For performance reasons, classes below utilize fitsio package for reading fits files."

import fitsio as fits
import numpy as np
from astropy.time import Time
import os
import pandas as pd
from glob import glob

class Spectrum:
    def __init__(self, filename, number = 0):
        '''
        initializes the class instance
        '''
        self.filename = filename
        self.readDataFromHeader(self.__readDataHeader(filename))
        
    def __readDataHeader(self, filename):
        '''
        reads and returns data header
        '''
        return fits.FITS(filename)[1]
    
    def readDataFromHeader(self, dataSectionOfFitsFile):
        '''
        reads data from FITS file header
        assumes, that data header was read before
        '''
        header = dataSectionOfFitsFile.read_header()
        self.lhcTab, self.rhcTab = self.__readPolData(dataSectionOfFitsFile.read())
        # --- reading key values ---
        self.sourcename = header['OBJECT']
        # -- DOPPLER TRACKING --
        self.Vlsr = header['VSYS'] # systemic velocity
        self.__freqRang = header['FRQ_RANG'] # frequency range
        self.__restFreq = header['FREQ'] / 1000000.0 # rest frequency (in FITS file it is in Hz, we convert it to MHz)

        try: # TSYS
            self.__tsys1 = header['TSYS1'] 
            self.__tsys2 = header['TSYS2']
        except:
            self.__tsys1 = header['TSYS']
            self.__tsys2 = header['TSYS']
        self.tsys = np.mean([self.__tsys1, self.__tsys2])

        # ADDITIONAL INFORMATIONS (should be mostly 0.0)
        self.__dopp_vto = header['DOPP_VTO']

        # -- DATE AND TIME --
        self.isotime = header['DATE-OBS']
        t = Time(self.isotime, format='isot', scale='utc')
        self.mjd = t.mjd

        del t # to save memory

        # -- IV STOKES PARAMS --
        self.iTab, self.vTab = self.__makeIVTabs(self.lhcTab, self.rhcTab)
        # -- RMS --
        self.rmsLhc = self.__calculateRMS(self.lhcTab, [25, 300], [-300, -25])
        self.rmsRhc = self.__calculateRMS(self.rhcTab, [25, 300], [-300, -25])
        self.rmsIhc = self.__calculateRMS(self.iTab, [25, 300], [-300, -25])
        self.rmsVhc = self.__calculateRMS(self.vTab, [25, 300], [-300, -25])
        # -- doppler tracking --
        self.velocityTable = self.__generateVeltab(self.Vlsr, self.__dopp_vto, self.__restFreq, self.__freqRang)
        

    def __generateVeltab(self, Vlsr, dopp_vto, restFreq, freqRang):
        '''
        generates doppler-tracked velocity table
        '''
        full_velocity = Vlsr + dopp_vto
        c = 299792.458 # km/s
        beta = full_velocity / c
        gamma = 1.0 / np.sqrt(1.0 - beta*beta)
        fCentr = restFreq * (gamma * (1.0 - beta))
        fBeg = fCentr - (freqRang / 2.0)
        # --
        freqs = np.linspace(fBeg, fBeg + freqRang, len(self.lhcTab))
        vels = -c * ( (freqs / restFreq) - 1.0)
        # --
        vels = vels[::-1] # reversing velocities array to get it in increasing order
        return vels

    def __calculateRMS(self, tab, minLims, maxLims):
        '''
        calculates RMS based on data within specified range
        '''
        suma = 0.0
        #lhc
        suma += np.sum(tab[minLims[0]:minLims[1]] * tab[minLims[0]:minLims[1]] )
        suma += np.sum(tab[maxLims[0]:maxLims[1]] * tab[maxLims[0]:maxLims[1]])
        suma /= (abs(minLims[1] - minLims[0]) + abs(maxLims[1] - maxLims[0])) - 1.0
        rms = np.sqrt(suma)
        return rms
        
    def __readPolData(self, header):
        '''
        reads polarization data from header
        also zeroes indexes 0:25, -25:0
        This is due to fact, that bandpass characteristics sometimes escapes to very high values
        on the edge of the spectrum
        '''
        self.lhcTab = np.asarray(header['Pol 1'])[::-1] # reversing to get velocity in ascending order
        self.rhcTab = np.asarray(header['Pol 2'])[::-1]
        self.lhcTab[0:25] =  0.0#np.zeros(25)
        self.lhcTab[-25:-1] = 0.0#np.zeros(25)
        self.rhcTab[0:25] = 0.0#np.zeros(25)
        self.rhcTab[-25:-1] = 0.0#np.zeros(25)

        self.rhcTab[-1] = 0.0
        self.lhcTab[-1] = 0.0
        return self.lhcTab, self.rhcTab

    def __makeIVTabs(self, lhcTab, rhcTab):
        '''
        Makes the I and V stokes parameters tables, based on input of the LHC and RHC tables data
        '''
        I = (lhcTab + rhcTab) / 2.0
        V = (rhcTab - lhcTab) / 2.0
        return I,V

    def get_dataframe(self) -> pd.DataFrame:
        '''
        returns a pandas dataframe with this spectrum
        '''
        return pd.DataFrame(np.column_stack((self.velocityTable, self.iTab, self.lhcTab, self.rhcTab, self.vTab)), columns=["Velocity", "I", "LHC", "RHC", "V"])

    def get_integrated_flux_density(self, min_chan: int, max_chan: int) -> np.ndarray:
        '''
        Returns the integrated flux density of the obs, based on min and max channels
        '''
        channels = np.asarray(range(1,len(self.iTab)+1))
        indices = np.logical_and(channels > min_chan, channels < max_chan)
        Velocity = self.velocityTable[indices]

        I = np.trapz(self.iTab[indices], Velocity)
        V = np.trapz(self.vTab[indices], Velocity)
        LHC = np.trapz(self.lhcTab[indices], Velocity)
        RHC = np.trapz(self.rhcTab[indices], Velocity)
        return np.asarray([self.mjd, I,V,LHC,RHC])

    def __str__(self):
        return repr(self.mjd)


class setOfSpec:
    def __init__(self, catWithSource):
        '''
        initializes the class
        '''
        self.dataCatalog = catWithSource
        self.flagged_obs = self.__read_flagged_obs(self.dataCatalog)
        self.spectra = self.__loadSpectraFromCat(self.dataCatalog)
        
        self.__bubblesort()

    def __loadSpectraFromCat(self, catWithSource):
        '''
        loads all the .FITS files from given directory
        '''
        e = os.walk(catWithSource, topdown=True)
        fitsFiles = []
        for root, dir, file in e:
            for filename in file:
                if filename.endswith('.fits') or filename.endswith('.FITS'):
                    if os.path.basename(filename) not in self.flagged_obs:
                        fitsFiles.append(os.path.join(catWithSource, filename))
        return [Spectrum(fitsFileName, no) for no, fitsFileName in enumerate(fitsFiles)]

    def __bubblesort(self):
        '''
        bubble-sorts spectra
        '''
        # -- main sorting loop --
        for i in range(len(self.spectra) - 1):
            for j in range(0, len(self.spectra) - i - 1):
                if self.spectra[j].mjd > self.spectra[j+1].mjd:
                    self.spectra[j], self.spectra[j+1] = self.spectra[j+1], self.spectra[j]

    def __read_flagged_obs(self, directory: str) -> np.ndarray:
        '''
        Reads the flagged filenames
        '''
        try:
            return np.loadtxt(os.path.join(directory, 'flagged_obs.dat'), dtype=str)
        except FileNotFoundError:
            return []
        
    def __str__(self):
        if len(self.spectra) > 0:
            return f"{len(self.spectra) }"
        else:
            return "None"
    # =======================
    # === UTILITY METHODS ===
    # =======================

    def get2DdataArray(self, pol = 'I'):
        '''
        returns the 2D container, containing data from the loaded spectra
        '''
        if len(self.spectra) == 0:
            raise BufferError("No data loaded!")
        else:
            z_array = np.empty( (len(self.spectra[0].iTab), len(self.spectra)) )
            # -- iterating to get the data
            for i in range(len(self.spectra)):
                if pol == 'V':
                    z_array[:,i] = self.spectra[i].vTab
                elif pol == 'LHC':
                    z_array[:,i] = self.spectra[i].lhcTab
                elif pol == 'RHC':
                    z_array[:,i] = self.spectra[i].rhcTab
                else:
                    z_array[:,i] = self.spectra[i].iTab
            return z_array
    
    def getVelArray(self):
        '''
        returns the VELOCITY table from the first spectrum
        '''
        if len(self.spectra) == 0:
            raise BufferError("No data loaded!")
        else:
            return self.spectra[0].velocityTable

    def getMjdArray(self):
        '''
        returns the dates array from the loaded dataset
        '''
        if len(self.spectra) == 0:
            raise BufferError("No data loaded!")
        else:
            return np.asarray([s.mjd for s in self.spectra])
    
    def get_mean_spectrum(self):
        '''
        Returns the mean spectrum as a data frame
        '''
        velocity = self.getVelArray()
        I = np.mean(np.asarray([sp.iTab for sp in self.spectra]), axis=0)
        V = np.mean(np.asarray([sp.vTab for sp in self.spectra]), axis=0)
        LHC = np.mean(np.asarray([sp.lhcTab for sp in self.spectra]), axis=0)
        RHC = np.mean(np.asarray([sp.rhcTab for sp in self.spectra]), axis=0)
        return pd.DataFrame(np.column_stack( (velocity, I, V, LHC, RHC)), columns=["Velocity", "I", "V", "LHC", "RHC"] )
    
    def get_integrated_flux_density(self, min_chan: int, max_chan:int, df=False) -> np.ndarray:
        '''
        Returns the integrated flux density for whole time series
        '''
        array = np.asarray([sp.get_integrated_flux_density(min_chan, max_chan) for sp in self.spectra])
        if df:
            return pd.DataFrame(array, columns=["MJD", "I", "V", "LHC", "RHC"])
        else:
            return array

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    cat = os.path.dirname(__file__)
    data = setOfSpec(os.path.join(cat, 'sample_data'))
    plt.pcolormesh(data.getMjdArray(), data.getVelArray() ,data.get2DdataArray(pol='I'), cmap='jet')
    plt.show()