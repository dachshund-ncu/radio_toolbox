'''
This is a definition fo a class, that reads .IPAC table with WISE photometric data
and stores them in easy-to-read tables
Data is also converted to flux density from MAG
'''

# -- importing necessary modules --
import numpy as np
from astropy.io import ascii
import math
import abc
# ---------------------------------

class photWiseData(abc.ABC):
    def __init__(self):
        '''
        Initializes the class instance and reads file contents upon initialization
        '''

    @abc.abstractmethod
    def read(self, filename: str):
        '''
        This method simply reads from file, provided in the
        << filename >> argument
        '''
        return
    
    @abc.abstractmethod
    def filter_by_flags(self, *args):
        '''
        This method aims to filter the data using flags, provided in the ipac table
        '''

    def flux_from_wise(self, mag, band):
        '''
        This is simple method for converting from MAG to Jy
        "BAND" can be only w(1-4)
        Formula comes from https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2flux
        '''
        if band == "w1":
            Fv0 = 309.540
        elif band == "w2":
            Fv0 = 171.787
        elif band == "w3":
            Fv0 = 31.674
        elif band == "w4":
            Fv0 = 8.363
        else:
            print("Did not use correct band for flux convertion!")
            return 0.0

        # - we calculate flux density -
        Fv = Fv0 * 10.0 ** (-mag / 2.5)

        # - and return corrected value -
        return Fv

    def remove_nans_from_arrays(self, verbose=False):
        '''
        This method simply removes np.nan values from our beautiful arrays
        In order to avoid problems, we need to iterate through arrays backwards
        '''
        # print message v.1
        if verbose:
            print('# ---- FLAGGING RESULTS ----')
            print(f'# Before: W1:{len(self.W[0])}, W2:{len(self.W[1])}')

        for indexWin in range(len(self.W)-1, -1, -1): # from W4 to W1
            for indexItem in range(len(self.W[indexWin])-1, -1, -1 ): # from end of W to the begin of W
                if math.isnan(self.W[indexWin][indexItem]):
                    self.W[indexWin] = np.delete(self.W[indexWin], indexItem)
                    self.W_error[indexWin] = np.delete(self.W_error[indexWin], indexItem)
                    self.mjd[indexWin] = np.delete(self.mjd[indexWin], indexItem)
                    self.Jy[indexWin] = np.delete(self.Jy[indexWin], indexItem)
        
        # print message v.2
        if verbose:
            print(f'# After:  W1:{len(self.W[0])}, W2:{len(self.W[1])}')
            print('# --------------------------')

    def flag_item(self, indexWin, indexItem):
        '''
        replaces item in self.magTab[indexWin][indexItem] witn np.nan
        '''
        self.W[indexWin][indexItem] = np.nan
        self.W_error[indexWin][indexItem] = np.nan
        self.mjd[indexWin][indexItem] = np.nan
        self.Jy[indexWin][indexItem] = np.nan
    
    def extract_wise_epochs(self, mjdTab, binWidth):
        '''
        extracts mean from wise epochs based on the binWidth 
        '''
        mjdReturn = [mjdTab[0]]

        for date in mjdTab:
            newEpochFlag = False

            for e in mjdReturn:
                if date < e - (binWidth / 2.0) or date > e + (binWidth / 2.0):
                    pass
                else:
                    newEpochFlag = True
                    
            if not newEpochFlag:
                mjdReturn.append(date)

        return np.asarray(mjdReturn)

    def get_mean_epoch_values(self, epochs, mjd, mag, flux, binWidth):
        '''
        Assumption: len(mjd) == len(mag) == len(flux)
        this simply finds mean for photometric observations gathered around "epochs"
        '''
        mjdReturn = []
        magReturn = []
        magErrReturn = []
        fluxReturn = []
        fluxErrReturn = []
        
        for border in epochs:
            tmpMjd = []
            tmpMag = []
            tmpFlux = []
            for i, date in enumerate(mjd):
                if date > border - binWidth / 2.0 and date < border + binWidth / 2.0:
                    tmpMjd.append(date)
                    tmpMag.append(mag[i])
                    tmpFlux.append(flux[i])
            if len(tmpMjd) > 0:
                mjdReturn.append(np.mean(tmpMjd))
                magReturn.append(np.mean(tmpMag))
                magErrReturn.append(np.std(tmpMag))
                fluxReturn.append(np.mean(tmpFlux))
                fluxErrReturn.append(np.std(tmpFlux))
        
        return np.asarray(mjdReturn), np.asarray(magReturn), np.asarray(magErrReturn), np.asarray(fluxReturn), np.asarray(fluxErrReturn)
    
    def get_mean_epochs(self, binWidth = 60.0):
        '''
        This method returns mean values for every WISE epoch
        WISE epoch is here defined as a set of photometric measurements, performed
        while source was inside the WISE field-of-view
        Typically there is one epoch per half a year, but due to flagging etc. we might
        be forced to omit some of them - thus, this method will generate separate arrays
        for every band (W1, W2, W3, W4)
        we will return:
        mjdMeanTab: [w1MjdMeanTab, w2MjdMeanTab, w3MjdMeanTab, w4MjdMeanTab]
        magMeanTab: [w1MagMeanTab, w2MagMeanTab, w3MagMeanTab, w4MagMeanTab]
        magErrTab: [w1MagErrTab, w2MagErrTab, w3MagErrTab, w4MagErrTab]
        jyMeanTab: [w1JymeanTab, w2JymeanTab, w3JymeanTab, w4JymeanTab]
        jyErrTab: [w1JyErrTab, w2JyErrTab, w3JyErrTab, w4JyErrTab]
        '''
        # --- returned arrays ---
        mjdMeanTab = []
        magMeanTab = []
        magErrTab = []
        jyMeanTab = []
        jyErrTab = []
        # -----------------------
        # --- main loop ---
        for index, window in enumerate(self.W):
            epochs = self.extract_wise_epochs(self.mjd[index], binWidth)
            mjd_tmp, mag_tmp, mag_err_tmp, jy_tmp, jy_err_tmp = self.get_mean_epoch_values(epochs, self.mjd[index], window, self.Jy[index], binWidth)
            mjdMeanTab.append(mjd_tmp)
            magMeanTab.append(mag_tmp)
            magErrTab.append(mag_err_tmp)
            jyMeanTab.append(jy_tmp)
            jyErrTab.append(jy_err_tmp)
        return np.asarray(mjdMeanTab), np.asarray(magMeanTab), np.asarray(magErrTab), np.asarray(jyMeanTab), np.asarray(jyErrTab)

class neowise_reader(photWiseData):
    def __init__(self):
        '''
        Initializes the neowise reader class
        '''
        super().__init__()
    
    def read(self, filename: str):
        '''
        reads ASCII file using ASTROPY.IO package and stores its contents
        in the suitable arrays
        temporary arrays are garbage collected after this method ends
        '''
        data = ascii.read(filename)
        # magnitudes
        self.W = []
        self.W_error = []
        self.Jy = []
        self.mjd = []
        mjd = np.array(data.columns['mjd'])
        for i in range(2):
            self.W.append(np.array(data.columns[f"w{i+1}mpro"]))
            self.W_error.append(np.array(data.columns[f"w{i+1}sigmpro"]))
            self.Jy.append(self.flux_from_wise(self.W[i], f"w{i+1}"))
            self.mjd.append(mjd.copy())
        # flags
        self.ccFlags = list(data.columns['cc_flags'])
        self.qualFrame = list(data.columns['qual_frame'])
        self.phQual = list(data.columns['ph_qual'])
    
    def filter_by_flags(self, ccFlags = [], qualFrame = 0, phQual = [], verbose=False):
        '''
        This is mwthod that erases broken data, based on quality flags
        proivded with photometry results
        Below methods simply fill places in arrays with np.nan
        and remove_nans_from_arrays() removes these
        '''
        if len(ccFlags) != 0:
            self.filter_cc_flags(ccFlags)
        if qualFrame != 0:
            self.filter_qual_frame(qualFrame)
        if len(phQual) != 0:
            self.filter_ph_qual(phQual)
        
        self.remove_nans_from_arrays(verbose)
    
    def filter_cc_flags(self, ccFlags):
        '''
        if input in self.ccFlags does not mach any input in ccFlags
        this element of array is replaced with np.nan
        '''
        ccFlags = list(map(str.upper, ccFlags))
        for indexWin, W in enumerate(self.W):
            for itemIndex, flag in enumerate(self.ccFlags):
                actualFlag = flag[indexWin].upper()
                if actualFlag not in ccFlags:
                    self.flag_item(indexWin, itemIndex)

    def filter_qual_frame(self, qualFrame):
        '''
        if input in self.qualFrame does not mach any input in qualFrame
        this element of array is replaced with np.nan
        qualFrame -> flags accepted to be ok [5] etc
        self.qualFrame -> flags of the data: 0, 5 or 10
        '''
        for indexWin, W in enumerate(self.W):
            for itemIndex, flag in enumerate(self.qualFrame):
                if flag < qualFrame:
                    self.flag_item(indexWin, itemIndex)

    def filter_ph_qual(self, phQual):
        '''
        if input in self.phQual does not mach any input in phQual
        this element of array is replaced with np.nan
        '''
        # making it easier by uppering every letter
        phQual = list(map(str.upper, phQual))

        # double nested for loops: 1: window 2: epoch
        for indexWin, W in enumerate(self.W):

            # phqual tends to be shorter, than 4 windows. Thus we introduce safety feature:
            if indexWin > len(self.phQual[0])-1:
                break
            
            # easy to understand: if flag does not mach any of the proposed flags in phQual (argument), then this item gets converted to np.nan
            for itemIndex, flag in enumerate(self.phQual):
                actualFlag = flag[indexWin].upper()
                if actualFlag not in phQual:
                    self.flag_item(indexWin, itemIndex)

class allwise_reader(photWiseData):
    def __init__(self):
        '''
        Initializes the neowise reader class
        '''
        super().__init__()
    
    def read(self, filename: str):
        '''
        Reads data from the ipac table
        '''
        data = ascii.read(filename)
        # magnitudes
        self.W = []
        self.W_error = []
        self.Jy = []
        self.mjd = []
        mjd = np.array(data.columns['mjd'])
        for i in range(2):
            self.W.append(np.array(data.columns[f"w{i+1}mpro_ep"]))
            self.W_error.append(np.array(data.columns[f"w{i+1}sigmpro_ep"]))
            self.Jy.append(self.flux_from_wise(self.W[i], f"w{i+1}"))
            self.mjd.append(mjd.copy())
        # epochs
        mjd = np.array(data.columns['mjd'])
        # flags
        self.qi_fact = list(data.columns['qi_fact'])

    def filter_by_flags(self, min_value: float, verbose=False):
        '''
        This is mathod that erases broken data, based on quality flags
        proivded with photometry results
        Below methods simply fill places in arrays with np.nan
        and remove_nans_from_arrays() removes these
        '''
        self.filter_qi_fact(min_value)
        self.remove_nans_from_arrays(verbose)

    def filter_qi_fact(self, maxValue: float):
        '''
        if input in self.qi_fact does not mach any input in qualFrame
        this element of array is replaced with np.nan
        qualFrame -> flags accepted to be ok [5] etc
        self.qualFrame -> flags of the data: 0, 5 or 10
        '''
        for indexWin, W in enumerate(self.W):
            for itemIndex, flag in enumerate(self.qi_fact):
                if flag < maxValue:
                    self.flag_item(indexWin, itemIndex)
