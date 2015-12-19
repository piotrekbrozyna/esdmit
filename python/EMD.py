import numpy as np
import scipy.interpolate as spi
import scipy.signal as sps
from Hilbert import Hilbert

class EMD(object):

    SD_LIMIT_VAL = 0.3

    def __init__(self, ecgSignal, samplingFrequency=360.0, tWindowSize=120.0, imfsNo=1):
        self.ecgSignal = np.array(ecgSignal)
        self.ecgSignalSize = len(ecgSignal)
        self.samplingFrequency = samplingFrequency
        self.tWindowSize = tWindowSize
        self.nWindowSize = int(tWindowSize * samplingFrequency)
        self.rPeaks = []
        self.imfsNo = imfsNo

    def process(self):
        self.imfs = self.findIMFs(self.imfsNo)
        self.hilbertAlgorithm = Hilbert(self.ecgSignal, self.tWindowSize, self.samplingFrequency)
        self.imfsSum = np.sum(self.imfs, axis=0)
        self.rPeaks = self.hilbertAlgorithm.processImfs(self.imfsSum)
        return self.rPeaks

    def findIMFs(self, n=3):
        imfs = np.empty([n, self.ecgSignalSize])

        residue = False

        N = self.ecgSignalSize

        x_part = np.copy(self.ecgSignal[:N])
        c = np.copy(x_part)

        for i in range(0, n):
            h = np.copy(c)
            sd = np.inf

            prev_h = np.copy(h)

            while sd > EMD.SD_LIMIT_VAL:
                maxima, minima = self.findExtrema(h)

                # is it already a residue?
                if len(maxima) + len(minima) < 2:
                    residue = True
                    break

                maxima = np.insert(maxima, 0, 0)
                minima = np.insert(minima, 0, 0)
                maxima = np.append(maxima, len(h) - 1)
                minima = np.append(minima, len(h) - 1)

                maxima_val = [h[j] for j in maxima]
                minima_val = [h[j] for j in minima]

                upper_bound_func = spi.InterpolatedUnivariateSpline(maxima,
                                                                    maxima_val,
                                                                    k=3)
                lower_bound_func = spi.InterpolatedUnivariateSpline(minima,
                                                                    minima_val,
                                                                    k=3)

                x_new = np.arange(0, N, 1)
                upper_bound = upper_bound_func(x_new)
                lower_bound = lower_bound_func(x_new)

                m = np.add(lower_bound, upper_bound) * 0.5

                h = np.subtract(prev_h, m)

                sd = self.calculateSD(prev_h, h, eps=0.0000001)
                prev_h = np.copy(h)

            imfs[i] = h

            if residue:
                break

            c = np.subtract(c, h)

        return imfs

    def calculateSD(self, previous, current, eps=0.0):
        eps_array = np.empty(len(current))
        eps_array.fill(eps)
        numerator = np.sum(np.square(np.subtract(previous, current)))
        denominator = np.sum(np.add(np.square(previous), eps_array))
        return np.divide(numerator, denominator)

    def findExtrema(self, signal):
        maxima = sps.argrelextrema(signal, np.greater)
        minima = sps.argrelextrema(signal, np.less)

        return maxima, minima
