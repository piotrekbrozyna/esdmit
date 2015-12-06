import numpy as np
from numpy import array
import scipy.interpolate as spi
import math
from Hilbert import Hilbert

class EMD(object):

    SD_FACTOR = 0.3

    def __init__(self, ecgSignal, samplingFrequency=360.0, imfsNo=3):
        self.ecgSignal = array(ecgSignal)
        self.ecgSignalSize = len(ecgSignal)
        self.samplingFrequency = samplingFrequency
        self.rPeaks = []
        self.imfsNo = imfsNo

    def process(self):
        self.imfs = self.findIMFs(self.imfsNo)
        self.hilbertAlgorithm = Hilbert(np.sum(self.imfs, axis=0), self.samplingFrequency)
        # self.hilbertAlgorithm = Hilbert(self.imfs, self.samplingFrequency)
        self.rPeaks = self.hilbertAlgorithm.process()
        return self.rPeaks

    def findIMFs(self, n=3):
        part_size = 1000
        max_parts = int(math.ceil(float(self.ecgSignalSize) / float(part_size)))
        # imfs = np.array([[],[],[]])
        # imfs = np.array([])
        imfs = [[] for x in xrange(n)]

        for p in range(0, max_parts * part_size, part_size):

            N = self.ecgSignalSize - p if p == (max_parts - 1) * part_size else part_size

            x_part = np.copy(self.ecgSignal[p:(p+N)])
            c = np.copy(x_part)

            imfs_parts = np.empty([n, N])

            for i in range(0, n):
                h = np.copy(c)
                sd = 1
                maxmin_length = 0

                prev_h = np.empty(N)

                while sd > EMD.SD_FACTOR:
                    maxima, minima = self.findExtrema(x_part, addBoundaries=True)

                    # is it already a residue?
                    maxmin_length = len(maxima) + len(minima)
                    if maxmin_length < 2:
                        break

                    maxima_val = [x_part[j] for j in maxima]
                    minima_val = [x_part[j] for j in minima]

                    upper_bound_func = spi.interp1d(maxima, maxima_val, kind='cubic', copy=True)
                    lower_bound_func = spi.interp1d(minima, minima_val, kind='cubic', copy=True)

                    x_new = np.arange(0, N, 1)
                    upper_bound = upper_bound_func(x_new)
                    lower_bound = lower_bound_func(x_new)

                    m = np.add(lower_bound, upper_bound) * 0.5

                    prev_h = np.copy(h)
                    h = x_part - m

                    sd = self.calculateSD(prev_h, h, eps=0.0000001)

                imfs_parts[i] = h

                if maxmin_length < 2:
                    break

                c = np.subtract(c, h)

            imfs = np.append(imfs, imfs_parts, axis=1)
            # self.imfs = np.append(self.imfs, imfs_parts)

        return imfs

    def calculateSD(self, previous, current, eps=0.0):
        eps_array = np.empty(len(current))
        eps_array.fill(eps)
        # (previous-current)^2/(current^2+eps)
        numerator = np.square(np.subtract(previous, current))
        denominator = np.add(np.square(current), eps_array)
        summands = np.divide(numerator, denominator)
        return np.sum(summands)

    def findExtrema(self, signal, addBoundaries=True):

        diffSignal = np.diff(signal)

        maxima = np.array([])
        minima = np.array([])

        for i in range(1, len(diffSignal) - 1):
            if diffSignal[i] == 0:
                if np.sign(diffSignal[i-1]) == 1 and np.sign(diffSignal[i+1]) == -1:
                    maxima = np.append(maxima, i)
                elif np.sign(diffSignal[i-1]) == -1 and np.sign(diffSignal[i+1]) == 1:
                    minima = np.append(minima, i)
            elif np.sign(diffSignal[i]) == 1 and np.sign(diffSignal[i+1]) == -1:
                maxima = np.append(maxima, i + 1)
            elif np.sign(diffSignal[i]) == -1 and np.sign(diffSignal[i+1]) == 1:
                minima = np.append(minima, i + 1)

        if addBoundaries:
            maxima = np.insert(maxima, 0, 0)
            minima = np.insert(minima, 0, 0)
            maxima = np.append(maxima, len(signal) - 1)
            minima = np.append(minima, len(signal) - 1)

        return maxima, minima


