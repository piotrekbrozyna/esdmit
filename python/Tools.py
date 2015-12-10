import numpy as np
import math

class Tools(object):

    @staticmethod
    def zeroIfUnderThreshold(signal, threshold):
        for i, x in enumerate(signal):
            signal[i] = 0.0 if (x <= threshold) else x

    @staticmethod
    def integrateInMovingWindow(signal, samplingFrequency, tIntegrationWindow):
        nIntegrationWindow = int(math.floor(tIntegrationWindow * samplingFrequency))
        kernel = [1.0 / nIntegrationWindow] * nIntegrationWindow
        return np.convolve(signal, kernel, mode='same')

    # gives the same results as numpy.convolve
    @staticmethod
    def integrateInMovingWindow2(signal, samplingFrequency, tIntegrationWindow):
        signalLen = len(signal)
        nIntegrationWindow = int(math.floor(tIntegrationWindow * samplingFrequency))
        kernel = [1.0 / nIntegrationWindow] * nIntegrationWindow
        kernelLen = len(kernel)
        kmin = 0
        kmax = 0
        result = []
        for i in range(signalLen + kernelLen - 1):
            result.append(0.0)
            kmin = (i - (kernelLen - 1)) if (i >= kernelLen - 1) else 0
            kmax = i if (i < signalLen - 1) else (signalLen - 1)
            for j in range(kmin, kmax + 1):
                result[i] += signal[j] * kernel[i-j]

        return result[26:-27]

    @staticmethod
    def findMaximumWithinNeigborhood(signal, index, radius):
        signalSize = len(signal)
        start = end = 0
        if (index > radius) and (index + radius < signalSize):
            start = index - radius
            end = index + radius + 1
        elif (index - radius) < 0:
            start = 0
            end = index + radius + 1
        elif (index + radius) > signalSize:
            start = index - radius
            end = -1
        return start + np.argmax(signal[start:end])

    @staticmethod
    def adjust(signal, peaks, samplingFrequency, radius):
        nRadius = int(math.ceil(radius * samplingFrequency))
        adjustedPeaks = []
        for peak in peaks:
            newPos = Tools.findMaximumWithinNeigborhood(signal, peak, nRadius)
            adjustedPeaks.append(newPos)
            print '{} - {}'.format(peak, newPos)

        return adjustedPeaks

    @staticmethod
    def preallocateList(length, value=None):
        return [value] * length