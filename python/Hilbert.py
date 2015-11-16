import numpy as np
import scipy.signal as ss
from Tools import Tools

class Hilbert(object):

    DIFF_KERNEL = [-0.125, -0.250, 0, 0.250, 0.125]

    def __init__(self, ecgSignal, samplingFrequency=1000.0):
        self.ecgSignal = ecgSignal
        self.ecgSignalSize = len(ecgSignal)
        self.samplingFrequency = samplingFrequency
        self.rPeaks = []

    def process(self):
        self.dSignal = self.differentiateEcgSignal()
        self.hSignal = ss.hilbert(self.dSignal)
        self.calculateEnvelope()

        startIndex = 0
        intervalSize = 300 * self.samplingFrequency
        while startIndex < self.ecgSignalSize:
            if startIndex + intervalSize >= self.ecgSignalSize:
                self.threshold(startIndex, -1)
            else:
                self.threshold(startIndex, startIndex + intervalSize)
            startIndex += intervalSize

    def threshold(self, startIndex, endIndex):
        max = np.amax(startIndex, startIndex + 5 * self.samplingFrequency)
        threshold = 0.8 * max
        self.rPeaks.append(0)


    def differentiateEcgSignal(self):
        return np.convolve(self.ecgSignal, Hilbert.DIFF_KERNEL, mode='same')

    def calculateEnvelope(self):
        self.envelope = np.sqrt(np.square(self.dSignal) + np.square(self.hSignal))