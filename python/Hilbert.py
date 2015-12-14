import numpy as np
import scipy.signal as ss
import math
from Tools import Tools

class Hilbert(object):

    DIFF_KERNEL = [-0.125, -0.250, 0, 0.250, 0.125]
    RR_SIZE = 8

    def __init__(self, ecgSignal, tWindowSize=300.0, samplingFrequency=360.0):
        self.ecgSignal = np.copy(ecgSignal)
        self.ecgSignalSize = len(ecgSignal)
        self.samplingFrequency = samplingFrequency
        self.nWindowSize = int(tWindowSize * samplingFrequency)
        self.rPeaks = []

    # signalForEnvelope can be e.g. sum o IMFs
    def process(self, signalForEnvelope=None):

        if signalForEnvelope is None:
            signalForEnvelope = self.ecgSignal

        self.diffSignal = np.convolve(signalForEnvelope,
                                      Hilbert.DIFF_KERNEL,
                                      mode='same')
        self.analyticalSignal = ss.hilbert(self.diffSignal)
        self.envelope = np.abs(self.analyticalSignal)

        startIndex = 0

        while startIndex < self.ecgSignalSize:

            if startIndex + self.nWindowSize >= self.ecgSignalSize:
                output = self.doThresholding(self.ecgSignal[startIndex:],
                                             self.envelope[startIndex:])
            else:
                output = self.doThresholding(self.ecgSignal[startIndex:(startIndex + self.nWindowSize)],
                                             self.envelope[startIndex:(startIndex + self.nWindowSize)])

            if output:
                for item in output:
                    self.rPeaks.append(item + startIndex)

            startIndex += self.nWindowSize

        return self.rPeaks

    def doThresholding(self, ecgInterval, envelopeInterval):
        max = np.amax(envelopeInterval[:int(5 * self.samplingFrequency)])
        threshold = 0.8 * max
        output = [0]
        rr = [math.floor(0.7 * self.samplingFrequency)] * Hilbert.RR_SIZE
        tmpEcg = [0.0] * 5

        for i, data in enumerate(ecgInterval):
            dt = i - output[-1]
            rrMean = np.mean(rr)
            outputSize = len(output)

            if data > threshold:
                if dt <= 0.2 * self.samplingFrequency or dt <= 0.55 * rrMean:
                    if data > ecgInterval[output[-1]]:
                        output[-1] = i
                        tmpEcg[outputSize % 5] = ecgInterval[output[-1]]
                    if len(output) > 1:
                        rr[outputSize % 8] = output[-1] - output[-2]
                else:
                    output.append(i)
                    outputSize = len(output)
                    tmpEcg[outputSize % 5] = ecgInterval[output[-1]]
                    if outputSize > 1:
                        rr[outputSize % 8] = output[-1] - output[-2]

            if (outputSize > 4):
                threshold = ((sum(tmpEcg) - np.amax(tmpEcg)) / 4.0) * 0.55

        if output[0] == 0:
            output = np.delete(output, 0)
            print output[0]

        output = Tools.adjust(ecgInterval, output, self.samplingFrequency, 0.5)

        return output