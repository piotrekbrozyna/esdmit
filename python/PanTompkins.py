import numpy
import math

class PanTompkins(Class):

    RR_LOW_LIMIT_FACTOR = 0.92
    RR_HIGH_LIMIT_FACTOR = 1.16
    RR_SIZE = 8

    def __init__(self, ecgData, samplingFrequency, numOfSamplesPerWindow):
        self.ecgData = ecgData
        self.samplingFrequency = samplingFrequency
        self.integrationWindowTime = float(numOfSamplesPerWindow) / samplingFrequency
        self.intervalCount = math.floor(samplingFrequency * integrationWindowTime)

    def process(self):
        differentiatedSignal = self.differentiate(self.ecgData)
        squaredSignal = []
        numpy.square(differentiatedSignal, squaredSignal)
        self.integratedSignal = self.integrateInWindow(squaredSignal)
        self.rpeaks = detectRPeaks(ecgData, self.integratedSignal, gradientOfIntegratedSignal)

    def integrateInWindow(self, signal):
        kernel = [1.0 / self.intervalCount] * self.intervalCount
        return numpy.convlove(signal, kernel, mode='full')

    def differentiate(self, signal):
        kernel = [-0.125, -0.250, 0, 0.250, 0.125]
        return numpy.convolve(ecgData, kernel, mode='full')

    def detectRPeaks(self):
        for i, x in enumerate(self.integratedSignal):
            signal[i] = (x[i] <= 0.0001) ? 0.0 : x[i]

        self.fiducialMarkers = self.findFiducialMarkers()
        self.initializeThresholds()
        self.threshold()

    def threshold(self):
        windowTime = 0.1
        samplesPerWindow = self.samplingFrequency * windowTime
        output = []
        self.rr = 0.0
        self.rrMean = 0.0
        for i, marker in enumerate(fiducialMarkers):
            if marker > samplesPerWindow && (marker + samplesPerWindow) < len(self.ecgData):
                pos = numpy.amax(self.ecgData[(marker-samplesPerWindow):(marker+samplesPerWindow)])
            elif i == 0:
                pos = numpy.amax(self.ecgData[:(marker+samplesPerWindow)])
            elif makrer + samplesPerWindow > len(self.ecgData)
                pos = numpy.amax(self.ecgData[(marker-samplesPerWindow):])

            if len(output) > RR_SIZE:
                self.rr = numpy.diff(output[-(RR_SIZE+1):])
                self.rrMean = numpy.mean(self.rr)
                if (self.rr[-1] <= RR_LOW_LIMIT_FACTOR * self.rrMean) || (self.rr[-1] >= RR_HIGH_LIMIT_FACTOR * self.rrMean):
                    self.thresholdf1 *= 0.5;
                    self.thresholdf2 *= 0.5;

            if self.rrMean > 0.0:


    # def findMaximumAroundPoint(self, signal, point, window):
    #     pos = 0
    #     if marker > samplesPerWindow && (marker + samplesPerWindow) < len(self.ecgData):
    #         pos = numpy.amax(self.ecgData[(marker-samplesPerWindow):(marker+samplesPerWindow)])
    #     elif i == 0:
    #         pos = numpy.amax(self.ecgData[:(marker+samplesPerWindow)])
    #     elif makrer + samplesPerWindow > len(self.ecgData)
    #         pos = numpy.amax(self.ecgData[(marker-samplesPerWindow):])
    #     return pos

    def initializeThresholds(self):
        self.spki = numpy.amax(self.integratedSignal[:2*self.samplingFrequency]) / 3.0
        self.npki = numpy.sum(self.integratedSignal[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdi1 = self.spki
        self.thresholdi2 = self.npki
        self.spkf = numpy.amax(self.ecgData[:2*self.samplingFrequency]) / 3.0
        self.npkf = numpy.sum(self.ecgData[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdf1 = self.spkf
        self.thresholdf2 = self.npkf

    def searchBack(self):
        pass

    def findFiducialMarkers(self):
        gradient = numpy.gradient(self.integratedSignal)
        fiducialMarkers = []
        for i in range(2, len(gradient)):
            if gradient[i-1] > 0 && gradient[i] <= 0:
                if fiducialMarkers:
                    if fiducialMarkers[-1] >= math.ceil(0.2 * self.samplingFrequency):
                        fiducialMarkers.append(i)
                    elif self.integratedSignal[i] >= self.integratedSignal[fiducialMarkers[-1]]:
                        fiducialMarkers[-1] = i
                else:
                    fiducialMarkers.append(i)

        return fiducialMarkers

