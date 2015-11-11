from math import floor, ceil
from numpy import amax, convolve, square, diff, gradient, mean, multiply

class PanTompkins(Class):

    RR_LOW_LIMIT_FACTOR = 0.92
    RR_HIGH_LIMIT_FACTOR = 1.16
    RR_MISSED_LIMIT_FACTOR = 1.66
    RR_SIZE = 8
    DIFF_KERNEL = [-0.125, -0.250, 0, 0.250, 0.125]
    SIGNAL_LOW_LIMIT = 0.0001

    # 360 ms
    T_WAVE_DETECTION_RR_HIGH_LIMIT = 0.36

    # thresholding factors
    SECOND_THRESHOLD_SPKF_PEAKF_FACTOR = 0.25
    SECOND_THRESHOLD_SPKF_SPKF_FACTOR = 0.75
    SECOND_THRESHOLD_SPKI_PEAKI_FACTOR = 0.25
    SECOND_THRESHOLD_SPKI_SPKI_FACTOR = 0.75

    SPKI_UPDATE_FACTOR_1 = 0.125
    SPKI_UPDATE_FACTOR_2 = 0.875
    NPKI_UPDATE_FACTOR_1 = 0.125
    NPKI_UPDATE_FACTOR_2 = 0.875
    SPKF_UPDATE_FACTOR_1 = 0.125
    SPKF_UPDATE_FACTOR_2 = 0.875
    NPKF_UPDATE_FACTOR_1 = 0.125
    NPKF_UPDATE_FACTOR_2 = 0.875
    THRESHOLDI_UPDATE_FACTOR_1 = 0.25
    THRESHOLDI_UPDATE_FACTOR_2 = 0.5
    THRESHOLDF_UPDATE_FACTOR_1 = 0.25
    THRESHOLDF_UPDATE_FACTOR_2 = 0.5


    def __init__(self, ecgData, samplingFrequency, numOfSamplesPerWindow):
        self.ecgData = ecgData
        self.samplingFrequency = samplingFrequency
        self.integrationWindowTime = float(numOfSamplesPerWindow) / samplingFrequency
        self.intervalCount = math.floor(samplingFrequency * integrationWindowTime)

    def process(self):
        differentiatedSignal = self.differentiate(self.ecgData)
        squaredSignal = square(differentiatedSignal)
        self.integratedSignal = self.integrateInWindow(squaredSignal)
        self.rPeaks = detectRPeaks(ecgData, self.integratedSignal, gradientOfIntegratedSignal)

    def integrateInWindow(self, signal):
        kernel = [1.0 / self.intervalCount] * self.intervalCount
        return convlove(signal, kernel, mode='full')

    def differentiate(self, signal):
        return convolve(ecgData, DIFF_KERNEL, mode='full')

    def detectRPeaks(self):
        for i, x in enumerate(self.integratedSignal):
            signal[i] = (x[i] <= SIGNAL_LOW_LIMIT) ? 0.0 : x[i]

        self.fiducialMarkers = self.findFiducialMarkers()
        self.initializeThresholds()
        self.threshold()

    def threshold(self):
        windowTime = 0.1
        samplesPerWindow = self.samplingFrequency * windowTime
        self.rPeaks = []
        self.rr = 0.0
        self.rrMean = 0.0
        skip = False
        ecgDataLen = len(self.ecgData)

        for i, marker in enumerate(fiducialMarkers):
            if marker > samplesPerWindow && (marker + samplesPerWindow) < ecgDataLen:
                pos = amax(self.ecgData[(marker-samplesPerWindow):(marker+samplesPerWindow)])
            elif i == 0:
                pos = amax(self.ecgData[:(marker+samplesPerWindow)])
            elif makrer + samplesPerWindow > len(self.ecgData)
                pos = amax(self.ecgData[(marker-samplesPerWindow):])

            if len(self.rPeaks) > RR_SIZE:
                self.rr = diff(self.rPeaks[-(RR_SIZE+1):])
                self.rrMean = mean(self.rr)
                if (self.rr[-1] <= RR_LOW_LIMIT_FACTOR * self.rrMean) || (self.rr[-1] >= RR_HIGH_LIMIT_FACTOR * self.rrMean):
                    self.thresholdf1 *= 0.5;
                    self.thresholdf2 *= 0.5;

            if (self.rrMean > 0.0) && (marker - self.rPeaks[-1]) > (RR_MISSED_LIMIT_FACTOR * self.rrMean):
                # self.searchBack()
                windowSamples = floor(0.2*self.samplingFrequency)
                fiducialMarkerTmp = amax(self.integratedSignal[(self.rPeaks[-1]+windowSamples):(marker-windowSamples)]) 
                # fiducialMarkerTmp -= self.rPeaks[-1] + windowSamples
                # fiducialMarkerTmp += self.rPeaks[-1] + windowSamples
                if (fiducialMarkerTmp > samplesPerWindow) && (fiducialMarkerTmp + samplesPerWindow < ecgDataLen):
                    posTmp = amax(self.ecgData[():(fiducialMarkerTmp + samplesPerWindow)])
                elif i == 0:
                    posTmp = amax(self.ecgData[:(fiducialMarkerTmp + samplesPerWindow)])
                elif (fiducialMarkerTmp + samplesPerWindow) > ecgDataLen:
                    posTmp = ....

                if self.integratedSignal[fiducialMarkerTmp] > self.thresholdi2:
                    if self.ecgData[posTmp] > self.thresholdf2:
                        self.rPeaks.append(posTmp)
                        self.spkf = SECOND_THRESHOLD_SPKF_PEAKF_FACTOR * self.ecgData[posTmp] + SECOND_THRESHOLD_SPKF_SPKF_FACTOR * self.spkf
                    self.spki = SECOND_THRESHOLD_SPKI_PEAKI_FACTOR * self.integratedSignal[fiducialMarkerTmp] + SECOND_THRESHOLD_SPKI_SPKI_FACTOR * self.spki

            if (self.integratedSignal[marker] >= self.thresholdi1):
                if len(self.rPeaks) > 2:
                    # T wave?
                    if (marker - self.rPeaks[-1]) <= ceil(T_WAVE_DETECTION_RR_HIGH_LIMIT * self.samplingFrequency):
                        # slopes
                        vector1 = self.integratedSignal[(marker-4):(marker-3)] + self.integratedSignal[(marker-1):marker]
                        vector2 = self.rPeaks[-5:-4] + self.rPeaks[-2:-1]
                        kernel = [-2, -1, 1, 2]
                        slope1 = multiply(vector1, kernel)
                        slope2 = multiply(vector2, kernel)
                        if (abs(slope1) <= 0.5 * abs(slope2)):
                            self.npki = NPKI_UPDATE_FACTOR_1 * self.integratedSignal[marker] + NPKI_UPDATE_FACTOR_2 * self.npki
                            self.npkf = NPKF_UPDATE_FACTOR_1 * self.ecgData[pos] + NPKF_UPDATE_FACTOR_2 * self.npkf
                            skip = True
                        else:
                            skip = False
                if not skip:
                    if self.ecgData[pos] > self.thresholdf1:
                        self.rPeaks.append(pos)
                        self.spkf = SPKF_UPDATE_FACTOR_1 * self.ecgData[pos] + SPKF_UPDATE_FACTOR_2 * self.spkf
                    self.spki = SPKI_UPDATE_FACTOR_1 * self.integratedSignal[marker] + SPKI_UPDATE_FACTOR_2 * self.spki

            elif (self.thresholdf2 <= self.integratedSignal[marker] && self.integratedSignal[marker] < self.thresholdi1) || self.integratedSignal[marker] < self.thresholdi2:
                self.npki = NPKI_UPDATE_FACTOR_1 * self.integratedSignal[marker] + NPKI_UPDATE_FACTOR_2 * self.npki
                self.npkf = NPKF_UPDATE_FACTOR_1 * self.ecgData[pos] + NPKF_UPDATE_FACTOR_2 * self.npkf

            self.updateThresholds()

            skip = False

    def possibleTWave(self, )

    def findMaximumWithinInterval(self, signal, startIndex, endIndex):
        pass

    def initializeThresholds(self):
        self.spki = amax(self.integratedSignal[:2*self.samplingFrequency]) / 3.0
        self.npki = sum(self.integratedSignal[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdi1 = self.spki
        self.thresholdi2 = self.npki
        self.spkf = amax(self.ecgData[:2*self.samplingFrequency]) / 3.0
        self.npkf = sum(self.ecgData[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdf1 = self.spkf
        self.thresholdf2 = self.npkf

    def updateThresholds(self):
        self.thresholdi1 = self.npki + THRESHOLDI_UPDATE_FACTOR_1 * (self.spki - self.npki)
        self.thresholdi2 = THRESHOLDI_UPDATE_FACTOR_2 * self.thresholdi1
        self.thresholdf1 = self.npkf + THRESHOLDF_UPDATE_FACTOR_1 * (self.spkf - self.npkf)
        self.thresholdf2 = THRESHOLDF_UPDATE_FACTOR_2 * self.thresholdf1

    def updateApproximations(self):
        pass

    def searchBack(self, peaks):
        pass

    def findFiducialMarkers(self):
        gradient = gradient(self.integratedSignal)
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

