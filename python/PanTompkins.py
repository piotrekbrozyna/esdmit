from math import floor, ceil
from numpy import amax, convolve, square, diff, gradient, mean, multiply
import Tools

class PanTompkins(object):

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


    def __init__(self, ecgSignal, samplingFrequency=1000, tIntegrationWindow=150, tThRadius=0.1):
        self.ecgSignal = ecgSignal
        self.ecgSignalSize = len(ecgSignal)
        self.samplingFrequency = samplingFrequency
        self.tIntegrationWindow = tIntegrationWindow
        self.nThRadius = samplingFrequency * tThRadius
        self.fiducialMarks = []
        self.rPeaks = []
        self.rr = 0.0
        self.rrMean = 0.0

    def process(self):
        dSignal = self.differentiateEcgSignal()
        sSignal = square(dSignal)
        self.iSignal = Tools.integrateInMovingWindow(sSignal, self.samplingFrequency, self.tIntegrationWindow)
        self.detectRPeaks()

    def differentiateEcgSignal(self):
        return convolve(self.ecgSignal, PanTompkins.DIFF_KERNEL, mode='full')

    def detectRPeaks(self):
        Tools.zeroIfUnderThreshold(self.iSignal, PanTompkins.SIGNAL_LOW_LIMIT)
        self.findfiducialMarks()
        self.initializeThresholds()
        self.threshold()

    def threshold(self):
        skip = False

        for i, marker in enumerate(self.fiducialMarks):
            pos = Tools.findMaximumWithinNeigborhood(self.ecgSignal, marker, self.nThRadius)

            if len(self.rPeaks) > PanTompkins.RR_SIZE:
                self.rr = diff(self.rPeaks[-(PanTompkins.RR_SIZE+1):])
                self.rrMean = mean(self.rr)
                if (self.rr[-1] <= PanTompkins.RR_LOW_LIMIT_FACTOR * self.rrMean) or (self.rr[-1] >= PanTompkins.RR_HIGH_LIMIT_FACTOR * self.rrMean):
                    self.thresholdf1 *= 0.5;
                    self.thresholdf2 *= 0.5;

            if (self.rrMean > 0.0) and (marker - self.rPeaks[-1]) > (PanTompkins.RR_MISSED_LIMIT_FACTOR * self.rrMean):
                # self.searchBack()
                windowSamples = floor(0.2*self.samplingFrequency)
                fiducialMarkTmp = amax(self.iSignal[(self.rPeaks[-1]+windowSamples):(marker-windowSamples)])
                posTmp = Tools.findMaximumWithinNeigborhood(self.ecgSignal, fiducialMarkTmp, self.nThRadius)

                if self.iSignal[fiducialMarkTmp] > self.thresholdi2:
                    if self.ecgSignal[posTmp] > self.thresholdf2:
                        self.rPeaks.append(posTmp)
                        self.spkf = PanTompkins.SECOND_THRESHOLD_SPKF_PEAKF_FACTOR * self.ecgSignal[posTmp] + PanTompkins.SECOND_THRESHOLD_SPKF_SPKF_FACTOR * self.spkf
                    self.spki = PanTompkins.SECOND_THRESHOLD_SPKI_PEAKI_FACTOR * self.iSignal[fiducialMarkTmp] + PanTompkins.SECOND_THRESHOLD_SPKI_SPKI_FACTOR * self.spki

            if (self.iSignal[marker] >= self.thresholdi1):
                if len(self.rPeaks) > 2:
                    # T wave?
                    if (marker - self.rPeaks[-1]) <= ceil(PanTompkins.T_WAVE_DETECTION_RR_HIGH_LIMIT * self.samplingFrequency):
                        # slopes
                        vector1 = self.iSignal[(marker-4):(marker-3)] + self.iSignal[(marker-1):marker]
                        vector2 = self.rPeaks[-5:-4] + self.rPeaks[-2:-1]
                        kernel = [-2, -1, 1, 2]
                        slope1 = multiply(vector1, kernel)
                        slope2 = multiply(vector2, kernel)
                        if (abs(slope1) <= 0.5 * abs(slope2)):
                            self.npki = PanTompkins.NPKI_UPDATE_FACTOR_1 * self.iSignal[marker] + PanTompkins.NPKI_UPDATE_FACTOR_2 * self.npki
                            self.npkf = PanTompkins.NPKF_UPDATE_FACTOR_1 * self.ecgSignal[pos] + PanTompkins.NPKF_UPDATE_FACTOR_2 * self.npkf
                            skip = True
                        else:
                            skip = False
                if not skip:
                    if self.ecgSignal[pos] > self.thresholdf1:
                        self.rPeaks.append(pos)
                        self.spkf = PanTompkins.SPKF_UPDATE_FACTOR_1 * self.ecgSignal[pos] + PanTompkins.SPKF_UPDATE_FACTOR_2 * self.spkf
                    self.spki = PanTompkins.SPKI_UPDATE_FACTOR_1 * self.iSignal[marker] + PanTompkins.SPKI_UPDATE_FACTOR_2 * self.spki

            elif (self.thresholdf2 <= self.iSignal[marker] and self.iSignal[marker] < self.thresholdi1) or (self.iSignal[marker] < self.thresholdi2):
                self.npki = PanTompkins.NPKI_UPDATE_FACTOR_1 * self.iSignal[marker] + PanTompkins.NPKI_UPDATE_FACTOR_2 * self.npki
                self.npkf = PanTompkins.NPKF_UPDATE_FACTOR_1 * self.ecgSignal[pos] + PanTompkins.NPKF_UPDATE_FACTOR_2 * self.npkf

            self.updateThresholds()

            skip = False

    # def possibleTWave(self, )



    def initializeThresholds(self):
        self.spki = amax(self.iSignal[:2*self.samplingFrequency]) / 3.0
        self.npki = sum(self.iSignal[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdi1 = self.spki
        self.thresholdi2 = self.npki
        self.spkf = amax(self.ecgSignal[:2*self.samplingFrequency]) / 3.0
        self.npkf = sum(self.ecgSignal[:2*self.samplingFrequency]) / (2*self.samplingFrequency*2.0)
        self.thresholdf1 = self.spkf
        self.thresholdf2 = self.npkf

    def updateThresholds(self):
        self.thresholdi1 = self.npki + PanTompkins.THRESHOLDI_UPDATE_FACTOR_1 * (self.spki - self.npki)
        self.thresholdi2 = PanTompkins.THRESHOLDI_UPDATE_FACTOR_2 * self.thresholdi1
        self.thresholdf1 = self.npkf + PanTompkins.THRESHOLDF_UPDATE_FACTOR_1 * (self.spkf - self.npkf)
        self.thresholdf2 = PanTompkins.THRESHOLDF_UPDATE_FACTOR_2 * self.thresholdf1

    def updateApproximations(self):
        pass

    def searchBack(self, peaks):
        pass

    def findfiducialMarks(self):
        gSignal = gradient(self.iSignal)
        for i in range(2, len(gSignal)):
            if gSignal[i-1] > 0 and gSignal[i] <= 0:
                if self.fiducialMarks:
                    if self.fiducialMarks[-1] >= ceil(0.2 * self.samplingFrequency):
                        self.fiducialMarks.append(i)
                    elif self.iSignal[i] >= self.iSignal[self.fiducialMarks[-1]]:
                        self.fiducialMarks[-1] = i
                else:
                    self.fiducialMarks.append(i)

        return self.fiducialMarks

