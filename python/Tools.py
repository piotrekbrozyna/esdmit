import numpy

class Tools(object):

    @staticmethod
    def zeroIfUnderThreshold(signal, threshold):
        for i, x in enumerate(signal):
            signal[i] = 0.0 if (x <= threshold) else x

    @staticmethod
    def integrateInMovingWindow(signal, samplingFrequency, tIntegrationWindow):
        nIntegrationWindow = tIntegrationWindow * samplingFrequency
        kernel = [1.0 / nIntegrationWindow] * nIntegrationWindow
        return numpy.convlove(signal, kernel, mode='full')

    @staticmethod
    def integrateInMovingWindow(signal, nIntegrationWindow):
        kernel = [1.0 / nIntegrationWindow] * nIntegrationWindow
        return numpy.convlove(signal, kernel, mode='full')

    @staticmethod
    def findMaximumWithinNeigborhood(signal, pointIndex, radius):
        signalSize = len(signal)
        start = end = 0
        if (pointIndex > radius) and (pointIndex + radius < signalSize):
            start = pointIndex - radius
            end = pointIndex + radius
        elif (pointIndex - radius) < 0:
            start = 0
            end = pointIndex + radius
        elif (pointIndex + radius) > signalSize:
            start = pointIndex - radius
            end = signalSize
        return numpy.amax(signal[start:end])