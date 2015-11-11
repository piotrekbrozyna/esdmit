import PanTompkins

class PanTompkinsTest(object):

    def __init__(self, refInFileName, refOutFileName, resOutFileName):
        self.refInFile = open(refInFileName, 'r')
        self.refOutFile = open(refOutFileName, 'r')
        self.resultFile = open(resOutFileName, 'w')

    def run(self):
        self.loadData()
        self.startAlgorithm()
        self.compareWithTestVectors()
        self.closeFiles()

    def loadData(self):
        pass

    def startAlgorithm(self):
        pass

    def compareWithTestVectors(self):
        pass

    def closeFiles(self):
        pass