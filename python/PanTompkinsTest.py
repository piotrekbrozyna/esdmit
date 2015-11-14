from PanTompkins import PanTompkins
import numpy as np
import unittest

class PanTompkinsTest(unittest.TestCase):

    # def setUp(self, refInFileName, refOutFileName, resOutFileName):
    def setUp(self):
        refInFileName = 'PanTompkinsInput.csv'
        refOutFileName = 'PanTompkinsOutput.csv'
        resOutFileName = 'results.txt'

        self.refInFile = open(refInFileName, 'r')
        self.refOutFile = open(refOutFileName, 'r')
        self.resultFile = open(resOutFileName, 'w')
        self.refInVector = []
        self.refOutVector = []
        self.resultVector = []

    def tearDown(self):
        self.refInFile.close()
        self.refOutFile.close()
        self.resultFile.close()

    def loadReferenceData(self):
        for line in self.refInFile:
            val = float(line.rstrip('\n').strip('\''))
            self.refInVector.append(val)

        for line in self.refOutFile:
            val = int(float(line.rstrip('\n').strip('\'')))
            self.refOutVector.append(val)

    def runAlgorithm(self):
        self.alg = PanTompkins(self.refInVector, 1000)
        self.resultVector = self.alg.process()
        for val in self.resultVector:
            self.resultFile.write('%d\n' % val)

    def isResultIdentical(self):
        refLength = len(self.refOutVector)
        resLength = len(self.resultVector)

        self.assertTrue(refLength == resLength,
                        'Lengths of reference vector ({}) and result vector ({}) are not equal'.format(refLength, resLength))

        for ref, res in zip(self.refOutVector, self.resultVector):
            print '{} {}'.format(ref, res)
            self.assertTrue(ref == res,
                            'reference ({}) and result ({}) are not equal'.format(ref, res))

        return self.refOutVector == self.resultVector

    def testAlgorithm(self):
        self.loadReferenceData()
        self.runAlgorithm()

        self.assertTrue(self.isResultIdentical(), 'Reference vector and result vector are not identical')

if __name__ == '__main__':
    unittest.main()