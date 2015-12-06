from EMD import EMD
import numpy as np
import unittest
import os

class EMDTest(unittest.TestCase):

    def setUp(self):
        refDir = os.path.dirname(os.getcwd())
        refNumStr = '101'
        refInFileDir = os.path.join(refDir, refNumStr, 'PanTompkinsInput.csv')
        refOutFileDir = os.path.join(refDir, refNumStr, 'PanTompkinsOutput.csv')
        resOutFileDir = os.path.join(refDir, refNumStr, 'EMDResultsPython.csv')

        self.samplingFreq = 360.0
        self.imfsNo = 1

        self.refInFile = open(refInFileDir, 'r')
        self.refOutFile = open(refOutFileDir, 'r')
        # self.refAdditionalFile = open(refAdditionalFileDir, 'r')
        self.resultFile = open(resOutFileDir, 'w')
        self.refInVector = []
        self.refOutVector = []
        self.resultVector = []
        self.refIntSignalVector = []
        self.refMarksVector = []

    def tearDown(self):
        self.refInFile.close()
        self.refOutFile.close()
        # self.refAdditionalFile.close()
        self.resultFile.close()

    def loadReferenceData(self):
        refInVectorStr = self.refInFile.readline().split(',')
        del refInVectorStr[-1]
        for item in refInVectorStr:
            val = float(item)
            self.refInVector.append(val)

        refOutVectorStr = self.refOutFile.readline().split(',')
        del refOutVectorStr[-1]
        for item in refOutVectorStr:
            val = int(float(item))
            self.refOutVector.append(val)

    def runAlgorithm(self):
        self.alg = EMD(self.refInVector, samplingFrequency=self.samplingFreq, imfsNo=self.imfsNo)
        self.resultVector = self.alg.process()
        for val in self.resultVector or []:
            self.resultFile.write('%d\n' % val)

    def isResultIdentical(self):
        refLength = len(self.refOutVector)
        resLength = len(self.resultVector)

        # self.assertTrue(refLength == resLength,
        #                 'Lengths of reference vector ({}) and result vector ({}) are not equal'.format(refLength, resLength))

        for ref, res in zip(self.refOutVector, self.resultVector):
            print '{} {}'.format(ref, res)
            # self.assertTrue(ref == res,
            #                 'reference ({}) and result ({}) are not equal'.format(ref, res))

    def testAlgorithm(self):
        self.loadReferenceData()
        self.runAlgorithm()

        # self.comapareOtherValues()
        self.assertTrue(self.isResultIdentical(), 'Reference vector and result vector are not identical')

if __name__ == '__main__':
    unittest.main()