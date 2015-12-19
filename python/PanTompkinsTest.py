from PanTompkins import PanTompkins
import numpy as np
import unittest
import os

class PanTompkinsTest(unittest.TestCase):

    def setUp(self):
        refDir = os.path.dirname(os.getcwd())
        refNumStr = '101'
        refInFileDir = os.path.join(refDir, refNumStr, 'Input.csv')
        refOutFileDir = os.path.join(refDir, refNumStr, 'PanTompkinsOutput.csv')
        resOutFileDir = os.path.join(refDir, refNumStr, 'PanTompkinsResultsPython.csv')

        self.samplingFreq = 360.0

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
        self.alg = PanTompkins(self.refInVector, samplingFrequency=self.samplingFreq)
        self.resultVector = self.alg.process()
        for val in self.resultVector:
            self.resultFile.write('%d\n' % val)

    def comapareOtherValues(self):
        # 1st line
        intSignalStr = self.refAdditionalFile.readline().split(',')
        del intSignalStr[-1]
        for item in intSignalStr:
            self.refIntSignalVector.append(float(item))

        for ref, res in zip(self.refIntSignalVector[28:], self.alg.nSignal):
            diff = ref - res
            if diff > 0.0000001:
                print diff

        # 2nd line
        marksStr = self.refAdditionalFile.readline().split(',')
        del marksStr[-1]
        for item in marksStr:
            self.refMarksVector.append(int(float(item)))

        for ref, res in zip(self.refMarksVector, self.alg.fiducialMarks):
            self.assertEqual(ref, res)
        pass

    def isResultIdentical(self):
        result = True

        for ref, res in zip(self.refOutVector, self.resultVector):
            print '{} {}'.format(ref, res)
            if ref != res:
                print 'reference ({}) and result ({}) are not equal'.format(ref, res)
                result = False

        return result

    def testAlgorithm(self):
        self.loadReferenceData()
        self.runAlgorithm()

        # self.comapareOtherValues()
        self.assertTrue(self.isResultIdentical(), 'Reference vector and result vector are not identical')

if __name__ == '__main__':
    unittest.main()
