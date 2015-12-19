from EMD import EMD
import cProfile
import pstats
import unittest
import os
import numpy as np

class EMDTest(unittest.TestCase):

    def setUp(self):
        refDir = os.path.dirname(os.getcwd())
        refNumStr = '103'
        refInFileDir = os.path.join(refDir, refNumStr, 'Input.csv')
        refOutFileDir = os.path.join(refDir, refNumStr, 'EMDOutput.csv')
        resOutFileDir = os.path.join(refDir, refNumStr, 'EMDResultsPython.csv')
        self.imfOutFileDir = os.path.join(refDir, refNumStr, 'Imf.csv')
        envOutFileDir = os.path.join(refDir, refNumStr, 'Envelope.csv')

        self.samplingFreq = 360.0
        self.imfsNo = 3
        self.tWindowSize = 120.0

        self.refInFile = open(refInFileDir, 'r')
        self.refOutFile = open(refOutFileDir, 'r')
        self.resultFile = open(resOutFileDir, 'w')
        # self.imfOutFileFile = open(imfOutFileDir, 'w')
        self.envOutFileFile = open(envOutFileDir, 'w')
        self.refInVector = []
        self.refOutVector = []
        self.resultVector = []
        self.refIntSignalVector = []
        self.refMarksVector = []

    def tearDown(self):
        self.refInFile.close()
        self.refOutFile.close()
        self.resultFile.close()
        # self.imfOutFileFile.close()
        self.envOutFileFile.close()

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
        self.alg = EMD(self.refInVector,
                       tWindowSize=self.tWindowSize,
                       samplingFrequency=self.samplingFreq,
                       imfsNo=self.imfsNo)
        self.resultVector = self.alg.process()
        self.env = self.alg.hilbertAlgorithm.envelope
        self.imfs = self.alg.imfs
        # self.spline_upper = self.alg.spline_upper

        for val in self.resultVector or []:
            self.resultFile.write('%d\n' % val)

        # for val in self.imf:
        #     self.imfOutFileFile.write('%f\n' % val)
        np.savetxt(self.imfOutFileDir, self.imfs, delimiter=',')

        for val in self.env:
            self.envOutFileFile.write('%f\n' % val)

    def isResultIdentical(self):
        result = True

        # self.assertTrue(len(self.resultVector) > 0)

        print 'Output size: {}'.format(len(self.resultVector))

        for ref, res in zip(self.refOutVector, self.resultVector):
            # print '{} {}'.format(ref, res)
            if ref != res:
                print 'reference ({}) and result ({}) are not equal'.format(ref, res)
                result = False

        return result

    def testAlgorithm(self):
        self.loadReferenceData()
        self.runAlgorithm()

        self.assertTrue(self.isResultIdentical(), 'Reference vector and result vector are not identical')

if __name__ == '__main__':
    unittest.main()