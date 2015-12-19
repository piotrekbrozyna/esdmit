import os
from EMD import EMD

refDir = os.path.dirname(os.getcwd())
refNumStr = '101'
refInFileDir = os.path.join(refDir, refNumStr, 'Input.csv')
resOutFileDir = os.path.join(refDir, refNumStr, 'EMDResultsPython.csv')

refInFile = open(refInFileDir, 'r')
resultFile = open(resOutFileDir, 'w')

refInVector = []

refInVectorStr = refInFile.readline().split(',')
del refInVectorStr[-1]
for item in refInVectorStr:
    val = float(item)
    refInVector.append(val)

alg = EMD(refInVector,
          tWindowSize=120.0,
          samplingFrequency=360,
          imfsNo=3)
resultVector = alg.process()

for val in resultVector or []:
    resultFile.write('%d\n' % val)

refInFile.close()
resultFile.close()