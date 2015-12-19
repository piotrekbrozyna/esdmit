import os
from PanTompkins import PanTompkins

refDir = os.path.dirname(os.getcwd())
refNumStr = '102'
refInFileDir = os.path.join(refDir, refNumStr, 'Input.csv')
resOutFileDir = os.path.join(refDir, refNumStr, 'PanTompkinsResultsPython.csv')

refInFile = open(refInFileDir, 'r')
resultFile = open(resOutFileDir, 'w')

refInVector = []

refInVectorStr = refInFile.readline().split(',')
del refInVectorStr[-1]
for item in refInVectorStr:
    val = float(item)
    refInVector.append(val)

alg = PanTompkins(refInVector,
                  samplingFrequency=360)
resultVector = alg.process()

for val in resultVector or []:
    resultFile.write('%d\n' % val)

refInFile.close()
resultFile.close()