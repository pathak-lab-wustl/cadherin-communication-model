

import time
import csv
import os

from cc3d import CompuCellSetup

start = time.time()
from Spheroid_invasionSteppables import CalculateARandFrontLeaders
CompuCellSetup.register_steppable(steppable=CalculateARandFrontLeaders(frequency=1))

from Spheroid_invasionSteppables import Calculate_P
CompuCellSetup.register_steppable(steppable=Calculate_P(frequency=1))

from Spheroid_invasionSteppables import Spheroid_invasionSteppable
CompuCellSetup.register_steppable(steppable=Spheroid_invasionSteppable(frequency=1))

CompuCellSetup.run()

end = time.time()
print(end - start, 'this is the time for simulation')
# print(dir(CompuCellSetup.persistent_globals._PersistentGlobals__output_dir))
# print(CompuCellSetup.persistent_globals._PersistentGlobals__output_dir)

# outputDir = os.path.dirname(CompuCellSetup.persistent_globals._PersistentGlobals__output_dir)
outputDirPath = os.path.dirname(CompuCellSetup.persistent_globals._PersistentGlobals__output_dir)

simTimeFile = 'simTime.csv'
simTimeFilePath = os.path.join(outputDirPath, simTimeFile)

with open(simTimeFilePath, 'a') as fout:
    writer = csv.writer(fout, delimiter=',')
    writer.writerow(['simTime_sec', (end-start)])
