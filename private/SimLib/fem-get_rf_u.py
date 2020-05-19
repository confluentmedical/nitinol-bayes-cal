#
# Extract displacement at TENSION-CONTROL-PT of the last step
#
# "A Probabilistic Approach with Built-in Uncertainty Quantification for 
# the Calibration of a Superelastic Constitutive Model from Full-field 
# Strain Data"
#
# https://github.com/confluentmedical/nitinol-bayes-cal
#
#    Copyright 2020 Harshad M. Paranjape
# 
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
from odbAccess import *
import sys
# Get frame numbers from cmd arguments
odbname = sys.argv[8]
outname  = 'rf_u_' + os.path.splitext(odbname)[0] + '.data'

f = open(outname, 'w')
# Open odb
odb = openOdb(path=odbname)
# Last step
step = odb.steps[odb.steps.keys()[-1]]
# Loop over frames
for frame in step.frames:
    # Get RF2 and U2 at TENSION-CONTROL-PT in this frame
    rf2 = frame.fieldOutputs['RF'].getSubset(position=NODAL, region=odb.rootAssembly.nodeSets['TENSION-CONTROL-PT']).getScalarField(componentLabel='RF2').values
    u2  = frame.fieldOutputs['U'].getSubset(position=NODAL, region=odb.rootAssembly.nodeSets['TENSION-CONTROL-PT']).getScalarField(componentLabel='U2').values
    # Print RF2, U2 to file
    f.write('%14.4e %14.4e' %
            (float(u2[0].data), float(rf2[0].data)))
    f.write('\n')
    f.flush()
# Done
f.close()
odb.close()

