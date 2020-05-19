#
# Extract strains and coordinates in the GAGE element set at the first and last frame of the last step
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
outname_first  = 's_e_maxima_' + os.path.splitext(odbname)[0] + '.data'
outname_last   = 's_e_minima_' + os.path.splitext(odbname)[0] + '.data'
f_first = open(outname_first, 'w')
f_last = open(outname_last, 'w')
# Open odb
odb = openOdb(path=odbname)
step = odb.steps[odb.steps.keys()[-1]]
# Get strain values at the first frame of last step
le11     = step.frames[0].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE11').values
le22     = step.frames[0].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE22').values
le12     = step.frames[0].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE12').values
le_princ = step.frames[0].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(invariant=MAX_PRINCIPAL).values
coord1   = step.frames[0].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[0]).values
coord2   = step.frames[0].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[1]).values
coord3   = step.frames[0].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[2]).values
# Loop over all integration points
for jj in range(len(le11)):
    f_first.write('%14.4e %14.4e %14.4e %14.4e %14.4e %14.4e' %
            (float(coord1[jj].data), float(coord2[jj].data), float(le11[jj].data), float(le22[jj].data), float(le12[jj].data), float(le_princ[jj].data)))
    f_first.write('\n')
    f_first.flush()
# Done
f_first.close()
# Get strain values at the last frame of last step
le11     = step.frames[-1].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE11').values
le22     = step.frames[-1].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE22').values
le12     = step.frames[-1].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel='LE12').values
le_princ = step.frames[-1].fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(invariant=MAX_PRINCIPAL).values
coord1   = step.frames[-1].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[0]).values
coord2   = step.frames[-1].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[1]).values
coord3   = step.frames[-1].fieldOutputs['COORD'].getSubset(position=INTEGRATION_POINT, region=odb.rootAssembly.instances['NDC-56-03018_REV3-1'].elementSets['GAGE-FACE']).getScalarField(componentLabel=step.frames[-1].fieldOutputs['COORD'].componentLabels[2]).values
# Loop over all integration points
for jj in range(len(le11)):
    f_last.write('%14.4e %14.4e %14.4e %14.4e %14.4e %14.4e' %
            (float(coord1[jj].data), float(coord2[jj].data), float(le11[jj].data), float(le22[jj].data), float(le12[jj].data), float(le_princ[jj].data)))
    f_last.write('\n')
    f_last.flush()
# Done
f_first.close()

odb.close()

