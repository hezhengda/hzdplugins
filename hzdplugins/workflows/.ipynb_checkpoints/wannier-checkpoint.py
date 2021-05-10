# This script can help us calculate the Wannier function
import os
import sys

from hzdplugins.aiidaplugins.utilities import pwInputToDict, dictToPwInput
from hzdplugins.aiidaplugins.info import getOptimizedStructure, getXMLFromPW

# prepare the input files
dict_input = pwInputToDict(loc_file=sys.argv[1])

# get optimized structure from xml files
data_xml = getXMLFromPW(xml_file=sys.argv[2])
newCell, newAtomicPositions = getOptimizedStructure(data_xml)

dict_input['CELL_PARAMETERS'] = newCell
dict_input['ATOMIC_POSITIONS'] = newAtomicPositions

# prepare the scf input file: inp_scf

# do the scf simulation


# do the pmw.x calculation

# use the wannier function as projector

# calculate the pdos
